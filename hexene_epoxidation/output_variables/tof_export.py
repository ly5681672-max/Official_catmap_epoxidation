
"""CatMAP: 14-site TOF, coverage, DRC, reaction flux and evidence audit.

Drop-in replacement for hexene_epoxidation/output_variables/tof_export.py.
Called by existing mkm_job.py: export_tof_tables(model, output_dir, product='C6H12O_g').
Does not change the reaction model, parameters, input energies, or conditions.
"""
from pathlib import Path
import csv
from mpmath import mp as _catmap_shared_mp

# Independent numeric context: never change CatMAP solver's global mp.dps.
mp = _catmap_shared_mp.clone()
mp.dps = 160

TS1_NAME = 'OOH-C6H12_s'
OOH_NAME = 'OOH_s'
ALKENE_NAME = 'C6H12_s'


def _number(x, what):
    try:
        value = mp.mpf(str(x))
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{what}: invalid number {x!r}') from exc
    if not mp.isfinite(value):
        raise ValueError(f'{what}: nonfinite number {x!r}')
    return value


def _text(x):
    return str(x)


def _save_csv(path, header, rows):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as stream:
        out = csv.writer(stream)
        out.writerow(header)
        out.writerows(rows)


def _labels_index(labels, name, context):
    values = list(labels)
    if values.count(name) != 1:
        raise ValueError(f'{context}: {name!r} not uniquely present in {values!r}')
    return values.index(name)


def _rel_mismatch(a, b):
    # Scale relative to both fluxes; 0/0 is treated as matched, not as proof of activity.
    denominator = max(abs(a), abs(b))
    return abs(a-b) / denominator if denominator else mp.mpf('0')


def _fit_xy(pairs):
    # pairs: [(theta_as_str, TOF_as_str, surface), ...]
    if len(pairs) < 3:
        return None
    with mp.workdps(80):
        xy = [(_number(x, 'theta_OOH'), mp.log10(_number(rate, 'TOF')))
              for x, rate, _ in pairs]
        n = len(xy)
        mx = sum(x for x, _ in xy) / n
        my = sum(y for _, y in xy) / n
        sxx = sum((x-mx)**2 for x, _ in xy)
        syy = sum((y-my)**2 for _, y in xy)
        if sxx == 0 or syy == 0:
            return None
        a = sum((x-mx)*(y-my) for x, y in xy) / sxx
        b = my-a*mx
        sse = sum((y-a*x-b)**2 for x, y in xy)
        loo_sse = mp.mpf('0')
        for index in range(n):
            train = xy[:index] + xy[index+1:]
            tx = sum(x for x, _ in train)/(n-1)
            ty = sum(y for _, y in train)/(n-1)
            txx = sum((x-tx)**2 for x, _ in train)
            if txx == 0:
                return dict(n=n, slope=_text(a), intercept=_text(b), R2=_text(1-sse/syy),
                            RMSE_log10=_text(mp.sqrt(sse/n)), LOOCV_Q2='UNDEFINED')
            ta = sum((x-tx)*(y-ty) for x, y in train)/txx
            tb = ty-ta*tx
            loo_sse += (xy[index][1] - ta*xy[index][0] - tb)**2
        return dict(n=n, slope=_text(a), intercept=_text(b), R2=_text(1-sse/syy),
                    RMSE_log10=_text(mp.sqrt(sse/n)), LOOCV_Q2=_text(1-loo_sse/syy))


def _fit_plot(pairs, fit, path):
    if fit is None:
        return
    import matplotlib.pyplot as plt
    xx = [float(_number(x, 'coverage')) for x, _, _ in pairs]
    yy = [float(mp.log10(_number(rate, 'TOF'))) for _, rate, _ in pairs]
    a = float(fit['slope'])
    b = float(fit['intercept'])
    fig, ax = plt.subplots(figsize=(8, 5.4))
    try:
        ax.scatter(xx, yy, s=42)
        order = sorted(range(len(xx)), key=xx.__getitem__)
        ax.plot([xx[i] for i in order], [a*xx[i]+b for i in order], linewidth=1.2)
        for x, y, (_, _, surface) in zip(xx, yy, pairs):
            ax.annotate(surface, (x, y), xytext=(4, 3),
                        textcoords='offset points', fontsize=8)
        ax.set_xlabel('OOH* coverage')
        ax.set_ylabel('log10(TOF C6H12O / s^-1)')
        ax.set_title('All-catalyst fit (not a kinetic validation)')
        ax.text(.02, .02, f"n={fit['n']}, R2={float(fit['R2']):.4f}, "
                f"LOOCV Q2={float(fit['LOOCV_Q2']):.3g}" if fit['LOOCV_Q2'] != 'UNDEFINED'
                else f"n={fit['n']}, R2={float(fit['R2']):.4f}, LOOCV Q2 undefined",
                transform=ax.transAxes, fontsize=9, va='bottom')
        fig.tight_layout()
        fig.savefig(path)
    finally:
        plt.close(fig)


def _classify(theta, x_ooh, x_ts1, r4, tof, site_ok, flux_ok, drc_ok):
    """Multiple provisional evidence tags, never mutually exclusive labels."""
    tags = []
    if theta >= mp.mpf('0.9'):
        tags.append('OOH_ACCUMULATION')
    elif theta <= mp.mpf('0.1'):
        tags.append('OOH_LOW_COVERAGE')
    else:
        tags.append('OOH_INTERMEDIATE_COVERAGE')
    if r4 > 0 and tof > 0 and flux_ok:
        tags.append('OXYGEN_TRANSFER_PRECURSOR_PATH')
    if drc_ok and abs(x_ooh) >= mp.mpf('0.1'):
        tags.append('OOH_ENERGY_SENSITIVE_CANDIDATE')
    if drc_ok and x_ts1 >= mp.mpf('0.5'):
        tags.append('TS1_RATE_CONTROL_CANDIDATE')
    if not site_ok:
        tags.append('SITE_COVERAGE_CHECK_FAILED')
    if not flux_ok:
        tags.append('FLUX_CLOSURE_FAILED')
    if not drc_ok:
        tags.append('DRC_NUMERICAL_CHECK_REQUIRED')
    tags.append('STEADY_STATE_REACHABILITY_UNVERIFIED')
    return ';'.join(tags)


def _interpret_tags(tags):
    present = set(tags.split(';'))
    messages = []
    if 'OOH_ACCUMULATION' in present:
        messages.append('OOH*高度覆盖：表面累积明显')
    elif 'OOH_LOW_COVERAGE' in present:
        messages.append('OOH*低覆盖：不属于主要积累物种')
    else:
        messages.append('OOH*处于中等覆盖区间')
    if 'OXYGEN_TRANSFER_PRECURSOR_PATH' in present:
        messages.append('氧转移反应前驱体路径：净通量与产物生成一致')
    if 'OOH_ENERGY_SENSITIVE_CANDIDATE' in present:
        messages.append('OOH*能量对速率有敏感性：待进一步数值复核')
    else:
        messages.append('尚无可靠证据说明OOH*能量显著控制速率')
    if 'TS1_RATE_CONTROL_CANDIDATE' in present:
        messages.append('氧转移TS1为候选速率控制态')
    if any(t in present for t in ('SITE_COVERAGE_CHECK_FAILED',
                                  'FLUX_CLOSURE_FAILED', 'DRC_NUMERICAL_CHECK_REQUIRED')):
        messages.append('数值核验未通过，以上均不可视为可靠机理结论')
    messages.append('稳态物理可达性尚未验证')
    return '；'.join(messages)


def export_tof_tables(model, output_dir, product='C6H12O_g'):
    """Export per-catalyst diagnostic CSV without changing CatMAP physics.

    Thresholds are descriptive screening flags, NOT universal physical criteria:
      coverage >= 0.9; |X_OOH| >= 0.1; X_TS1 >= 0.5;
      8-step net flux relative mismatch <= 1e-3; site tolerance 1e-8.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    descriptors = list(model.descriptor_names)
    surfaces = list(model.surface_names)
    if len(descriptors) != 2:
        raise ValueError('Exactly 2 descriptors required for legacy CSVs')
    if not surfaces or len(set(surfaces)) != len(surfaces):
        raise ValueError('surface_names empty or duplicated')
    # CatMAP site_names also includes the gas reservoir 'g'.
    gas_sites = set(getattr(model, '_gas_sites', ('g',)))
    ads_sites = set(model.site_names) - gas_sites
    if ads_sites != {'s'}:
        raise ValueError(f'Site balance logic expects one surface site s; found {ads_sites}')
    site_total = _number(model.species_definitions['s']['total'], 'site total')
    if site_total != 1:
        raise ValueError('Only a normalized single-site model with total=1 is supported')
    rxns = model.elementary_rxns
    if len(rxns) != 8 or not (OOH_NAME in rxns[3][0] and
                              ALKENE_NAME in rxns[3][0] and
                              TS1_NAME in rxns[3][1] and
                              'C6H12O_s' in rxns[3][-1] and
                              product in rxns[4][-1]):
        raise ValueError('Current reaction 4/5 does not match expected oxygen-transfer / desorption mechanism')

    gas_labels = list(model.output_labels.get('production_rate', []))
    prod_idx = _labels_index(gas_labels, product, 'production_rate')
    grid = getattr(model, 'production_rate_map', None)
    if not grid:
        raise ValueError('CatMAP production_rate_map is missing or empty')
    grid_rows = []
    for pos, rates in grid:
        grid_rows.append([*map(_text, pos), _text(rates[prod_idx])])

    all_rows = []
    legacy_tof_rows, legacy_pairs, pairs = [], [], []
    warnings = []
    old_output_variables = list(model.output_variables)
    try:
        # CatMAP computes DRC for the product via finite energy perturbations.
        # rate and coverage are emitted for the SAME exact catalyst coordinates.
        model.output_variables = ['coverage', 'production_rate', 'rate', 'rate_control']
        for surface_i, surface in enumerate(surfaces):
            coords = []
            for descriptor in descriptors:
                values = model.species_definitions[descriptor].get('formation_energy')
                if not isinstance(values, (list, tuple)) or len(values) != len(surfaces):
                    raise ValueError(f'{descriptor}: invalid per-surface formation energy list')
                coords.append(float(_number(values[surface_i], f'{surface} {descriptor}')))
            try:
                model.mapper.get_point_output(coords)
            except Exception as exc:
                raise RuntimeError(f'CatMAP single-point DRC/TOF FAILED at {surface}: {coords}') from exc

            cov_labels = list(model.output_labels.get('coverage', []))
            prod_labels = list(model.output_labels.get('production_rate', []))
            rc_labels = model.output_labels.get('rate_control', [])
            if len(rc_labels) != 2:
                raise ValueError(f'{surface}: malformed DRC labels: {rc_labels!r}')
            rc_gases, rc_species = map(list, rc_labels)
            idx_oo = _labels_index(cov_labels, OOH_NAME, 'coverage')
            idx_c6 = _labels_index(cov_labels, ALKENE_NAME, 'coverage')
            idx_to = _labels_index(prod_labels, product, 'production_rate')
            idx_rc_gas = _labels_index(rc_gases, product, 'rate_control gas axis')
            idx_rc_ooh = _labels_index(rc_species, OOH_NAME, 'rate_control energy axis')
            idx_rc_ts1 = _labels_index(rc_species, TS1_NAME, 'rate_control energy axis')
            coverages = list(model._coverage)
            rates = list(model._rate)
            drc = model._rate_control
            if len(coverages) != len(cov_labels) or len(rates) != len(rxns):
                raise ValueError(f'{surface}: unexpected coverage or step-rate dimensions')
            if len(drc) != len(rc_gases) or len(drc[idx_rc_gas]) != len(rc_species):
                raise ValueError(f'{surface}: DRC array / label dimensions disagree')

            cvg_values = [_number(x, f'{surface} coverage {name}')
                          for name, x in zip(cov_labels, coverages)]
            theta = cvg_values[idx_oo]
            theta_c6 = cvg_values[idx_c6]
            sum_theta = sum(cvg_values)
            empty = site_total - sum_theta  # inferred, NOT independently measured
            r2, r4, r5 = [_number(rates[i], f'{surface} reaction {i+1}') for i in (1, 3, 4)]
            tof = _number(model._production_rate[idx_to], f'{surface} product TOF')
            x_ooh = _number(drc[idx_rc_gas][idx_rc_ooh], f'{surface} X_OOH')
            x_ts1 = _number(drc[idx_rc_gas][idx_rc_ts1], f'{surface} X_TS1')
            all_drc = [_number(q, f'{surface} DRC member') for q in drc[idx_rc_gas]]
            drc_ok = any(q != 0 for q in all_drc) and all(abs(q) <= 10 for q in all_drc)

            site_ok = (min(cvg_values) >= -mp.mpf('1e-8') and
                       max(cvg_values) <= 1+mp.mpf('1e-8') and
                       empty >= -mp.mpf('1e-8') and empty <= 1+mp.mpf('1e-8'))
            # Three independent rate comparisons, NOT the tautological site sum.
            gap24 = _rel_mismatch(r2, r4)
            gap45 = _rel_mismatch(r4, r5)
            gap5tof = _rel_mismatch(r5, tof)
            max_gap = max(gap24, gap45, gap5tof)
            flux_ok = max_gap <= mp.mpf('1e-3') and r2 >= 0 and r4 >= 0 and r5 >= 0
            logtof = _text(mp.log10(tof)) if tof > 0 else 'UNDEFINED_TOF_ZERO'
            label = _classify(theta, x_ooh, x_ts1, r4, tof, site_ok, flux_ok, drc_ok)
            if not site_ok or not flux_ok or not drc_ok:
                warnings.append(f'{surface}: {label}')
            if tof > 0:
                pairs.append((_text(theta), _text(tof), surface))
            legacy_tof_rows.append([surface, *map(_text, coords), _text(tof)])
            legacy_pairs.append([surface, *map(_text, coords), _text(theta), _text(tof), logtof])
            all_rows.append([
                surface, *map(_text, coords), _text(theta), _text(theta_c6),
                _text(sum_theta), _text(empty), _text(tof), logtof,
                _text(x_ooh), _text(x_ts1), _text(r2), _text(r4), _text(r5),
                _text(r4-tof), _text(gap24), _text(gap45), _text(gap5tof),
                _text(max_gap), 'PASS' if site_ok else 'FAIL',
                'PASS' if flux_ok else 'FAIL',
                'CHECK' if not drc_ok else 'FINITE_NONZERO', _interpret_tags(label), label,
            ])
    finally:
        model.output_variables = old_output_variables

    # Preserve the previously published filenames and columns.
    _save_csv(output_dir / f'TOF_{product}_grid.csv',
              [*descriptors, f'{product}_production_rate_per_site_s-1'], grid_rows)
    _save_csv(output_dir / f'TOF_{product}_catalysts.csv',
              ['surface_name', *descriptors, f'{product}_production_rate_per_site_s-1'], legacy_tof_rows)
    _save_csv(output_dir / f'OOH_coverage_log10TOF_{product}_catalysts.csv',
              ['surface_name', *descriptors, 'theta_OOH_s',
               f'{product}_production_rate_per_site_s-1', 'log10_TOF_per_site_s-1'], legacy_pairs)
    _save_csv(output_dir / f'OOH_kinetic_evidence_{product}_14sites.csv', [
        'surface_name', *descriptors, 'theta_OOH_s', 'theta_C6H12_s',
        'sum_adsorbate_coverages', 'theta_star_inferred',
        'C6H12O_product_TOF_per_site_s-1', 'log10_TOF',
        'X_OOH_s_DRC', 'X_TS1_OOH-C6H12_s_DRC',
        'r2_OOH_formation_net_s-1', 'r4_oxygen_transfer_net_s-1',
        'r5_epoxide_desorption_net_s-1', 'r4_minus_product_TOF_s-1',
        'rel_gap_r2_r4', 'rel_gap_r4_r5', 'rel_gap_r5_TOF', 'max_rel_flux_gap',
        'site_occupancy_check', 'flux_closure_check', 'DRC_numeric_screen',
        'mechanistic_interpretation_cn', 'mechanistic_evidence_tags'], all_rows)

    fit = _fit_xy(pairs)
    fit_path = output_dir / f'OOH_coverage_log10TOF_{product}_regression.csv'
    if fit:
        keys = ['n', 'slope', 'intercept', 'R2', 'RMSE_log10', 'LOOCV_Q2']
        _save_csv(fit_path, keys, [[fit[k] for k in keys]])
        try:
            _fit_plot(pairs, fit, output_dir / f'OOH_coverage_log10TOF_{product}.pdf')
        except Exception as exc:
            warnings.append(f'Coverage-logTOF PDF: {exc}')
    else:
        _save_csv(fit_path, ['status'], [['FIT_UNDEFINED']])

    count_accum = sum('OOH_ACCUMULATION' in row[-1] for row in all_rows)
    count_sens = sum('OOH_ENERGY_SENSITIVE_CANDIDATE' in row[-1] for row in all_rows)
    count_ts1 = sum('TS1_RATE_CONTROL_CANDIDATE' in row[-1] for row in all_rows)
    count_ok = sum(row[-5] == 'PASS' and row[-4] == 'PASS' for row in all_rows)
    report = [
        'CatMAP OOH kinetic evidence — preliminary, NOT a proof of physical accessibility',
        f'Catalysts: {len(all_rows)}',
        f'OOH* coverage >= 0.9: {count_accum}',
        f'|X_OOH| >= 0.1 with nonzero, finite, bounded DRC row: {count_sens}',
        f'X_TS1 >= 0.5 with nonzero, finite, bounded DRC row: {count_ts1}',
        f'Passed simple coverage plausibility + relative reaction flux checks: {count_ok}',
        f"OOH coverage-logTOF all-site R2: {fit['R2'] if fit else 'UNDEFINED'}",
        f"OOH coverage-logTOF LOOCV Q2: {fit['LOOCV_Q2'] if fit else 'UNDEFINED'}",
        '',
        'INTERPRETATION:',
        '* OOH_ACCUMULATION denotes abundance, NOT rate control.',
        '* OXYGEN_TRANSFER_PRECURSOR_PATH denotes reaction-network participation and flux match, NOT exclusive control.',
        '* OOH_ENERGY_SENSITIVE_CANDIDATE flags numerical CatMAP DRC sensitivity, NOT robust independent proof.',
        '* TS1_RATE_CONTROL_CANDIDATE means TS1 rate sensitivity; it is not X_OOH.',
        '* Multiple labels are allowed; the three roles are not mutually exclusive.',
        '* theta_star_inferred = 1 - sum_adsorbate_coverages; this is a definition, not an independent site balance verification.',
        '* Site screen checks physical ranges. Flux comparisons r2/r4, r4/r5, r5/TOF give additional checks.',
        '* Steady-state reachability and sensitivity to initial inventory/perturbations remain UNVERIFIED.',
        '* DRC numerical screen only checks gross anomalies (all zeros or abs>10), not perturbation convergence.',
        '* TOF here is CatMAP product production_rate, not a measured mass-normalized rate.',
        '', 'WARNINGS:', *warnings]
    (output_dir / f'OOH_kinetic_evidence_{product}_summary.txt').write_text(
        '\n'.join(report)+'\n', encoding='utf-8')
    print('\n'.join(report[:9]))
    print('Combined table:', output_dir / f'OOH_kinetic_evidence_{product}_14sites.csv')
    return (output_dir / f'TOF_{product}_grid.csv',
            output_dir / f'TOF_{product}_catalysts.csv')
