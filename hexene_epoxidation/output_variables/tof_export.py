"""CatMAP product-rate export with paired OOH* coverage and honest linear diagnostics.

Drop-in replacement for hexene_epoxidation/output_variables/tof_export.py.
The already-connected mkm_job.py calls export_tof_tables(model, output_dir, product='C6H12O_g').
No modification to DFT data, reaction network, scaler, solver or model parameters.
"""
from pathlib import Path
import csv
from math import isfinite

from mpmath import mp


def _write_csv(path, header, rows):
    with Path(path).open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def _linear_diagnostics(records):
    """Fit log10(rate / s^-1) = slope * theta_OOH + intercept, all valid surfaces."""
    if len(records) < 3:
        return None, "Insufficient positive-TOF points (< 3)."
    with mp.workdps(80):
        xy = [(mp.mpf(str(cvg)), mp.log10(mp.mpf(str(rate)))) for cvg, rate, _ in records]
        n = len(xy)
        xm = sum(x for x, _ in xy) / n
        ym = sum(y for _, y in xy) / n
        sxx = sum((x - xm) ** 2 for x, _ in xy)
        syy = sum((y - ym) ** 2 for _, y in xy)
        if sxx == 0 or syy == 0:
            return None, "Zero coverage or log-TOF variance; linear fit undefined."
        slope = sum((x - xm) * (y - ym) for x, y in xy) / sxx
        intercept = ym - slope * xm
        err = sum((y - slope * x - intercept) ** 2 for x, y in xy)
        r2 = 1 - err / syy
        rmse = mp.sqrt(err / n)
        # Leave-one-out predictions: all surfaces count as the assessment cohort.
        loo_sse = mp.mpf("0")
        for i in range(n):
            train = xy[:i] + xy[i + 1:]
            tx = sum(p[0] for p in train) / (n - 1)
            ty = sum(p[1] for p in train) / (n - 1)
            txx = sum((x - tx) ** 2 for x, _ in train)
            if txx == 0:
                return None, "LOOCV cannot be evaluated: identical training coverages."
            ta = sum((x - tx) * (y - ty) for x, y in train) / txx
            tb = ty - ta * tx
            loo_sse += (xy[i][1] - ta * xy[i][0] - tb) ** 2
        q2 = 1 - loo_sse / syy
        return {
            "n": n,
            "slope": mp.nstr(slope, 22),
            "intercept": mp.nstr(intercept, 22),
            "R2": mp.nstr(r2, 14),
            "RMSE_log10": mp.nstr(rmse, 14),
            "LOOCV_Q2": mp.nstr(q2, 14),
        }, None


def _plot_relation(pairs, fit, file_path):
    """Plot only observed points and their fit; display x-axis scaling is not a model change."""
    import matplotlib.pyplot as plt
    with mp.workdps(80):
        x = [mp.mpf(str(item[0])) for item in pairs]
        y = [float(mp.log10(mp.mpf(str(item[1])))) for item in pairs]
        xmax = max(abs(value) for value in x)
        if xmax == 0:
            return
        exponent = int(mp.floor(mp.log10(xmax))) if xmax < mp.mpf("0.000001") else 0
        factor = mp.power(10, exponent)
        xp = [float(t / factor) for t in x]
        a, b = mp.mpf(fit["slope"]), mp.mpf(fit["intercept"])
        order = sorted(range(len(xp)), key=xp.__getitem__)
        fit_x = [xp[i] for i in order]
        fit_y = [float(a * x[i] + b) for i in order]
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    try:
        ax.scatter(xp, y, s=36, zorder=3)
        ax.plot(fit_x, fit_y, linewidth=1.5, label="OLS (all valid catalysts)")
        for xx, yy, item in zip(xp, y, pairs):
            ax.annotate(item[2], (xx, yy), textcoords="offset points", xytext=(4, 3), fontsize=8)
        ax.set_xlabel(r"OOH* coverage ($\theta_{\mathrm{OOH}^*}$)" +
                      (f"   [display units: x $10^{{{exponent}}}$]" if exponent else ""))
        ax.set_ylabel(r"$\log_{10}[\mathrm{TOF}_{\mathrm{C_6H_{12}O}}/(1\,\mathrm{s}^{-1})]$")
        ax.text(.02, .98, f"n={fit['n']}  R²={float(fit['R2']):.4f}  LOOCV Q²={float(fit['LOOCV_Q2']):.4f}",
                ha="left", va="top", transform=ax.transAxes, fontsize=9)
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        fig.savefig(file_path)
    finally:
        plt.close(fig)


def export_tof_tables(model, output_dir, product="C6H12O_g", coverage_species="OOH_s"):
    """Export CatMAP grid, exact catalyst TOF and exact matching OOH* coverage."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    names = list(model.descriptor_names)
    if len(names) != 2:
        raise ValueError("Expected two descriptor_names")
    gas_labels = list(model.output_labels.get("production_rate", []))
    if product not in gas_labels:
        raise ValueError(f"Product {product} not found in production_rate labels {gas_labels}")
    gas_index = gas_labels.index(product)
    rate_map = getattr(model, "production_rate_map", None)
    if not rate_map:
        raise ValueError("CatMAP production_rate_map is missing or empty")

    grid_rows = []
    for coords, rates in rate_map:
        if len(coords) != 2 or gas_index >= len(rates):
            raise ValueError("Malformed production_rate_map")
        grid_rows.append([str(coords[0]), str(coords[1]), str(rates[gas_index])])

    selected_surfaces = list(model.surface_names)
    catalyst_rows, paired_rows, valid_pairs = [], [], []
    old_outputs = list(model.output_variables)
    try:
        model.output_variables = ["coverage", "production_rate"]
        for surf_i, surf in enumerate(selected_surfaces):
            coords = []
            for name in names:
                values = model.species_definitions[name].get("formation_energy")
                if not isinstance(values, (tuple, list)) or len(values) != len(selected_surfaces):
                    raise ValueError(f"Invalid formation energies for {name}")
                coordinate = float(values[surf_i])
                if not isfinite(coordinate):
                    raise ValueError(f"Nonfinite descriptor {name} for {surf}")
                coords.append(coordinate)
            try:
                model.mapper.get_point_output(coords)  # same physical state for TOF and OOH coverage
            except Exception as exc:
                raise RuntimeError(f"CatMAP single-point solve FAILED: {surf}, {coords}") from exc
            current_gas_labels = list(model.output_labels["production_rate"])
            current_coverage_labels = list(model.output_labels["coverage"])
            if product not in current_gas_labels or coverage_species not in current_coverage_labels:
                raise ValueError("Requested product or coverage species missing from CatMAP labels")
            rate = model._production_rate[current_gas_labels.index(product)]
            cvg = model._coverage[current_coverage_labels.index(coverage_species)]
            rate_mp, cvg_mp = mp.mpf(str(rate)), mp.mpf(str(cvg))
            if not mp.isfinite(rate_mp) or not mp.isfinite(cvg_mp):
                raise ValueError(f"Nonfinite TOF/coverage for {surf}")
            if cvg_mp < -mp.mpf("1e-15") or cvg_mp > 1 + mp.mpf("1e-15"):
                raise ValueError(f"Unphysical OOH coverage for {surf}: {cvg}")
            rate_s, cvg_s = str(rate), str(cvg)
            catalyst_rows.append([surf, str(coords[0]), str(coords[1]), rate_s])
            if rate_mp > 0:
                logtof = mp.nstr(mp.log10(rate_mp), 25)
                valid_pairs.append((cvg_s, rate_s, surf))
            else:
                logtof = "UNDEFINED_TOF_ZERO"
            paired_rows.append([surf, str(coords[0]), str(coords[1]), cvg_s, rate_s, logtof])
    finally:
        model.output_variables = old_outputs

    # Preserve original CSV names, contents and ordering.
    grid_file = output_dir / f"TOF_{product}_grid.csv"
    catalyst_file = output_dir / f"TOF_{product}_catalysts.csv"
    paired_file = output_dir / f"OOH_coverage_log10TOF_{product}_catalysts.csv"
    fit_file = output_dir / f"OOH_coverage_log10TOF_{product}_regression.csv"
    figure_file = output_dir / f"OOH_coverage_log10TOF_{product}.pdf"
    _write_csv(grid_file, [*names, f"{product}_production_rate_per_site_s-1"], grid_rows)
    _write_csv(catalyst_file, ["surface_name", *names,
                               f"{product}_production_rate_per_site_s-1"], catalyst_rows)
    _write_csv(paired_file, ["surface_name", *names, f"theta_{coverage_species}",
                            f"{product}_production_rate_per_site_s-1", "log10_TOF_per_site_s-1"], paired_rows)
    fit, reason = _linear_diagnostics(valid_pairs)
    if fit is None:
        _write_csv(fit_file, ["status", "reason", "number_of_positive_tof_points"],
                   [["UNAVAILABLE", reason, len(valid_pairs)]])
        print(f"OOH-coverage versus logTOF fit unavailable: {reason}")
    else:
        _write_csv(fit_file, ["n", "slope", "intercept", "R2", "RMSE_log10", "LOOCV_Q2"],
                   [[fit[k] for k in ("n", "slope", "intercept", "R2", "RMSE_log10", "LOOCV_Q2")]])
        print(f"OOH coverage/log10 TOF fit (all valid surfaces): R2={fit['R2']}, LOOCV_Q2={fit['LOOCV_Q2']}")
        try:
            _plot_relation(valid_pairs, fit, figure_file)
        except Exception as exc:
            print(f"WARNING: paired CSV and fit CSV saved, but scatter PDF failed: {exc}")
    print(f"CatMAP TOF exports: {grid_file}, {catalyst_file}, {paired_file}, {fit_file}")
    return grid_file, catalyst_file
