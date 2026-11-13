#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Analyze CatMAP 0.3.1 transient coverages without changing model inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import traceback
from datetime import datetime
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
EXPECTED_REACTIONS = [
    [["H2O2_g", "s", "s"], ["Ha_s", "OOH_s"]],
    [["C6H12_g", "s"], ["C6H12_s"]],
    [["OOH_s", "C6H12_s", "s"], ["OOH-C6H12_s", "s", "s"],
     ["C6H12O_s", "O_s", "Hb_s"]],
    [["C6H12O_s"], ["C6H12O_g", "s"]],
    [["Hb_s", "O_s"], ["Hb-O_s", "s"], ["OH_s", "s"]],
    [["Ha_s", "OH_s"], ["s", "H2O_s"]],
    [["H2O_s"], ["H2O_g", "s"]],
]


def _add_local_catmap_source(project_dir: Path) -> None:
    """Use the project's checked-in CatMAP source when available."""
    candidates = [
        project_dir.parent.parent / "catmap_official_source",
        SCRIPT_DIR.parent.parent / "catmap_official_source",
    ]
    env_source = os.environ.get("CATMAP_SOURCE")
    if env_source:
        candidates.insert(0, Path(env_source))
    for candidate in candidates:
        if (candidate / "catmap" / "__init__.py").is_file():
            sys.path.insert(0, str(candidate))
            return


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return math.nan


def _unique_path(parent: Path, prefix: str) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = parent / f"{prefix}_{stamp}"
    suffix = 1
    while candidate.exists():
        candidate = parent / f"{prefix}_{stamp}_{suffix:02d}"
        suffix += 1
    candidate.mkdir(parents=True)
    return candidate


def _load_maps(model, pkl_path: Path):
    import pickle

    with pkl_path.open("rb") as handle:
        data = pickle.load(handle)
    needed = ("coverage_map", "turnover_frequency_map", "rate_map")
    missing = [key for key in needed if key not in data]
    if missing:
        raise ValueError(f"CatMAP pkl 缺少完整网格映射: {missing}")
    if not all(len(data[key]) == 2500 for key in needed):
        raise ValueError("预期 50×50 原始网格，但 pkl 映射数量不是 2500")
    if len(model.gas_names) != len(data["turnover_frequency_map"][0][1]):
        raise ValueError("turnover_frequency_map 通道数与 CatMAP gas_names 不一致")
    coverage_width = len(data["coverage_map"][0][1])
    n_surface_sites = len([site for site in model.site_names if site != "g"])
    if coverage_width not in (len(model.adsorbate_names),
                              len(model.adsorbate_names) + n_surface_sites):
        raise ValueError("coverage_map 通道数无法由 CatMAP 吸附物种和位点标签解释")
    if coverage_width > len(model.adsorbate_names):
        # CatMAP's numbers solver appends site-population variables after the
        # adsorbate vector. For this one-site LH model, that final value is the
        # vacancy coverage (theta_star); validate it against site balance.
        for point, values in data["coverage_map"]:
            ads = sum(_float(v) for v in values[:len(model.adsorbate_names)])
            sites = sum(_float(v) for v in values[len(model.adsorbate_names):])
            if not math.isclose(ads + sites, 1.0, rel_tol=1e-6, abs_tol=1e-6):
                raise ValueError(f"CatMAP pkl 附加位点覆盖度未通过位点守恒核验: {point}")
    return data


def _prepare_official_solver(model):
    """Apply CatMAP 0.3.1 ReactionModel.run numeric setup without running a map."""
    import mpmath as mp

    if model.numerical_representation != "mpmath":
        raise ValueError("910 的模型未使用 mpmath；不静默切换数值表示")
    mp.mp.dps = model.decimal_precision + 25
    model._math = mp
    model._mpfloat = mp.mpf
    model._matrix = mp.matrix
    model._kB = mp.mpf("8.617332478e-5")
    model._h = mp.mpf("4.135667516e-15")
    model._math.infnorm = lambda values: mp.norm(values, "inf")
    model._math.lsqnorm = lambda values: mp.norm(values, 2) ** 2
    # CatMAP compiles its own kinetic expressions on first mapper use. Compile
    # them directly to avoid calling model.run(), which would remap/write files.
    model.solver.compile()


def _build_grid(model, data, product_name):
    import numpy as np

    if product_name not in model.gas_names:
        raise ValueError(f"{product_name} 不在 CatMAP gas_names={model.gas_names}")
    channel = model.gas_names.index(product_name)
    entries = data["turnover_frequency_map"]
    xs = np.array(sorted({float(point[0]) for point, _ in entries}), dtype=float)
    ys = np.array(sorted({float(point[1]) for point, _ in entries}), dtype=float)
    if (len(xs), len(ys)) != (50, 50):
        raise ValueError(f"原始描述符网格尺寸为 {len(xs)}×{len(ys)}，预期 50×50")
    x_index = {value: index for index, value in enumerate(xs)}
    y_index = {value: index for index, value in enumerate(ys)}
    tof = np.full((len(ys), len(xs)), np.nan)
    for point, values in entries:
        tof[y_index[float(point[1])], x_index[float(point[0])]] = _float(values[channel])
    if not np.isfinite(tof).all():
        raise ValueError("原始 TOF 网格存在缺失或非有限数据")
    return xs, ys, tof, x_index, y_index, channel


def _line_cells(a, b):
    import numpy as np

    # All candidate points are native grid nodes; this samples only existing
    # cells along their straight grid segment and creates no interpolated data.
    count = max(abs(a[0] - b[0]), abs(a[1] - b[1]), 1) * 4 + 1
    cells = set()
    for t in np.linspace(0.0, 1.0, count):
        row = int(round(a[0] + t * (b[0] - a[0])))
        col = int(round(a[1] + t * (b[1] - a[1])))
        cells.add((row, col))
    return sorted(cells)


def _auto_points(xs, ys, tof):
    import numpy as np

    positive = np.isfinite(tof) & (tof > 0)
    if positive.sum() == 0:
        raise ValueError("C6H12O_g 原始 TOF 网格没有正值，无法识别双峰")
    log_tof = np.full(tof.shape, -np.inf)
    log_tof[positive] = np.log10(tof[positive])
    peaks = []
    rows, cols = tof.shape
    for row in range(rows):
        for col in range(cols):
            if not positive[row, col]:
                continue
            local = log_tof[max(0, row - 1):min(rows, row + 2),
                            max(0, col - 1):min(cols, col + 2)]
            if log_tof[row, col] >= np.max(local):
                peaks.append((float(log_tof[row, col]), row, col))
    peaks.sort(reverse=True)

    pair = None
    best_key = None
    for i, first in enumerate(peaks):
        for second in peaks[i + 1:]:
            _, r1, c1 = first
            _, r2, c2 = second
            if abs(c1 - c2) < 2 or (r1 == r2 and c1 == c2):
                continue
            cells = _line_cells((r1, c1), (r2, c2))
            valley_cell = min(cells, key=lambda rc: log_tof[rc])
            valley_log = float(log_tof[valley_cell])
            if valley_log >= min(first[0], second[0]):
                continue
            # Prefer two distinct, high native maxima with a real sampled dip.
            prominence = min(first[0], second[0]) - valley_log
            key = (min(first[0], second[0]), first[0] + second[0], prominence)
            if best_key is None or key > best_key:
                pair, best_key = (first, second, valley_cell), key
    if pair is None:
        raise ValueError("原始 50×50 TOF 网格中不能确认两个独立局部峰及其间谷底")

    first, second, valley_cell = pair
    left, right = sorted((first, second), key=lambda item: xs[item[2]])
    points = []
    for label, peak in (("peak_left", left), ("peak_right", right)):
        _, row, col = peak
        points.append({"label": label, "row": row, "col": col,
                       "ha": float(xs[col]), "o": float(ys[row]),
                       "tof": float(tof[row, col])})
    row, col = valley_cell
    points.append({"label": "valley", "row": row, "col": col,
                   "ha": float(xs[col]), "o": float(ys[row]),
                   "tof": float(tof[row, col])})
    return points


def _manual_points(spec, xs, ys, tof):
    points = []
    values = [part.strip() for part in spec.split(",")]
    if len(values) != 6:
        raise ValueError("--manual-points 需要六个数值: Ha1,O1,Ha2,O2,Ha3,O3")
    pairs = [(float(values[i]), float(values[i + 1])) for i in range(0, 6, 2)]
    labels = ("peak_left", "peak_right", "valley")
    for label, (ha, oxygen) in zip(labels, pairs):
        col = int(abs(xs - ha).argmin())
        row = int(abs(ys - oxygen).argmin())
        points.append({"label": label, "row": row, "col": col,
                       "ha": ha, "o": oxygen,
                       "tof": float(tof[row, col]),
                       "nearest_ha": float(xs[col]), "nearest_o": float(ys[row]),
                       "manual_coordinate": True})
    return points


def _index_lookup(mapping):
    return {tuple(float(v) for v in point): values for point, values in mapping}


def _grid_reference(data, point, xs, ys):
    coordinate = (float(xs[point["col"]]), float(ys[point["row"]]))
    return coordinate


def _integrate_point(model, descriptors, args):
    import numpy as np
    from scipy.integrate import solve_ivp

    solver = model.solver
    rxn_parameters = model.scaler.get_rxn_parameters(list(descriptors))
    # CatMAP numbers-solver state = adsorbates followed by non-gas site
    # populations. This model has one site variable, initially theta_star=1.
    n_site_variables = len([site for site in model.site_names if site != "g"])
    if n_site_variables != 1:
        raise ValueError(f"洁净初值实现仅按当前单表面位点模型核验，解析到位点变量数={n_site_variables}")
    y0 = np.zeros(len(model.adsorbate_names) + n_site_variables, dtype=float)
    y0[len(model.adsorbate_names)] = 1.0
    solver._rxn_parameters = rxn_parameters
    solver._descriptors = list(descriptors)

    def coverage_and_derivative(y):
        coverage_mp = [solver._mpfloat(float(value)) for value in y]
        rate_constants = solver.get_rate_constants(rxn_parameters, coverage_mp)
        # get_rate_constants updates CatMAP's _kf/_kr; the CatMAP-generated
        # steady-state function is exactly the d(theta)/dt vector.
        derivative = solver.ideal_steady_state_function(coverage_mp)
        return coverage_mp, rate_constants, np.asarray([float(v) for v in derivative], dtype=float)

    def sampled_rates(coverage_mp, rate_constants):
        rates = solver.get_rxn_rates(coverage_mp, rate_constants)
        tof = solver.get_turnover_frequency(
            rxn_parameters, rates=rates, verify_coverages=False)
        return (np.asarray([float(v) for v in rates], dtype=float),
                np.asarray([float(v) for v in tof], dtype=float))

    def rhs(_time, y):
        _, _, derivative = coverage_and_derivative(y)
        return derivative

    def jacobian(_time, y):
        coverage_mp, _, _ = coverage_and_derivative(y)
        jac = solver.ideal_steady_state_jacobian(coverage_mp)
        return np.asarray([[float(jac[i, j]) for j in range(jac.cols)]
                           for i in range(jac.rows)], dtype=float)

    positive_times = np.geomspace(args.first_time, args.max_time,
                                  max(2, args.samples - 1))
    t_eval = np.concatenate(([0.0], positive_times))
    result = solve_ivp(rhs, (0.0, args.max_time), y0, method="Radau", jac=jacobian,
                       t_eval=t_eval, rtol=args.rtol, atol=args.atol,
                       max_step=args.max_step if args.max_step else np.inf)
    if result.t.size == 0:
        raise RuntimeError("Radau 没有返回任何时间点")
    if not result.success:
        raise RuntimeError(f"Radau 积分失败: {result.message}")

    derivs, rates_rows, tof_rows = [], [], []
    for coverage in result.y.T:
        coverage_mp, rate_constants, derivative = coverage_and_derivative(coverage)
        rates, tof = sampled_rates(coverage_mp, rate_constants)
        derivs.append(derivative)
        rates_rows.append(rates)
        tof_rows.append(tof)
    return result, np.asarray(derivs), np.asarray(rates_rows), np.asarray(tof_rows)


def _diagnostics(model, result, derivatives, rates, tof_rows, point, raw_coverage,
                 raw_tof, args):
    import numpy as np

    y = result.y.T
    rate_scale = np.max(np.abs(rates), axis=1)
    rate_scale = np.maximum(rate_scale, np.finfo(float).tiny)
    mean_rate = np.mean(rates, axis=1)
    flux_error = np.max(np.abs(rates - mean_rate[:, None]), axis=1) / rate_scale
    max_derivative = np.max(np.abs(derivatives), axis=1)
    normalized_derivative = max_derivative / rate_scale
    site_sum = np.sum(y, axis=1)
    hb_idx = model.adsorbate_names.index("Hb_s")
    o_idx = model.adsorbate_names.index("O_s")
    hb_o_diff = y[:, hb_idx] - y[:, o_idx]
    site_error = np.abs(site_sum - np.minimum(site_sum, 1.0))
    stable = "UNVERIFIED"
    stability_eigenvalues = []
    if result.success:
        try:
            solver = model.solver
            final = y[-1]
            cvg_mp = [solver._mpfloat(float(v)) for v in final]
            constants = solver.get_rate_constants(solver._rxn_parameters, cvg_mp)
            # Ensure _kf/_kr correspond to the final coverage before Jacobian.
            jac = solver.ideal_steady_state_jacobian(cvg_mp)
            matrix = np.asarray([[float(jac[i, j]) for j in range(jac.cols)]
                                 for i in range(jac.rows)], dtype=float)
            eig = np.linalg.eigvals(matrix)
            stability_eigenvalues = [complex(v) for v in eig]
            scale = max(1.0, float(np.max(np.abs(eig))))
            if np.max(np.real(eig)) < -1e-9 * scale:
                stable = "STABLE"
            elif np.max(np.real(eig)) > 1e-9 * scale:
                stable = "UNSTABLE"
            else:
                stable = "NEUTRAL_OR_MARGINAL"
        except Exception:
            stable = "UNVERIFIED: Jacobian calculation failed"

    endpoint = y[-1]
    endpoint_tof = tof_rows[-1]
    # CatMAP output maps store gas TOFs; this array was channel-aligned earlier.
    steady = bool(normalized_derivative[-1] <= args.steady_rate_tol and
                  flux_error[-1] <= args.flux_tol and
                  np.min(endpoint) >= -args.bound_tol and
                  np.max(site_sum) <= 1.0 + args.bound_tol and
                  np.min(site_sum) >= -args.bound_tol and
                  np.min(y) >= -args.bound_tol and
                  abs(hb_o_diff[-1]) <= args.bound_tol)
    return {
        "time": result.t, "coverage": y, "rates": rates, "tofs": tof_rows,
        "derivatives": derivatives, "max_abs_dtheta_dt": max_derivative,
        "normalized_dtheta_dt": normalized_derivative,
        "flux_closure_error": flux_error, "hb_o_difference": hb_o_diff,
        "site_sum": site_sum, "site_error": site_error,
        "steady": steady, "stability": stable,
        "stability_eigenvalues": stability_eigenvalues,
        "endpoint_tof": endpoint_tof,
    }


def _write_point_outputs(output_dir, point, model, diagnostics, raw_coverage,
                         raw_tof, gas_names, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    output_dir.mkdir(parents=True, exist_ok=False)
    time = diagnostics["time"]
    coverage = diagnostics["coverage"]
    rates = diagnostics["rates"]
    csv_path = output_dir / "transient_data.csv"
    header = ["time", "theta_star"]
    header += ["theta_" + name.removesuffix("_s") for name in model.adsorbate_names]
    header += [f"R{i}" for i in range(1, len(model.elementary_rxns) + 1)]
    header += ["max_abs_dtheta_dt", "normalized_dtheta_dt", "flux_closure_error",
               "hb_o_coverage_difference", "site_coverage_sum"]
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for i, t in enumerate(time):
            writer.writerow([f"{t:.17g}", f"{coverage[i, len(model.adsorbate_names)]:.17g}"] +
                            [f"{v:.17g}" for v in coverage[i, :len(model.adsorbate_names)]] +
                            [f"{v:.17g}" for v in rates[i]] +
                            [f"{diagnostics['max_abs_dtheta_dt'][i]:.17g}",
                             f"{diagnostics['normalized_dtheta_dt'][i]:.17g}",
                             f"{diagnostics['flux_closure_error'][i]:.17g}",
                             f"{diagnostics['hb_o_difference'][i]:.17g}",
                             f"{diagnostics['site_sum'][i]:.17g}"])

    labels = [name.replace("_s", "*") for name in model.adsorbate_names]
    fig, ax = plt.subplots(figsize=(10, 6))
    for i, label in enumerate(labels):
        ax.plot(time, coverage[:, i], label=label)
    ax.set(xlabel="Time (s)", ylabel="Coverage", ylim=(0, 1),
           title=f"{point['label']}: coverage vs time")
    ax.grid(True, alpha=.25); ax.legend(ncol=2, fontsize=8)
    fig.tight_layout(); fig.savefig(output_dir / "coverage_vs_time.png", dpi=180); plt.close(fig)
    if time[-1] / max(time[1], np.finfo(float).tiny) > 1e3:
        fig, ax = plt.subplots(figsize=(10, 6))
        for i, label in enumerate(labels):
            ax.plot(time[1:], coverage[1:, i], label=label)
        ax.set_xscale("log"); ax.set(xlabel="Time (s, log scale)", ylabel="Coverage",
                                      ylim=(0, 1), title=f"{point['label']}: log-time coverage")
        ax.grid(True, which="both", alpha=.25); ax.legend(ncol=2, fontsize=8)
        fig.tight_layout(); fig.savefig(output_dir / "coverage_vs_time_log.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(rates.shape[1]):
        ax.plot(time, rates[:, i], label=f"R{i + 1}")
    ax.set(xlabel="Time (s)", ylabel="Net reaction rate", title=f"{point['label']}: net rates")
    ax.set_yscale("symlog", linthresh=1e-8)
    ax.grid(True, alpha=.25); ax.legend(ncol=2)
    fig.tight_layout(); fig.savefig(output_dir / "reaction_rates_vs_time.png", dpi=180); plt.close(fig)
    if time[-1] / max(time[1], np.finfo(float).tiny) > 1e3:
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(rates.shape[1]):
            ax.plot(time[1:], rates[1:, i], label=f"R{i + 1}")
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1e-8)
        ax.set(xlabel="Time (s, log scale)", ylabel="Net reaction rate",
               title=f"{point['label']}: net rates on log time")
        ax.grid(True, which="both", alpha=.25); ax.legend(ncol=2)
        fig.tight_layout(); fig.savefig(output_dir / "reaction_rates_vs_time_log.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time, diagnostics["max_abs_dtheta_dt"], label="max(abs(dtheta/dt))")
    ax.plot(time, diagnostics["normalized_dtheta_dt"], label="normalized max(abs(dtheta/dt))")
    ax.plot(time, diagnostics["flux_closure_error"], label="normalized flux closure")
    ax.set(xlabel="Time (s)", ylabel="Residual / normalized error",
           title=f"{point['label']}: steady-state diagnostics")
    ax.set_yscale("symlog", linthresh=1e-30)
    ax.grid(True, which="both", alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(output_dir / "steady_state_residual_vs_time.png", dpi=180); plt.close(fig)

    max_cover_diff = float(np.max(np.abs(diagnostics["coverage"][-1] - raw_coverage)))
    min_cov = float(np.min(coverage)); max_site = float(np.max(diagnostics["site_sum"]))
    max_hbo = float(np.max(np.abs(diagnostics["hb_o_difference"])))
    product_index = gas_names.index("C6H12O_g")
    raw_product_tof = float(raw_tof[product_index])
    return {
        "point": point, "raw_product_tof": raw_product_tof,
        "transient_product_tof": float(diagnostics["endpoint_tof"][product_index]),
        "max_abs_dtheta_dt": float(diagnostics["max_abs_dtheta_dt"][-1]),
        "normalized_dtheta_dt": float(diagnostics["normalized_dtheta_dt"][-1]),
        "flux_closure_error": float(diagnostics["flux_closure_error"][-1]),
        "site_conservation_error": float(max(0.0, max_site - 1.0)),
        "hb_o_conservation_error": max_hbo,
        "max_coverage_difference_vs_raw_catmap": max_cover_diff,
        "min_coverage": min_cov, "max_site_coverage_sum": max_site,
        "steady": diagnostics["steady"], "stability": diagnostics["stability"],
        "final_time": float(diagnostics["time"][-1]),
        "coverage_names": list(model.adsorbate_names) + ["theta_star"],
        "transient_coverage": diagnostics["coverage"][-1].tolist(),
        "raw_coverage": raw_coverage.tolist(),
        "raw_tof_all_gases": raw_tof.tolist(),
        "transient_tof_all_gases": diagnostics["endpoint_tof"].tolist(),
        "stability_eigenvalues": [str(v) for v in diagnostics["stability_eigenvalues"]],
    }


def _write_summary(output_dir, summaries, model, points, args, catmap_version,
                   pkl_path, mkm_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    csv_path = output_dir / "transient_summary.csv"
    columns = ["label", "descriptor_Ha", "descriptor_O", "grid_row", "grid_col",
               "raw_TOF_C6H12O_g", "transient_final_TOF_C6H12O_g",
               "max_abs_dtheta_dt", "flux_closure_error", "site_conservation_error",
               "hb_o_conservation_error", "max_coverage_difference_vs_raw_catmap",
               "final_time_s", "steady_state", "local_stability"]
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for item in summaries:
            p = item["point"]
            writer.writerow({"label": p["label"], "descriptor_Ha": p["ha"],
                             "descriptor_O": p["o"], "grid_row": p["row"],
                             "grid_col": p["col"], "raw_TOF_C6H12O_g": item["raw_product_tof"],
                             "transient_final_TOF_C6H12O_g": item["transient_product_tof"],
                             "max_abs_dtheta_dt": item["max_abs_dtheta_dt"],
                             "flux_closure_error": item["flux_closure_error"],
                             "site_conservation_error": item["site_conservation_error"],
                             "hb_o_conservation_error": item["hb_o_conservation_error"],
                             "max_coverage_difference_vs_raw_catmap": item["max_coverage_difference_vs_raw_catmap"],
                             "final_time_s": item["final_time"], "steady_state": item["steady"],
                             "local_stability": item["stability"]})

    labels = [item["point"]["label"] for item in summaries]
    fig, (ax, table_ax) = plt.subplots(1, 2, figsize=(12, 5),
                                       gridspec_kw={"width_ratios": [1, 1.25]})
    for item in summaries:
        p = item["point"]
        ax.scatter(p["ha"], p["o"], s=60, label=p["label"])
    ax.set(xlabel="Ha descriptor (eV)", ylabel="O descriptor (eV)", title="Selected raw-grid points")
    ax.grid(True, alpha=.25); ax.legend()
    table_ax.axis("off")
    table_rows = []
    for item in summaries:
        cov = dict(zip(item["coverage_names"], item["transient_coverage"]))
        table_rows.append([item["point"]["label"], f"{cov['OOH_s']:.3e}", f"{cov['O_s']:.3e}"])
    table_ax.table(cellText=table_rows, colLabels=["Point", "OOH* endpoint", "O* endpoint"],
                   loc="center", cellLoc="center")
    table_ax.set_title("Transient endpoint coverages")
    fig.tight_layout(); fig.savefig(output_dir / "comparison_coverage.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 6))
    raw = [item["raw_product_tof"] for item in summaries]
    transient = [item["transient_product_tof"] for item in summaries]
    index = np.arange(len(labels)); width = .38
    ax.scatter(index - width / 2, raw, s=70, marker="o", label="Raw CatMAP grid TOF")
    ax.scatter(index + width / 2, transient, s=70, marker="s", label="Transient endpoint TOF")
    ax.set_xticks(index, labels)
    all_values = np.asarray(raw + transient, dtype=float)
    nonzero_abs = np.abs(all_values[np.nonzero(all_values)])
    if np.all(all_values > 0):
        ax.set_yscale("log")
        ax.set_ylim(max(np.min(all_values) / 3.0, np.finfo(float).tiny), np.max(all_values) * 3.0)
    else:
        linear_threshold = max(float(np.min(nonzero_abs)) if nonzero_abs.size else 1e-30,
                               np.finfo(float).tiny)
        ax.set_yscale("symlog", linthresh=linear_threshold)
    ax.set(ylabel="C6H12O_g TOF", title="Raw CatMAP vs transient endpoint TOF")
    ax.grid(True, axis="y", alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(output_dir / "comparison_tof.png", dpi=180); plt.close(fig)

    def sha256(path):
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    energy_path = mkm_path.parent / "energies_dft_neb.txt"
    report = [
        "CatMAP transient analysis report",
        f"Created: {datetime.now().isoformat(timespec='seconds')}",
        f"CatMAP version: {catmap_version}",
        f"Source mkm: {mkm_path}", f"Source pkl: {pkl_path}",
        f"Source energy table: {energy_path}",
        f"SHA256 mkm: {sha256(mkm_path)}", f"SHA256 pkl: {sha256(pkl_path)}",
        f"SHA256 energies_dft_neb.txt: {sha256(energy_path)}",
        f"Temperature (K): {model.temperature}",
        f"Gas pressures/activities (CatMAP order): {list(model.gas_pressures)}",
        f"Adsorbates (CatMAP order): {list(model.adsorbate_names)}",
        f"Gas species (CatMAP order): {list(model.gas_names)}",
        f"Elementary reactions: {len(model.elementary_rxns)}",
        "Integrator: scipy.integrate.solve_ivp(method='Radau')",
        "RHS: CatMAP SteadyStateSolver.get_rate_constants + ideal_steady_state_function; rates: get_rxn_rates",
        "State arithmetic: SciPy float64; each state is passed to CatMAP as mpf(float64). No clipping or renormalization is applied.",
        f"rtol={args.rtol}; atol={args.atol}; first_time={args.first_time}; max_time={args.max_time}; samples={args.samples}",
        f"Steady criteria: max|dtheta/dt|/max(|R1..R7|) <= {args.steady_rate_tol}; normalized flux closure <= {args.flux_tol}; physical-bound tolerance={args.bound_tol}",
        "",
    ]
    for item in summaries:
        p = item["point"]
        raw_ref = item.get("raw_reference_coordinate", [p["ha"], p["o"]])
        report += [f"[{p['label']}] Ha={p['ha']:.12g}, O={p['o']:.12g}, grid=({p['row']},{p['col']}), raw reference=({raw_ref[0]:.12g},{raw_ref[1]:.12g}), raw TOF={item['raw_product_tof']:.8e}",
                   f"  final time={item['final_time']:.8e} s; status={'STEADY' if item['steady'] else 'NOT_CONVERGED'}; local stability={item['stability']}",
                   f"  transient endpoint TOF={item['transient_product_tof']:.8e}; max|dtheta/dt|={item['max_abs_dtheta_dt']:.8e}; normalized derivative={item['normalized_dtheta_dt']:.8e}; flux error={item['flux_closure_error']:.8e}",
                   f"  max site sum={item['max_site_coverage_sum']:.8e}; min coverage={item['min_coverage']:.8e}; max |Hb-O|={item['hb_o_conservation_error']:.8e}",
                   f"  max coverage difference vs raw CatMAP={item['max_coverage_difference_vs_raw_catmap']:.8e}",
                   "  coverages by CatMAP label:"]
        report += [f"    {name}: transient={tr:.8e}, raw={raw:.8e}"
                   for name, tr, raw in zip(item["coverage_names"], item["transient_coverage"], item["raw_coverage"])]
        report.append("")
    report += ["Physical interpretation: this report describes numerical integration of the unchanged CatMAP model.",
               "A returned endpoint is not called steady unless both derivative and flux criteria pass.",
               "Local Jacobian eigenvalues are a linear check; neutral modes are reported as marginal."]
    (output_dir / "analysis_report.txt").write_text("\n".join(report), encoding="utf-8")
    (output_dir / "analysis_summary.json").write_text(json.dumps(summaries, ensure_ascii=False, indent=2), encoding="utf-8")


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", help="CatMAP run directory; relative paths resolve from this script's folder")
    parser.add_argument("--manual-points", help="six comma-separated Ha,O coordinates replacing automatic peak detection")
    parser.add_argument("--product", default="C6H12O_g")
    parser.add_argument("--max-time", type=float, default=1e8)
    parser.add_argument("--first-time", type=float, default=1e-15)
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-14)
    parser.add_argument("--max-step", type=float, default=0.0)
    parser.add_argument("--steady-rate-tol", type=float, default=1e-4)
    parser.add_argument("--flux-tol", type=float, default=1e-3)
    parser.add_argument("--bound-tol", type=float, default=1e-8)
    return parser.parse_args()


def main():
    import numpy as np
    args = _parse_args()
    project_dir = SCRIPT_DIR
    runs_dir = project_dir / "runs"
    if args.run_dir:
        run_dir = Path(args.run_dir)
        if not run_dir.is_absolute():
            run_dir = project_dir / run_dir
    else:
        candidates = sorted((p for p in runs_dir.iterdir() if p.is_dir()),
                            key=lambda p: p.name, reverse=True)
        if not candidates:
            raise FileNotFoundError(f"runs 目录下没有运行记录: {runs_dir}")
        print("可用 CatMAP 运行记录:")
        for index, candidate in enumerate(candidates, 1):
            print(f"{index}. {candidate.name}")
        choice = int(input("请选择运行编号: "))
        if choice < 1 or choice > len(candidates):
            raise ValueError("运行编号超出范围")
        run_dir = candidates[choice - 1]
    run_dir = run_dir.resolve()
    mkm_path = run_dir / "alkene_epoxidation.mkm"
    pkl_path = run_dir / "alkene_epoxidation_dft.pkl"
    if not mkm_path.is_file() or not pkl_path.is_file():
        raise FileNotFoundError(f"运行目录缺少 alkene_epoxidation.mkm 或 pkl: {run_dir}")
    _add_local_catmap_source(project_dir)
    from catmap import ReactionModel, __version__
    if __version__ != "0.3.1":
        raise RuntimeError(f"本脚本按 CatMAP 0.3.1 源码核验，当前导入版本为 {__version__}")

    # CatMAP reads data_file relative to cwd. Change only the process cwd; no
    # model.run() is called and the selected source run remains read-only.
    os.chdir(run_dir)
    model = ReactionModel(setup_file=str(mkm_path))
    _prepare_official_solver(model)
    if len(model.elementary_rxns) != 7:
        raise ValueError(f"模型解析出 {len(model.elementary_rxns)} 步反应，预期 7；停止分析")
    if [list(map(list, reaction)) for reaction in model.elementary_rxns] != EXPECTED_REACTIONS:
        raise ValueError("当前 .mkm 解析出的 LH 反应网络与 910 约定不一致；停止分析")
    if tuple(model.adsorbate_names) != ("C6H12O_s", "C6H12_s", "H2O_s", "Ha_s", "Hb_s", "OH_s", "OOH_s", "O_s"):
        raise ValueError(f"解析到的吸附物种顺序/集合不符合本任务预期: {model.adsorbate_names}")
    data = _load_maps(model, pkl_path)
    xs, ys, tof_grid, x_index, y_index, product_channel = _build_grid(model, data, args.product)
    points = (_manual_points(args.manual_points, xs, ys, tof_grid) if args.manual_points
              else _auto_points(xs, ys, tof_grid))
    raw_coverage = _index_lookup(data["coverage_map"])
    raw_tof_map = _index_lookup(data["turnover_frequency_map"])
    output_parent = run_dir / "descriptor_Ha_O" / "transient_analysis"
    output_parent.mkdir(parents=True, exist_ok=True)
    output_dir = _unique_path(output_parent, "analysis")

    summaries = []
    try:
        for point in points:
            print(f"RUNNING {point['label']} at Ha={point['ha']:.6g}, O={point['o']:.6g}", flush=True)
            native = _grid_reference(data, point, xs, ys)
            if native not in raw_coverage or native not in raw_tof_map:
                raise ValueError(f"描述符节点没有原始覆盖度/TOF映射: {native}")
            point_dir = output_dir / point["label"]
            descriptors = [point["ha"], point["o"]]
            # For manual off-grid coordinates, the transient uses the exact
            # requested coordinates; the original CatMAP comparison is nearest node.
            if point.get("manual_coordinate"):
                native_descriptors = [native[0], native[1]]
            else:
                native_descriptors = descriptors
            raw_cvg = np.asarray([_float(v) for v in raw_coverage[native]], dtype=float)
            raw_tof = np.asarray([_float(v) for v in raw_tof_map[native]], dtype=float)
            result, derivs, rates, tof_rows = _integrate_point(model, descriptors, args)
            diag = _diagnostics(model, result, derivs, rates, tof_rows, point,
                                raw_cvg, raw_tof, args)
            summary = _write_point_outputs(point_dir, point, model, diag,
                                           raw_cvg, raw_tof, model.gas_names, args)
            if native_descriptors != descriptors:
                summary["raw_reference_coordinate"] = native_descriptors
            summaries.append(summary)
            print(f"FINISHED {point['label']}: {'STEADY' if summary['steady'] else 'NOT_CONVERGED'} "
                  f"at t={summary['final_time']:.4g} s", flush=True)
        _write_summary(output_dir, summaries, model, points, args, __version__, pkl_path, mkm_path)
    except Exception:
        (output_dir / "analysis_error_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise
    print(f"分析结果目录: {output_dir}")
    for item in summaries:
        print(f"{item['point']['label']}: {'STEADY' if item['steady'] else 'NOT_CONVERGED'}; "
              f"Ha={item['point']['ha']:.6g}, O={item['point']['o']:.6g}, "
              f"final time={item['final_time']:.4g} s, max|dtheta/dt|={item['max_abs_dtheta_dt']:.3e}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        failure = traceback.format_exc()
        try:
            args = _parse_args()
            selected = Path(args.run_dir) if args.run_dir else SCRIPT_DIR
            if not selected.is_absolute():
                selected = SCRIPT_DIR / selected
            error_dir = selected / "descriptor_Ha_O" / "transient_analysis"
            error_dir.mkdir(parents=True, exist_ok=True)
            name = "analysis_error_" + datetime.now().strftime("%Y%m%d_%H%M%S") + ".txt"
            (error_dir / name).write_text(failure, encoding="utf-8")
        except Exception:
            pass
        print(failure, file=sys.stderr)
        raise
