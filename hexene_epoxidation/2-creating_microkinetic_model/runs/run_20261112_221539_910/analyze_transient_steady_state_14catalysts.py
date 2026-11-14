#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""对 910 运行中的 14 个真实催化剂进行 CatMAP 瞬态覆盖度分析。"""

import argparse
import csv
import os
import pickle
import traceback
from datetime import datetime
from pathlib import Path

import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent
SPECIES = [
    "OOH_s", "O_s", "Hb_s", "Ha_s",
    "C6H12_s", "C6H12O_s", "OH_s", "H2O_s",
]


def load_maps(pkl_file):
    if not pkl_file.exists():
        return [], []
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)
    return data.get("coverage_map", []), data.get("turnover_frequency_map", [])


def nearest_map(mapping, descriptors):
    if not mapping:
        return None, None
    point, values = min(
        mapping,
        key=lambda item: np.linalg.norm(
            np.asarray(item[0], dtype=float) - np.asarray(descriptors, dtype=float)
        ),
    )
    distance = float(np.linalg.norm(
        np.asarray(point, dtype=float) - np.asarray(descriptors, dtype=float)
    ))
    return np.asarray(values, dtype=float), distance


def load_model(mkm_file, output_dir):
    from catmap import ReactionModel

    model = ReactionModel(setup_file=str(mkm_file))

    # 仅初始化 CatMAP 原生 Parser / Scaler / Solver，不执行描述符网格扫描。
    model.descriptor_ranges = None
    model.descriptors = None
    model.descriptor_values = None
    model.output_variables = []
    model.data_file = str(output_dir / "transient_init.pkl")
    model.log_file = str(output_dir / "transient_init.log")
    model.run()

    # CatMAP 只在 mapper（get_point_output / get_output_map）内部调用 solver.compile()。
    # 这里跳过网格扫描，必须显式编译动态函数；且必须放在 model.run() 之后，
    # 因为模板需要 run() 里设置的 _kB/_h/_mpfloat/_matrix。
    # 未编译时 ReactionModelWrapper.__getattr__ 会把这些名字返回 None，
    # 调用即报 TypeError: 'NoneType' object is not callable。
    model.solver.compile()

    # hasattr 对这些名字恒为真（包装器返回 None），必须用 callable 校验。
    required = ("rate_constants", "ideal_mean_field_steady_state",
                "elementary_rates", "constrain_coverage_function")
    missing = [name for name in required
               if not callable(getattr(model.solver, name, None))]
    if missing:
        raise RuntimeError(f"CatMAP Solver 动态函数编译失败：{missing}")
    return model


def transient_one(model, descriptors, args):
    solver = model.solver
    with mp.workdps(args.transient_dps):
        return _transient_one_impl(model, descriptors, args)


def _transient_one_impl(model, descriptors, args):
    solver = model.solver
    ads = list(model.adsorbate_names)
    missing = [s for s in SPECIES if s not in ads]
    if missing:
        raise ValueError(f"模型缺少吸附物种：{missing}")

    # use_numbers_solver=True（CatMAP 0.3.1 默认）时空位是自由变量：
    # 状态向量 = 吸附物种 + 非气相位点，编译后的 theta_length = n_ads + n_sites - 1，
    # 速率表达式中的空位写作 theta[n_ads]。洁净表面初值只能是 [0]*n_ads + [1]。
    n_ads = len(ads)
    if getattr(model, "use_numbers_solver", False):
        site_vars = [s for s in model.site_names if s != "g"]
        if len(site_vars) != 1:
            raise ValueError(
                f"洁净初值实现仅适配单一位点模型，位点变量数={len(site_vars)}")
        theta0 = [mp.mpf("0") for _ in ads] + [mp.mpf("1") for _ in site_vars]
    else:
        theta0 = [mp.mpf("0") for _ in ads]

    rxn_parameters = model.scaler.get_rxn_parameters(list(descriptors))
    solver._rxn_parameters = rxn_parameters
    solver._descriptors = list(descriptors)

    # 理想平均场、无吸附质相互作用：速率常数与覆盖度无关，只算一次并常驻 solver._kf/_kr。
    rate_constants = solver.get_rate_constants(rxn_parameters, theta0)

    def rhs_mp(theta):
        # 返回 dtheta/dt（含空位分量），保留 mpmath 精度用于稳态判定。
        return solver.ideal_steady_state_function(theta)

    # ---- 用刚性求解器 scipy Radau 替代 mpmath odefun ----
    # mpmath.odefun 是非刚性算法，面对微动力学刚性系统会退化成近乎无限的微小步长而“卡死”；
    # CatMAP 已编译解析雅可比 ideal_steady_state_jacobian，交给 scipy Radau 可用隐式牛顿迭代
    # + 雅可比，在少量步内跨过刚性瞬态。
    from scipy.integrate import solve_ivp

    def rhs_np(t, y):
        theta = [mp.mpf(v) for v in y]
        return np.array([float(v) for v in rhs_mp(theta)])

    def jac_np(t, y):
        theta = [mp.mpf(v) for v in y]
        J = solver.ideal_steady_state_jacobian(theta)
        return np.asarray(J.tolist(), dtype=float)

    y0 = np.array([float(v) for v in theta0])
    sample_times = np.geomspace(args.first_time, args.max_time, args.samples - 1)
    sol = solve_ivp(
        rhs_np, (0.0, float(args.max_time)), y0,
        method=args.ivp_method, jac=jac_np,
        rtol=args.ivp_rtol, atol=args.ivp_atol, t_eval=sample_times,
    )
    if not sol.success:
        raise RuntimeError(f"solve_ivp({args.ivp_method}) 未收敛：{sol.message}")

    times = [mp.mpf("0")] + [mp.mpf(str(t)) for t in sample_times]
    theta_series = [list(theta0)] + [
        [mp.mpf(str(float(v))) for v in sol.y[:, k]] for k in range(sol.y.shape[1])
    ]

    rows = []
    residuals_mp, flux_errors_mp = [], []
    last_theta = list(theta0)
    for t, theta in zip(times, theta_series):
        theta = list(theta)
        last_theta = list(theta)
        rates = solver.get_rxn_rates(theta, rate_constants)
        derivative = rhs_mp(theta)

        # 保留 mpmath 精度完成稳态判定，最后才转换为 CSV 数值。
        residual_mp = max(abs(v) for v in derivative)
        rates_mp = list(rates)
        rate_scale_mp = max([abs(v) for v in rates_mp] + [mp.mpf("1e-300")])
        flux_error_mp = max(
            abs(v - sum(rates_mp) / len(rates_mp)) for v in rates_mp
        ) / rate_scale_mp
        normalized_residual_mp = residual_mp / max(rate_scale_mp, mp.mpf("1"))

        ads_sum = float(sum(theta[:n_ads]))
        if len(theta) > n_ads:
            vacancy = float(theta[n_ads])
        else:
            vacancy = 1.0 - ads_sum
        residuals_mp.append(residual_mp)
        flux_errors_mp.append(flux_error_mp)
        row = {
            "time_s": float(t),
            "theta_vacant": vacancy,
            "max_abs_dtheta_dt": float(residual_mp),
            "normalized_dtheta_dt": float(normalized_residual_mp),
            "flux_closure_relative": float(flux_error_mp),
            "Hb_minus_O": float(
                theta[ads.index("Hb_s")] - theta[ads.index("O_s")]
            ),
            "site_occupancy_sum": ads_sum,
            "site_balance_error": abs(1.0 - (ads_sum + vacancy)),
        }
        for species in SPECIES:
            row["theta_" + species.replace("_s", "")] = float(theta[ads.index(species)])
        for i, rate in enumerate(rates_mp, 1):
            row[f"R{i}"] = float(rate)
        rows.append(row)

    # 稳态判据：相对残差（|dθ/dt| / max|r|）达标即可，绝对残差 <= tolerance 作为兜底。
    # 直接用 .mkm 的 tolerance（1e-140）作绝对判据在 mpmath 下几乎不可达。
    steady_tol_mp = mp.mpf(str(args.steady_tol))
    tail_n = min(3, len(rows))
    steady = all(
        (rows[-i]["normalized_dtheta_dt"] <= args.steady_tol
         or residuals_mp[-i] <= model.tolerance)
        and flux_errors_mp[-i] <= mp.mpf(str(args.flux_tol))
        and rows[-i]["theta_vacant"] >= -args.site_tol
        and min(rows[-i]["theta_" + s.replace("_s", "")] for s in SPECIES)
            >= -args.site_tol
        and rows[-i]["site_balance_error"] <= args.site_tol
        for i in range(1, tail_n + 1)
    )

    final_theta = [mp.mpf(v) for v in last_theta]
    final_constants = solver.get_rate_constants(rxn_parameters, final_theta)
    final_rates = solver.get_rxn_rates(final_theta, final_constants)
    final_tofs = solver.get_turnover_frequency(
        rxn_parameters, rates=final_rates, verify_coverages=False
    )
    product_index = list(model.gas_names).index("C6H12O_g")

    return rows, {
        "steady_status": "STEADY" if steady else "NOT_CONVERGED",
        "final_TOF_C6H12O_g": float(final_tofs[product_index]),
        "final_OOH_coverage": rows[-1]["theta_OOH"],
        "final_O_coverage": rows[-1]["theta_O"],
        "final_Hb_coverage": rows[-1]["theta_Hb"],
        "final_vacancy_coverage": rows[-1]["theta_vacant"],
        "final_residual": rows[-1]["max_abs_dtheta_dt"],
        "final_normalized_residual": rows[-1]["normalized_dtheta_dt"],
        "final_flux_closure_relative": rows[-1]["flux_closure_relative"],
        "minimum_vacancy_coverage": min(r["theta_vacant"] for r in rows),
        "minimum_adsorbate_coverage": min(
            r["theta_" + s.replace("_s", "")] for r in rows for s in SPECIES
        ),
        "maximum_site_occupancy": max(r["site_occupancy_sum"] for r in rows),
        "maximum_site_balance_error": max(r["site_balance_error"] for r in rows),
        "maximum_abs_Hb_minus_O": max(abs(r["Hb_minus_O"]) for r in rows),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # 归档后脚本已位于 run_XXX 目录内部，默认直接使用脚本所在目录；
    # 若需在外部指定其他 run，可显式传入绝对路径或相对于脚本目录的路径。
    parser.add_argument("--run-dir", default=".")
    parser.add_argument("--max-time", type=float, default=1e3)
    parser.add_argument("--first-time", type=float, default=1e-15)
    parser.add_argument("--samples", type=int, default=120)
    parser.add_argument("--ode-tol", default="1e-20")
    parser.add_argument("--transient-dps", type=int, default=40)
    parser.add_argument("--ivp-method", default="Radau", choices=["Radau", "BDF", "LSODA"])
    parser.add_argument("--ivp-rtol", type=float, default=1e-8)
    parser.add_argument("--ivp-atol", type=float, default=1e-12)
    parser.add_argument("--flux-tol", type=float, default=1e-6)
    parser.add_argument("--site-tol", type=float, default=1e-10)
    parser.add_argument("--steady-tol", type=float, default=1e-8)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_absolute():
        run_dir = ROOT / run_dir
    run_dir = run_dir.resolve()
    mkm_file = run_dir / "alkene_epoxidation.mkm"
    pkl_file = run_dir / "alkene_epoxidation_dft.pkl"
    energy_file = run_dir / "energies_dft_neb.txt"
    for path in (mkm_file, energy_file):
        if not path.is_file():
            raise FileNotFoundError(path)

    os.chdir(run_dir)
    output_dir = run_dir / "descriptor_Ha_O" / "transient_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir /= datetime.now().strftime("analysis_%Y%m%d_%H%M%S")
    output_dir.mkdir()

    coverage_map, tof_map = load_maps(pkl_file)
    model = load_model(mkm_file, output_dir)

    if len(model.elementary_rxns) != 7:
        raise ValueError(f"模型解析出 {len(model.elementary_rxns)} 步反应，预期为 LH 七步。")
    if len(model.surface_names) != 14:
        raise ValueError(f"模型读取到 {len(model.surface_names)} 个催化剂，预期 14 个。")
    if list(model.descriptor_names) != ["Ha_s", "O_s"]:
        raise ValueError(f"描述符与预期不符：{model.descriptor_names}")

    # 描述符直接来自当前 .mkm / energies_dft_neb.txt 解析出的真实 formation_energy。
    descriptors_by_surface = {}
    for i, surface in enumerate(model.surface_names):
        pair = []
        for name in model.descriptor_names:
            values = model.species_definitions[name].get("formation_energy")
            if values is None or len(values) != len(model.surface_names):
                raise ValueError(f"{name} 的 formation_energy 与催化剂数量不匹配。")
            pair.append(float(values[i]))
        descriptors_by_surface[surface] = pair

    pdf_file = output_dir / "transient_coverage_14catalysts.pdf"
    all_rows, summaries = [], []

    with PdfPages(pdf_file) as pdf:
        for number, (surface, descriptors) in enumerate(descriptors_by_surface.items(), 1):
            print(f"[{number}/14] {surface}: Ha={descriptors[0]:.6f}, O={descriptors[1]:.6f}")
            try:
                rows, result = transient_one(model, descriptors, args)
                for row in rows:
                    row["surface"] = surface
                all_rows.extend(rows)

                raw_cov, cov_distance = nearest_map(coverage_map, descriptors)
                raw_tof, tof_distance = nearest_map(tof_map, descriptors)
                if raw_cov is not None and len(raw_cov) == len(model.adsorbate_names):
                    final_cov = np.asarray(
                        [rows[-1]["theta_" + s.replace("_s", "")] for s in model.adsorbate_names]
                    )
                    result["nearest_pkl_coverage_max_abs_diff"] = float(
                        np.max(np.abs(final_cov - raw_cov))
                    )
                    result["nearest_pkl_descriptor_distance"] = cov_distance
                if raw_tof is not None:
                    result["nearest_pkl_TOF_C6H12O_g"] = float(
                        raw_tof[list(model.gas_names).index("C6H12O_g")]
                    )
                    result["nearest_pkl_TOF_descriptor_distance"] = tof_distance

                summary = {
                    "surface": surface, "Ha_eV": descriptors[0], "O_eV": descriptors[1],
                    **result,
                }
                summaries.append(summary)

                fig, ax = plt.subplots(figsize=(9, 5.5))
                t = np.asarray([r["time_s"] for r in rows])
                for species in SPECIES:
                    y = np.asarray([
                        r["theta_" + species.replace("_s", "")] for r in rows
                    ])
                    ax.plot(t[1:], y[1:], label=species.replace("_s", "*"))
                ax.set_xscale("log")
                ax.set_ylim(-0.02, 1.02)
                ax.set_xlabel("Time (s)")
                ax.set_ylabel("Surface coverage")
                ax.set_title(
                    f"{surface} | Ha={descriptors[0]:.3f} eV, O={descriptors[1]:.3f} eV | "
                    f"{result['steady_status']}"
                )
                ax.grid(alpha=0.25)
                ax.legend(ncol=4, fontsize=8)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

                print("  ", result["steady_status"],
                      "OOH =", result["final_OOH_coverage"],
                      "O =", result["final_O_coverage"])
            except Exception as exc:
                traceback.print_exc()
                summaries.append({
                    "surface": surface, "Ha_eV": descriptors[0], "O_eV": descriptors[1],
                    "steady_status": "FAILED", "error": f"{type(exc).__name__}: {exc}",
                })
                print(f"  FAILED: {type(exc).__name__}: {exc}")

    if all_rows:
        with open(output_dir / "transient_data.csv", "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_rows[0]))
            writer.writeheader()
            writer.writerows(all_rows)

    fields = sorted({key for row in summaries for key in row})
    with open(output_dir / "transient_summary.csv", "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summaries)

    print("\n输出目录：", output_dir)
    print("PDF：", pdf_file)
    print("CSV：transient_data.csv、transient_summary.csv")
    for row in summaries:
        print(row["surface"], row["steady_status"])


if __name__ == "__main__":
    main()
