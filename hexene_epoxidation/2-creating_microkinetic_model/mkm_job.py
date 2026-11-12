from catmap import ReactionModel, analyze
from pathlib import Path
from math import isfinite
import shutil
import random
import csv
import time
from datetime import datetime


# ─────────────────────────────────────────────
#  工具函数
# ─────────────────────────────────────────────
def auto_descriptor_ranges(model, margin=0.05):
    ranges = []
    for name in model.descriptor_names:
        data = model.species_definitions[name].get("formation_energy")
        if not isinstance(data, (list, tuple)):
            raise ValueError(f"{name}: 缺少逐催化剂能量数据")
        if len(data) != len(model.surface_names):
            raise ValueError(f"{name}: 能量数量与催化剂数量不一致")

        values = []
        for surface, energy in zip(model.surface_names, data):
            try:
                value = float(energy)
            except (ValueError, TypeError):
                raise ValueError(f"{surface}: {name} 缺少有效能量")
            if not isfinite(value):
                raise ValueError(f"{surface}: {name} 能量不是有限数值")
            values.append(value)

        low, high = min(values), max(values)
        if high <= low:
            raise ValueError(f"{name}: 数据没有有效跨度，无法自动确定范围")

        padding = (high - low) * margin
        ranges.append([low - padding, high + padding])
    return ranges


def save_run_to_archive(output_dir):
    runs_dir = Path(__file__).resolve().parent / 'runs'
    runs_dir.mkdir(parents=True, exist_ok=True)

    now = datetime.now()
    rand_suffix = random.randint(0, 999)
    run_name = f"run_{now.strftime('%Y%m%d_%H%M%S')}_{rand_suffix:03d}"
    run_dir = runs_dir / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    base = Path(__file__).resolve().parent
    script_file = Path(__file__).resolve()
    for f in base.iterdir():
        if f.is_file():
            shutil.copy2(f, run_dir / f.name)
    shutil.copytree(output_dir, run_dir / output_dir.name, dirs_exist_ok=True)
    print(f"运行结果已保存到: {run_dir}")
    print(f"脚本快照已保存: {run_dir / script_file.name}")


def to_float(v):
    try:
        return float(v)
    except (OverflowError, ValueError):
        return 0.0


def export_tof_and_coverage(model, output_dir):
    tof_cvg_file = output_dir / 'TOF_and_coverage.csv'

    surface_descriptors = {}
    for i, surface in enumerate(model.surface_names):
        desc_vals = []
        for desc_name in model.descriptor_names:
            energies = model.species_definitions[desc_name].get("formation_energy")
            if energies and i < len(energies):
                desc_vals.append(float(energies[i]))
            else:
                desc_vals.append(None)
        if all(v is not None for v in desc_vals):
            surface_descriptors[surface] = tuple(desc_vals)

    tof_lookup = {}
    for pt, tof_vals in model.turnover_frequency_map:
        key = tuple(to_float(x) for x in pt)
        tof_lookup[key] = [to_float(v) for v in tof_vals]
    cvg_lookup = {}
    for pt, cvg_vals in model.coverage_map:
        key = tuple(to_float(x) for x in pt)
        cvg_lookup[key] = [to_float(v) for v in cvg_vals]

    def find_nearest(target, lookup):
        if not lookup:
            return None
        best_key = None
        best_dist = float('inf')
        for key in lookup:
            dist = sum((a - b) ** 2 for a, b in zip(target, key))
            if dist < best_dist:
                best_dist = dist
                best_key = key
        return lookup[best_key]

    with open(tof_cvg_file, 'w', newline='') as f:
        writer = csv.writer(f)

        header = ['surface']
        header += ['descriptor_' + n.removesuffix('_s') for n in model.descriptor_names]
        header += ['TOF_' + g for g in model.gas_names]
        header += ['coverage_' + a for a in model.adsorbate_names]
        writer.writerow(header)

        count = 0
        for surface, desc_vals in surface_descriptors.items():
            row = [surface]
            row += [f"{v:.6f}" for v in desc_vals]

            tof_data = find_nearest(desc_vals, tof_lookup)
            cvg_data = find_nearest(desc_vals, cvg_lookup)

            if tof_data:
                row += [f"{v:.6e}" for v in tof_data]
            else:
                row += ['N/A'] * len(model.gas_names)

            if cvg_data:
                row += [f"{v:.6e}" for v in cvg_data]
            else:
                row += ['N/A'] * len(model.adsorbate_names)

            writer.writerow(row)
            count += 1

    print(f"\nTOF 和覆盖率数据已保存到: {tof_cvg_file}")
    print(f"共导出 {count} 个催化剂表面的数据")

    print("\n── 各催化剂 TOF 摘要 ──")
    for surface, desc_vals in surface_descriptors.items():
        tof_data = find_nearest(desc_vals, tof_lookup)
        if tof_data:
            parts = [f"{g}={v:.3e}" for g, v in zip(model.gas_names, tof_data) if v > 0]
            tof_str = ", ".join(parts) if parts else "全部为0"
            print(f"  {surface}: {tof_str}")

    print("\n── 各催化剂覆盖率摘要 ──")
    for surface, desc_vals in surface_descriptors.items():
        cvg_data = find_nearest(desc_vals, cvg_lookup)
        if cvg_data:
            parts = [f"{a}={v:.3e}" for a, v in zip(model.adsorbate_names, cvg_data) if v > 1e-30]
            cvg_str = ", ".join(parts) if parts else "可忽略"
            print(f"  {surface}: {cvg_str}")


# ─────────────────────────────────────────────
#  模型加载与配置
# ─────────────────────────────────────────────
mkm_file = 'alkene_epoxidation.mkm'
model = ReactionModel(setup_file=mkm_file)

descriptor_pair = "_".join(name.removesuffix("_s") for name in model.descriptor_names)
output_dir = Path(__file__).resolve().parent / f"descriptor_{descriptor_pair}"
output_dir.mkdir(parents=True, exist_ok=True)

print(f"当前描述符: {model.descriptor_names}")
print(f"图表输出目录: {output_dir}")

_t_total_start = time.perf_counter()
_t_checkpoints = [('模型加载', _t_total_start)]

def _t_mark(label):
    _t_checkpoints.append((label, time.perf_counter()))

model.descriptor_ranges = auto_descriptor_ranges(model, margin=0.05)
_t_mark('描述符范围')
print("自动描述符范围:")
for name, bounds in zip(model.descriptor_names, model.descriptor_ranges):
    print(f"  {name}: {bounds}")


# ─────────────────────────────────────────────
#  第一次运行：基础速率图
# ─────────────────────────────────────────────
model.run()
_t_mark('第1次运行(run)')

vm = analyze.VectorMap(model)

vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-140
vm.max = 1e-1
vm.plot(save=str(output_dir / 'rate.pdf'))

vm.unique_only = False
vm.plot(save=str(output_dir / 'all_rates.pdf'))
vm.unique_only = True
_t_mark('基础速率图')


# ─────────────────────────────────────────────
#  第二次运行：production_rate + turnover_frequency + rate_control
# ─────────────────────────────────────────────
model.output_variables += ['production_rate', 'turnover_frequency', 'rate_control']
model.run()
_t_mark('第2次运行(run)')

# ── 手动构建 rate_control_map（mapper 缺少 get_rate_control_map 方法）──
import numpy as _np
print("正在计算 rate_control_map（DRC）...")
_resolution = model.resolution
if hasattr(_resolution, '__len__'):
    _resx, _resy = int(_resolution[0]), int(_resolution[1])
else:
    _resx = _resy = int(_resolution)
_d1_vals = _np.linspace(model.descriptor_ranges[0][0], model.descriptor_ranges[0][1], _resx)
_d2_vals = _np.linspace(model.descriptor_ranges[1][0], model.descriptor_ranges[1][1], _resy)
_rate_control_map = []
_total_pts = _resx * _resy
for _i, _d1 in enumerate(_d1_vals):
    for _j, _d2 in enumerate(_d2_vals):
        _desc = [_d1, _d2]
        _params = model.scaler.get_rxn_parameters(_desc)
        _rc = model.solver.get_rate_control(_params)
        _rate_control_map.append([_desc, _rc])
        if (_i * _resy + _j + 1) % 100 == 0:
            print(f"  DRC 进度: {_i * _resy + _j + 1}/{_total_pts}")
model.rate_control_map = _rate_control_map
model.output_labels['rate_control'] = [model.adsorbate_names, model.elementary_rxns]
print(f"rate_control_map 构建完成: {len(_rate_control_map)} 个点")
_t_mark('DRC map构建')

vm.production_rate_map = model.production_rate_map
vm.threshold = 1e-30
vm.include_empty = True
vm.plot_variable = 'production_rate'
vm.min = 1e-80
vm.max = 1e-1
vm.log_scale = True
vm.include_descriptors = False
vm.plot(save=str(output_dir / 'production_rate.pdf'))

vm.descriptor_labels = [
    f"{name.removesuffix('_s')} reactivity [eV]"
    for name in model.descriptor_names
]
vm.subplots_adjust_kwargs = {'left': 0.2, 'right': 0.8, 'bottom': 0.15}
vm.plot(save=str(output_dir / 'pretty_production_rate.pdf'))

vm.plot_variable = 'coverage'
vm.log_scale = False
vm.min = 0
vm.max = 1
vm.plot(save=str(output_dir / 'coverage.pdf'))

coverage_descriptor = model.descriptor_names[0]
vm.include_labels = [coverage_descriptor]
vm.plot(save=str(output_dir / f"{coverage_descriptor.removesuffix('_s')}_coverage.pdf"))

sa = analyze.ScalingAnalysis(model)
sa.plot(save=str(output_dir / 'scaling.pdf'))
_t_mark('生产速率与覆盖率图')


# ─────────────────────────────────────────────
#  DRC 分析（Degree of Rate Control）
# ─────────────────────────────────────────────
mm = analyze.MatrixMap(model)
mm.plot_variable = 'rate_control'
mm.log_scale = False
mm.min = -2
mm.max = 2
mm.unique_only = False
mm.descriptor_labels = [
    f"{name.removesuffix('_s')} reactivity [eV]"
    for name in model.descriptor_names
]
mm.subplots_adjust_kwargs = {'left': 0.2, 'right': 0.8, 'bottom': 0.15}
mm.plot(save=str(output_dir / 'rate_control.pdf'))
_t_mark('DRC分析')


# ─────────────────────────────────────────────
#  数据导出 & 归档
# ─────────────────────────────────────────────
export_tof_and_coverage(model, output_dir)
_t_mark('数据导出')

save_run_to_archive(output_dir)
_t_mark('归档保存')

# ─────────────────────────────────────────────
#  运行时间统计（写入 .log 文件）
# ─────────────────────────────────────────────
_t_total_end = time.perf_counter()
_t_log_file = output_dir / f"{Path(mkm_file).stem}.log"

def _t_fmt(sec):
    h, rem = divmod(int(sec), 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}" if h > 0 else f"{m:02d}:{s:02d}"

with open(_t_log_file, 'a', encoding='utf-8') as _f:
    _f.write('\n' + '═' * 60 + '\n')
    _f.write('  运行时间统计\n')
    _f.write('═' * 60 + '\n')
    _f.write(f'  开始时间: {datetime.fromtimestamp(_t_total_start).strftime("%Y-%m-%d %H:%M:%S")}\n')
    _f.write(f'  结束时间: {datetime.fromtimestamp(_t_total_end).strftime("%Y-%m-%d %H:%M:%S")}\n')
    _f.write(f'  总耗时:   {_t_fmt(_t_total_end - _t_total_start)} ({_t_total_end - _t_total_start:.2f} 秒)\n')
    _f.write('─' * 60 + '\n')
    _f.write('  各阶段耗时:\n')
    for _i in range(1, len(_t_checkpoints)):
        _label, _t0 = _t_checkpoints[_i - 1]
        _,    _t1 = _t_checkpoints[_i]
        _dt = _t1 - _t0
        _f.write(f'    {_label:<22s}  {_t_fmt(_dt):>8s}  ({_dt:>7.2f} 秒)\n')
    _f.write('═' * 60 + '\n')

print(f"\n运行时间统计已追加到: {_t_log_file}")
print(f"总耗时: {_t_fmt(_t_total_end - _t_total_start)} ({_t_total_end - _t_total_start:.2f} 秒)")
for _i in range(1, len(_t_checkpoints)):
    _label, _t0 = _t_checkpoints[_i - 1]
    _,    _t1 = _t_checkpoints[_i]
    print(f"  {_label:<22s}  {_t1 - _t0:>7.2f} 秒")

