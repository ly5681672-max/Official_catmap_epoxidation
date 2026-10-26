from catmap import ReactionModel
from pathlib import Path
from math import isfinite


def auto_descriptor_ranges(model, margin=0.05):
    ranges = []

    for name in model.descriptor_names:
        data = model.species_definitions[name].get(
            "formation_energy"
        )

        if not isinstance(data, (list, tuple)):
            raise ValueError(f"{name}: 缺少逐催化剂能量数据")

        if len(data) != len(model.surface_names):
            raise ValueError(f"{name}: 能量数量与催化剂数量不一致")

        values = []
        for surface, energy in zip(model.surface_names, data):
            try:
                value = float(energy)
            except (ValueError, TypeError):
                raise ValueError(
                    f"{surface}: {name} 缺少有效能量"
                )

            if not isfinite(value):
                raise ValueError(
                    f"{surface}: {name} 能量不是有限数值"
                )

            values.append(value)

        low = min(values)
        high = max(values)

        if high <= low:
            raise ValueError(
                f"{name}: 数据没有有效跨度，无法自动确定范围"
            )

        padding = (high - low) * margin
        ranges.append([
            low - padding,
            high + padding
        ])

    return ranges


mkm_file = 'alkene_epoxidation.mkm'
model = ReactionModel(setup_file=mkm_file)

# 根据当前描述符自动建立结果文件夹
descriptor_pair = "_".join(
    name.removesuffix("_s")
    for name in model.descriptor_names
)

output_dir = (
    Path(__file__).resolve().parent
    / f"descriptor_{descriptor_pair}"
)

output_dir.mkdir(parents=True, exist_ok=True)

print(f"当前描述符: {model.descriptor_names}")
print(f"图表输出目录: {output_dir}")

model.descriptor_ranges = auto_descriptor_ranges(model)

print("自动描述符范围:")
for name, bounds in zip(
    model.descriptor_names, model.descriptor_ranges
):
    print(f"  {name}: {bounds}")

model.run()

from catmap import analyze
vm = analyze.VectorMap(model)
vm.plot_variable = 'rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-80 #minimum rate to plot
vm.max = 1e-1 #maximum rate to plot
vm.plot(save=str(output_dir /'rate.pdf')) #draw the plot and save it as "rate.pdf"

vm.unique_only = False
vm.plot(save=str(output_dir /'all_rates.pdf'))
vm.unique_only = True

model.output_variables += ['production_rate']
model.run()
vm.production_rate_map = model.production_rate_map #attach map
vm.threshold = 1e-30 #do not plot rates below this
vm.plot_variable = 'production_rate'
vm.plot(save=str(output_dir /'production_rate.pdf'))

vm.descriptor_labels = [
    f"{name.removesuffix('_s')} reactivity [eV]"
    for name in model.descriptor_names
]

vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save=str(output_dir /'pretty_production_rate.pdf'))

vm.plot_variable = 'coverage'
vm.log_scale = False
vm.min = 0
vm.max = 1
vm.plot(save=str(output_dir /'coverage.pdf'))

coverage_descriptor = model.descriptor_names[0]

vm.include_labels = [coverage_descriptor]
vm.plot(
    save=str(
        output_dir /
        f"{coverage_descriptor.removesuffix('_s')}_coverage.pdf"
    )
)

sa = analyze.ScalingAnalysis(model)
sa.plot(save=str(output_dir /'scaling.pdf'))