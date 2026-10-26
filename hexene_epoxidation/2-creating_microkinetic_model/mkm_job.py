from catmap import ReactionModel
from pathlib import Path

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

model.run()

from catmap import analyze
descriptor_names = ['OOH_s', 'O_s']
vm = analyze.VectorMap(model)
vm.plot_variable = 'rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-25 #minimum rate to plot
vm.max = 1e3 #maximum rate to plot
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

vm.include_labels = ['OOH_s']
vm.plot(save=str(output_dir /'OOH_coverage.pdf'))

sa = analyze.ScalingAnalysis(model)
sa.plot(save=str(output_dir /'scaling.pdf'))
