from pathlib import Path

from catmap import ReactionModel
from descriptor_output import configure_descriptor_output
from tof_export import export_tof_tables

mkm_file = 'alkene_epoxidation.mkm'
model = ReactionModel(setup_file=mkm_file)
# CatMAP writes its native model-state log at the end of model.run().
# Keep that artifact separate so the CMD transcript can own alkene_epoxidation.log.
model.log_file = str(Path(__file__).resolve().parent / 'alkene_epoxidation_catmap.log')
output_dir = configure_descriptor_output(model, __file__)

scaler_variables = ['rxn_parameter','electronic_energy',
    'free_energy','zero_point_energy','enthalpy','entropy']
# Warning: adding 'rate_control' or 'selectivity_control' to output_variables increases calculation time significantly
solver_variables = ['production_rate','rate_control','coverage','rate',
    'consumption_rate','selectivity','rxn_direction',
    'rate_constant','equilibrium_constant','carbon_selectivity',
    'selectivity_control','rxn_order','apparent_activation_energy',
    'interacting_energy','directional_rates','forward_rate','reverse_rate',
    'forward_rate_constant','reverse_rate_constant']
two_dim_variables = ['rate_control', 'selectivity_control', 'rxn_order']
one_dim_variables = list(set(scaler_variables + solver_variables) - set(two_dim_variables))
model.output_variables += solver_variables + scaler_variables

model.run()

export_tof_tables(
    model,
    output_dir,
    product='C6H12O_g'
)

from catmap import analyze
vm = analyze.VectorMap(model)
vm.log_scale = True  # not necessarily the right choice of parameters for all output_variables
vm.min = 1e-200
vm.max = 1e-1
vm.threshold = 1e-300
vm.unique_only = False

# vm.n_ticks = 10 # 增加全局色度条刻度密度

for out_var in one_dim_variables:
    vm.plot_variable = out_var
    vm.n_ticks = 4 if out_var == 'directional_rates' else 6 # 只针对该文件增加全局色度条刻度密度
    fig = vm.plot(save=str(output_dir / (out_var + '.pdf')))
    fig.clf()

mm = analyze.MatrixMap(model)
mm.log_scale = False
mm.min = -2
mm.max = 2
mm.unique_only = False
for out_var in two_dim_variables:
    mm.plot_variable = out_var
    fig = mm.plot(save=str(output_dir / (out_var + '.pdf')))
    fig.clf()
