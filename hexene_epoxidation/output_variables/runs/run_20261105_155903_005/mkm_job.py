from pathlib import Path
from math import isnan
import numpy as np

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
vm.threshold = 0
vm.unique_only = False

# vm.n_ticks = 10 # 增加全局色度条刻度密度

for out_var in one_dim_variables:
    vm.plot_variable = out_var
    print(f'PLOT_VARIABLE={out_var}', flush=True)
    raw_map = getattr(model, out_var + '_map', None)
    if raw_map is None or len(raw_map) == 0:
        print(f'PLOT_DATA_STATUS=SKIPPED_EMPTY_MAP variable={out_var} points=0', flush=True)
        continue

    finite_values = []
    nan_count = posinf_count = neginf_count = conversion_error_count = 0
    raw_count = 0
    for _point, output_values in raw_map:
        for raw_value in output_values:
            raw_count += 1
            try:
                value = float(raw_value)
            except (TypeError, ValueError, OverflowError):
                conversion_error_count += 1
                continue
            if isnan(value):
                nan_count += 1
            elif value == float('inf'):
                posinf_count += 1
            elif value == float('-inf'):
                neginf_count += 1
            else:
                finite_values.append(value)

    invalid_count = nan_count + posinf_count + neginf_count + conversion_error_count
    if not finite_values:
        print(
            f'PLOT_DATA_STATUS=SKIPPED_ALL_DATA_INVALID variable={out_var} total={raw_count} '
            f'invalid={invalid_count} NaN={nan_count} +Inf={posinf_count} '
            f'-Inf={neginf_count} conversion_error={conversion_error_count}; raw map unchanged',
            flush=True,
        )
        continue
    if invalid_count:
        print(
            f'PLOT_DATA_STATUS=SKIPPED_NONFINITE variable={out_var} total={raw_count} '
            f'finite={len(finite_values)} invalid={invalid_count} NaN={nan_count} +Inf={posinf_count} '
            f'-Inf={neginf_count} conversion_error={conversion_error_count}; raw map unchanged',
            flush=True,
        )
        continue

    data_min = min(finite_values)
    data_max = max(finite_values)
    zero_count = sum(value == 0.0 for value in finite_values)
    negative_count = sum(value < 0.0 for value in finite_values)
    positive_count = sum(value > 0.0 for value in finite_values)

    # Raw positivity is necessary for log plotting. CatMAP uses cubic interpolation
    # before its final log10; cubic overshoot can still create nonpositive grid values.
    use_log_scale = data_min > 0.0 and data_max > data_min
    interpolation_nonfinite_count = 0
    interpolation_nonpositive_count = 0
    if use_log_scale:
        log_interpolated_map = model.map_to_array(
            raw_map,
            model.descriptor_ranges,
            vm.resolution * vm.resolution_enhancement,
            log_interpolate=True,
            minval=data_min,
        )
        interpolation_nonfinite_count = int(np.count_nonzero(~np.isfinite(log_interpolated_map)))
        interpolation_nonpositive_count = int(np.count_nonzero(log_interpolated_map <= 0.0))
        if interpolation_nonfinite_count or interpolation_nonpositive_count:
            use_log_scale = False
            scale_reason = 'catmap_log_interpolation_nonfinite_or_nonpositive'
        else:
            scale_reason = 'raw_and_catmap_log_interpolation_strictly_positive'
    else:
        scale_reason = 'raw_data_not_strictly_positive_or_constant'

    vm.log_scale = use_log_scale
    vm.min = data_min
    vm.max = data_max
    scale_name = 'log' if use_log_scale else 'linear'
    print(
        f'PLOT_DATA_STATUS=READY variable={out_var} scale={scale_name} reason={scale_reason} '
        f'total={raw_count} finite={len(finite_values)} positive={positive_count} '
        f'zero={zero_count} negative={negative_count} min={data_min:.17g} max={data_max:.17g} '
        f'interpolation_nonfinite={interpolation_nonfinite_count} '
        f'interpolation_nonpositive={interpolation_nonpositive_count}',
        flush=True,
    )

    if data_min == data_max:
        print(
            f'PLOT_DATA_STATUS=SKIPPED_CONSTANT variable={out_var} value={data_min:.17g}; '
            'no contour levels generated',
            flush=True,
        )
        continue

    vm.n_ticks = 4 if out_var == 'directional_rates' else 6 # 只针对该文件增加全局色度条刻度密度
    try:
        fig = vm.plot(save=str(output_dir / (out_var + '.pdf')))
    except Exception as exc:
        print(
            f'PLOT_STATUS=FAILED variable={out_var} error={type(exc).__name__}: {exc}',
            flush=True,
        )
        raise
    else:
        fig.clf()
        print(f'PLOT_STATUS=SUCCESS variable={out_var}', flush=True)

mm = analyze.MatrixMap(model)
mm.log_scale = False
mm.min = -2
mm.max = 2
mm.unique_only = False
for out_var in two_dim_variables:
    mm.plot_variable = out_var
    fig = mm.plot(save=str(output_dir / (out_var + '.pdf')))
    fig.clf()
