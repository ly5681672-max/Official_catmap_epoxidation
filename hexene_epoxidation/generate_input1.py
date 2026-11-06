# 生成 CatMAP 输入文件：计算形成能、输出 TXT/CSV，并验证解析结果。
from ase.symbols import string2symbols
import csv
import os
from pathlib import Path

abinitio_energies = {
    # ========== 气相物种 ==========
    'H2O2_gas': -18.122091,
    'C6H12_gas': -95.031223,
    'C6H12O_gas': -101.037771,
    'H2O_gas': -14.204457,

    # ========== 吸附物种 ==========
    'Ha_111':     -868.06986,
    'Hb_111':     -867.91225,
    'OOH_111':    -877.1433925,
    'C6H12_111':  -958.3287775,
    'C6H12O_111': -964.45279,
    'O_111':      -866.997122,
    'OH_111':     -872.6800805,
    'H2O_111':    -877.747855,

    # ========== 过渡态 ==========
    'OOH-C6H12_111': -983.369042,
    'Hb-O_111':      -872.064286,

    # ========== 裸板 ==========
    'slab_111': -863.139429,
}
# 预计算原子参考能，单位：eV
# ref_dict = {
#     'O': -4.37,
#     'C': -9.28,
#     'H': -1.11,
#     '111': abinitio_energies['slab_111'],
# }
#该方法原子参考能对线性关系描述更合理
ref_dict = {}
ref_dict['O'] = abinitio_energies['H2O2_gas']-abinitio_energies['H2O_gas']
ref_dict['H'] = 0.5*(abinitio_energies['H2O_gas'] - ref_dict['O'])
ref_dict['C'] = 1/6*(abinitio_energies['C6H12_gas'] - 12*ref_dict['H'])
ref_dict['111'] = abinitio_energies['slab_111']

# 特殊物种名称 -> 实际化学组成。
# Ha/Hb 表示不同吸附位置上的 H，不是化学元素 Ha/Hb。
formula_map = {
    'Ha': 'H',
    'Hb': 'H',
    'Hb-O': 'HO',
}


def get_formation_energies(energy_dict, references):
    """根据原子参考能和裸板能计算形成能。"""
    formation_energies = {}

    for key, energy in energy_dict.items():
        if '_' not in key:
            raise ValueError(f'物种 key 缺少 site 信息: {key}')

        name, site = key.split('_', 1)

        if 'slab' in name:
            continue

        E0 = energy

        # 表面吸附态/过渡态先减去裸板能；气相不减裸板。
        if site == '111':
            E0 -= references['111']

        # 特殊名称显式映射，普通名称则移除 TS 连字符后按化学式解析。
        formula = formula_map.get(name, name.replace('-', ''))

        try:
            composition = string2symbols(formula)
        except Exception as exc:
            raise ValueError(
                f'无法解析物种 {name}，对应化学式 {formula}，原 key: {key}'
            ) from exc

        for atom in composition:
            if atom not in references:
                raise KeyError(f'物种 {key} 中的原子 {atom} 未在 ref_dict 中定义参考能')
            E0 -= references[atom]

        formation_energies[key] = round(E0, 3)

    return formation_energies


formation_energies = get_formation_energies(abinitio_energies, ref_dict)

# 打印形成能检查
for key, value in formation_energies.items():
    print(f'{key} {value}')

frequency_dict = {
    'H2O2_gas': [],
    'C6H12_gas': [],
    'C6H12O_gas': [],
    'H2O_gas': [],
    'Ha_111': [],
    'Hb_111': [],
    'OOH_111': [],
    'C6H12O_111': [],
    'O_111': [],
    'OH_111': [],
    'OOH-C6H12_111': [],
    'Hb-O_111': [],
    'slab_111': [],
}


def make_input_file(file_name, energy_dict, frequencies):
    """生成 CatMAP TableParser 使用的制表符分隔输入文件。"""
    header = '\t'.join([
        'surface_name',
        'site_name',
        'species_name',
        'formation_energy',
        'frequencies',
        'reference',
    ])

    lines = []

    for key, energy in energy_dict.items():
        if '_' not in key:
            continue

        name, site = key.split('_', 1)

        if 'slab' in name:
            continue

        frequency = frequencies.get(key, [])
        surface = None if site == 'gas' else 'ti'
        outline = [surface, site, name, energy, frequency, 'Input File Tutorial.']
        lines.append('\t'.join(str(value) for value in outline))

    lines.sort()
    input_file = '\n'.join([header] + lines)

    with open(file_name, 'w', encoding='utf-8') as file:
        file.write(input_file)

    print(f'Successfully created input file: {file_name}')


def make_csv_file(file_name, energy_dict, frequencies):
    """将同一套形成能数据输出为 CSV，便于表格软件检查。"""
    fieldnames = [
        'surface_name',
        'site_name',
        'species_name',
        'formation_energy',
        'frequencies',
        'reference',
    ]

    with open(file_name, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for key, energy in energy_dict.items():
            if '_' not in key:
                continue

            name, site = key.split('_', 1)

            if 'slab' in name:
                continue

            freq_list = frequencies.get(key, [])
            freq_str = ';'.join(
                str(int(freq))
                if isinstance(freq, (int, float)) and float(freq).is_integer()
                else str(freq)
                for freq in freq_list
            ) if freq_list else ''

            writer.writerow({
                'surface_name': None if site == 'gas' else 'tini',
                'site_name': site,
                'species_name': name,
                'formation_energy': energy,
                'frequencies': freq_str,
                'reference': 'Input File Tutorial.',
            })

    print(f'Successfully created CSV file: {file_name}')


# Prevent silently overwriting reviewed energy tables. Opt in explicitly if needed.
def assert_safe_output(paths):
    existing = [str(p) for p in paths if Path(p).exists()]
    if existing and os.environ.get('ALLOW_ENERGY_OVERWRITE') != '1':
        raise FileExistsError(
            'Refusing to overwrite existing energy tables: ' + ', '.join(existing)
            + '. Set ALLOW_ENERGY_OVERWRITE=1 only after reviewing the source values.'
        )


# 生成输出文件
txt_file_name = 'tini_energies.txt'
csv_file_name = 'tini_energies.csv'
assert_safe_output([txt_file_name, csv_file_name])
make_input_file(txt_file_name, formation_energies, frequency_dict)
make_csv_file(csv_file_name, formation_energies, frequency_dict)

# CatMAP 解析测试
try:
    from catmap.model import ReactionModel
    from catmap.parsers import TableParser

    rxm = ReactionModel()
    rxm.surface_names = ['tini']

    # CatMAP 内部物种 key 必须带 site 后缀。
    # 输入表中的 species_name 仍保持 H2O2/Ha/Hb/...，TableParser 会根据 site_name=111
    # 将这些名称匹配到 *_s 物种。
    rxm.adsorbate_names = (
         'Ha_s', 'Hb_s', 'OOH_s',
        'C6H12O_s', 'O_s', 'OH_s',
    )
    rxm.transition_state_names = ('OOH-C6H12_s', 'Hb-O_s')
    rxm.gas_names = ('H2O2_g', 'C6H12_g', 'C6H12O_g', 'H2O_g')
    rxm.site_names = ('s',)

    # Ha/Hb/Hb-O 不是标准化学式，必须在于模型一致的 *_s key 上显式指定组成。
    rxm.species_definitions = {
        's': {'site_names': ['111']},
        'Ha_s': {'composition': {'H': 1}},
        'Hb_s': {'composition': {'H': 1}},
        'Hb-O_s': {'composition': {'H': 1, 'O': 1}},
    }

    parser = TableParser(rxm)
    parser.input_file = txt_file_name
    parser.parse()

    print('CatMAP parsing successful.')
    for key, value in rxm.species_definitions.items():
        print(f'{key} {value}')

except Exception as exc:
    print('CatMAP 解析测试失败：')
    print(exc)
    # 不再吞掉异常：解析失败时让进程返回非 0 退出码，便于定位问题。
    raise
