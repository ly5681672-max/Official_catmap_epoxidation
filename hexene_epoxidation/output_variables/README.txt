CatMAP output_variables：TOF数值表最小增量
====================================

现有 rate_control.pdf 就是 DRC（速率控制程度）图，无需重跑另一脚本。
现有 production_rate.pdf 不是数值表；vm.threshold=1e-25 会隐藏低于阈值的速率曲线。

安装：
1. 将 tof_export.py 放在 hexene_epoxidation/output_variables/ 目录内。
2. 在该目录的 mkm_job.py 顶部添加：
   from tof_export import export_tof_tables
3. 在 mkm_job.py 末尾（绘图循环之后）添加：
   export_tof_tables(model, output_dir, product='C6H12O_g')
4. 正常运行已有 run_mkm_utf8.cmd（无需改动CMD）。

预期文件：
  descriptor_C6H12_OOH/TOF_C6H12O_g_grid.csv
  descriptor_C6H12_OOH/TOF_C6H12O_g_catalysts.csv

含义：
- grid.csv：当前描述符网格上CatMAP求解的产物生成速率；并非各催化剂真实点的结果。
- catalysts.csv：选定14个催化剂在各自真实描述符值处进行的CatMAP单点求解结果；不是最近邻取值。
- 均取CatMAP production_rate中的环氧化产物C6H12O_g，表示正向产物生成速率的非负部分；不是把单步反应速率误当TOF。
- 当前模型的TOF是按CatMAP有效位点归一化的模型值，不能直接当作实验催化剂质量归一化速率。
- 所有点必须求解成功后才写出CSV。若有任何点求解失败，明确报错，不补零、不插值。

执行状态：仅完成语法及模拟测试；尚未执行用户本地的CatMAP实算，也未修改GitHub仓库。
