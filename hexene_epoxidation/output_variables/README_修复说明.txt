CatMAP 14位点TOF/DRC导出修复（2026-10-27）

检查的最新提交：9154be680873dd441e8dbca4bd591fe40ad9436d
运行日志：ValueError: tife: unexpected coverage or step-rate dimensions（tof_export.py:242）
根因：numbers solver 启用时，CatMAP _coverage 包含 adsorbate_names 对应覆盖度 + 末尾显式空位覆盖度；output_labels[coverage] 只有吸附物种名。旧代码要求二者长度相等，错误地在首个位点报错。

修改：
1. 根据 model.use_numbers_solver 选择期望维度。
2. 只把前 len(output_labels[coverage]) 项作为吸附物种覆盖度。
3. 当存在显式空位时单独提取 theta_star_solver，核对 site_balance_residual。
4. 保留原先导出、DRC、通量核验、回归和图形逻辑；无物理参数更改。

使用：
1. 先备份 H:\catmap\Official_catmap_epoxidation\hexene_epoxidation\output_variables\tof_export.py
2. 解压 ZIP，把其中的 tof_export.py 覆盖该文件。
3. 用相同的 run_mkm_utf8.cmd 运行。
4. 检查 CatMAP execution exit code 和新的 OOH_kinetic_evidence_C6H12O_g_14sites.csv。
5. DRC 若出现新错误，请保留最新 alkene_epoxidation.log。

验证：Python语法与14位点模拟导出测试通过（numbers solver=True/False）；非用户本机实际CatMAP求解验证。
