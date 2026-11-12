
import pandas as pd
from pathlib import Path

base = Path(
    r"H:\catmap\Official_catmap_epoxidation"
    r"\hexene_epoxidation\2-creating_microkinetic_model"
    r"\runs\run_20261109_113149_536"
)

for filename in ["error_log.csv", "jacobian_norm.csv"]:

    path = base / filename

    print("\n==========", filename, "==========")

    # 关键：明确告诉 Pandas 文件没有表头
    df = pd.read_csv(
        path,
        header=None,
        dtype=str,
        nrows=10
    )

    print("实际列数：", df.shape[1])
    print("前10条原始记录：")
    print(df.to_string(index=False, header=False))

    print("文件大小：",
          round(path.stat().st_size / 1024**2, 2), "MB")
