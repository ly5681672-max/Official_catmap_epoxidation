
import csv
from pathlib import Path
from collections import defaultdict
from decimal import Decimal, InvalidOperation
import pandas as pd

BASE = Path(
    r"H:\catmap\Official_catmap_epoxidation"
    r"\hexene_epoxidation\2-creating_microkinetic_model"
    r"\runs\run_20261109_113149_536"
)

OUT = BASE / "diagnostics"
OUT.mkdir(exist_ok=True)


def scan_csv(filename, kind):
    path = BASE / filename

    # key = (Ha, O)
    # 每个描述符点保留多次独立求解记录
    records = defaultdict(list)

    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)

        for row in reader:
            if not row:
                continue

            row = [x.strip() for x in row]

            expected = 4 if kind == "error" else 5
            if len(row) != expected:
                continue

            try:
                ha = row[0]
                oxygen = row[1]
                iteration = int(row[2])

                values = [Decimal(x) for x in row[3:]]

            except (ValueError, InvalidOperation):
                continue

            key = (ha, oxygen)

            if not records[key] or iteration <= records[key][-1]["iteration"]:
                records[key].append({
                    "iteration": iteration,
                    "rows": 0,
                    "last": None
                })

            attempt = records[key][-1]
            attempt["iteration"] = iteration
            attempt["rows"] += 1
            attempt["last"] = values

    output = []

    for (ha, oxygen), attempts in records.items():

        for number, attempt in enumerate(attempts, start=1):

            values = attempt["last"]

            item = {
                "Ha": ha,
                "O": oxygen,
                "attempt": number,
                "record_count": attempt["rows"],
                "last_iteration": attempt["iteration"]
            }

            if kind == "error":
                item["final_residual"] = str(values[0])

            else:
                item["jacobian_col4"] = str(values[0])
                item["jacobian_col5"] = str(values[1])

            output.append(item)

    result = pd.DataFrame(output)

    result.to_csv(
        OUT / f"{kind}_attempt_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    print(f"\n========== {filename} ==========")
    print("描述符点数量:", len(records))
    print("求解批次数量:", len(output))
    print("平均每点求解次数:", round(len(output) / max(len(records), 1), 2))
    print("输出文件:", OUT / f"{kind}_attempt_summary.csv")

    return result


error = scan_csv("error_log.csv", "error")
jacobian = scan_csv("jacobian_norm.csv", "jacobian")

# 输出残差最大的求解批次
error["residual_float"] = pd.to_numeric(
    error["final_residual"],
    errors="coerce"
)

error_sorted = error.sort_values(
    "residual_float",
    ascending=False
)

error_sorted.head(100).to_csv(
    OUT / "largest_residual_attempts.csv",
    index=False,
    encoding="utf-8-sig"
)

# 输出每个描述符点最后一次求解记录
final_error = (
    error.sort_values("attempt")
    .groupby(["Ha", "O"], as_index=False)
    .tail(1)
)

final_error.to_csv(
    OUT / "final_error_by_grid.csv",
    index=False,
    encoding="utf-8-sig"
)

print("\n诊断完成。")
print("结果目录:", OUT)
print("\n残差最大的10条求解记录：")
print(error_sorted.head(10).to_string(index=False))
