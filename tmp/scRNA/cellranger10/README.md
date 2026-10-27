# Ceratopteris richardii：Cell Ranger 10 计数

独立于历史模板 `../celcount.sh`；不修改旧流程。适用于分别建库的 C1/C2 常规 10x GEX 数据，不将两个样本合并计数。这里提供可复用命令，不表示 2026-10-16 的实际任务已经完成。

## 与旧模板的差异

| 项目 | 历史模板 | 本目录脚本 |
| --- | --- | --- |
| 软件 | 注释为 7.0.0，实际调用 PATH 中的 `cellranger` | 默认 `cellranger10`，检查实际版本为 10.x |
| 核数 | `--localcores=30` | 默认 20，可由 `LOCALCORES` 覆盖 |
| BAM | 未显式指定 | `--create-bam=true` |
| 二次分析 | `--nosecondary` | 保留；不做聚类/降维，不代表跳过比对或计数 |
| 输入与失败 | 参数未引用，无参数检查 | 引用路径、检查参数、传播非零退出码 |
| 日志 | 示例为 `2>&1 > file` | `> file 2>&1`，同时保存 stdout/stderr |

Cell Ranger 8 起必须显式指定 `--create-bam`；旧命令不能原样搬到 10。`--create-bam true` 和本脚本的 `--create-bam=true` 等效。BAM 会增加磁盘和计算开销。`--localcores` 是本地任务的核数上限，不能代替内存限制；需要时设置 `LOCALMEM`（正整数 GB），并匹配调度器分配。

`cellranger10` 是本地命名约定；官方命令通常为 `cellranger`。脚本需要 PATH 中的可执行文件或完整路径，不能依赖交互式 shell alias。例如 `export CELLRANGER_BIN=/opt/cellranger-10.0.0/cellranger`。启动前可运行 `"${CELLRANGER_BIN:-cellranger10}" --version` 和 `"${CELLRANGER_BIN:-cellranger10}" count --help` 核对实际安装。

参考目录必须为兼容的 `cellranger mkref` 输出，不能直接传 FASTA/GTF。确认物种、组装和注释版本；如旧 STAR 索引不兼容，应使用 Cell Ranger 10 在新目录重建参考，保留原参考。脚本仅检查目录存在，完整参考及 FASTQ 兼容性由 Cell Ranger preflight 检查。`--sample` 是 FASTQ 文件名前缀，不是任意标签；确认 C1/C2 的 R1/R2 文件命名和前缀。保留自动 chemistry 检测，不根据物种猜试剂版本；10x 3′ v4 的试剂版本与软件 10.x 版本是不同概念。未新增 `--force-cells` 或更改内含子计数策略（10 默认包含内含子），跨版本比较还需核对算法与参考差异。

## C1/C2 串行运行

从 Ceratopteris 项目目录执行；将脚本路径换成本地仓库的绝对路径。输入目录应为 `dataVault/rawdata/C1/` 与 `dataVault/rawdata/C2/`。

```bash
RUNNER=/absolute/path/NGS-pipeline/tmp/scRNA/cellranger10/run_c1_c2.sh
export CELLRANGER_BIN=cellranger10
export LOCALCORES=20
# export LOCALMEM=120  # 示例；按实际分配设置，单位 GB
mkdir -p results/cellranger10
launch_log="results/cellranger10/launch-$(date +%Y%m%dT%H%M%S)-$$.log"
nohup bash "$RUNNER" \
  /data/wmc_data/index/cellrangerRNA/Ceratopteris_richardii \
  ./dataVault/rawdata \
  ./results/cellranger10 > "$launch_log" 2>&1 &
pid=$!
printf '%s\n' "$pid" > "${launch_log%.log}.pid"
printf 'PID=%s; launch log=%s\n' "$pid" "$launch_log"
```

C1 命令退出码为 0 后才启动 C2；C1 失败则停止，C2 不会运行。结果分别在 `results/cellranger10/C1/`、`C2/`。每次调用都会创建唯一的 `logs/run-*/`，包含 `C1.log`、`C2.log`（仅在启动后存在）与追加式 `status.tsv`；单样本日志记录版本、参数和计数命令退出码。启动日志记录本次具体日志目录。

单独计数也可复用（输出在当前目录下的 ID 文件夹）：

```bash
bash /absolute/path/NGS-pipeline/tmp/scRNA/cellranger10/count.sh \
  C1 C1 /absolute/path/reference /absolute/path/rawdata/C1 > C1.log 2>&1
```

## 状态与结果检查

```bash
# 在上述启动 shell 中；新 shell 请填入实际日志路径和 PID
tail -n 40 "$launch_log"
ps -p "$pid" -o pid,etime,args
# 将 RUN_LOG_DIR 替换为启动日志中打印的 Logs 路径
RUN_LOG_DIR=/absolute/path/results/cellranger10/logs/run-TIMESTAMP-SUFFIX
cat "$RUN_LOG_DIR/status.tsv"
tail -n 50 "$RUN_LOG_DIR/C1.log"
tail -n 50 "$RUN_LOG_DIR/C2.log"  # C2 尚未启动时文件不存在

# 检查各样本的运行记录及主要产物
for sample in C1 C2; do
  tail -n 30 "results/cellranger10/$sample/_log"
  ls -lh "results/cellranger10/$sample/outs/"{web_summary.html,metrics_summary.csv,filtered_feature_bc_matrix.h5,possorted_genome_bam.bam}
done
```

`SUCCEEDED` 表示本次命令退出码为 0，还需结合 `_log` 中的完成信息、主要输出及 `web_summary.html` 的细胞数/比对率等 QC 判断。进程消失、日志非空或 `outs/` 存在都不足以证明成功。`status.tsv` 是事件记录，以每个样本的最后一行为准；被强制终止或断电时可能停在 `RUNNING`，需结合进程和 Cell Ranger 日志判断，不能据此认定仍在运行。`ps` 仅为辅助检查，历史 PID 可能被复用。

失败后先查看该样本日志和 pipestance 中的错误记录，修正原因。同一输出目录、同一 ID 和相同参数重新运行时，Cell Ranger 可尝试恢复兼容的未完成任务；先确认旧进程已退出，不要并发复用同一 ID。更改参考或参数应改用新输出目录。若 C1 已完成而仅需重试 C2，使用单样本脚本即可，不必重跑整个串行入口。不要为重跑直接删除已有输出。

## 官方参数依据

- [Cell Ranger 10 命令行参数](https://www.10xgenomics.com/support/software/cell-ranger/10.0/resources/cr-command-line-arguments)
- [Cell Ranger count 与 create-bam 迁移说明](https://www.10xgenomics.com/support/software/cell-ranger/9.0/tutorials/cr-tutorial-ct)

本次验证范围：Bash 语法与模拟可执行程序的参数传递、路径空格、串行/失败行为；未在真实 C1/C2 数据或真实 Cell Ranger 10 上运行。
