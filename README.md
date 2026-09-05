# Sequencing read depth calculator

根据 FASTQ 文件中的 reads 数量估算理论测序深度。脚本支持单端、双端数据，以及未压缩或 gzip 压缩的 FASTQ 文件。

> 本工具计算的是原始测序数据的理论深度，不是比对后的覆盖深度。若需要统计 BAM 文件在染色体或区间上的实际覆盖深度，请参考 [plot_depth.md](plot_depth.md)。

## 计算公式

```text
depth = read_count x end_count x read_length / genome_size
```

- `read_count`：单个输入 FASTQ 中的记录数；每条 FASTQ 记录必须占 4 行。
- `end_count`：单端测序为 `1`，双端测序为 `2`。
- `read_length`：每个 end 的读长，单位为 bp。
- `genome_size`：单倍体基因组大小，单位为 bp。

对于双端数据，请为每个样本只传入一个 mate 文件（通常为 R1），因为 `end_count=2` 已经将 R1 和 R2 都计入公式。若同时传入同一样本的 R1 和 R2，脚本会将它们作为两个独立输入分别输出，而不是合并为一个样本。

## 依赖

- Bash 3.2 或更高版本
- `awk`
- `gzip`（读取 `.gz` 文件时需要）

## 使用方法

```bash
./read_depth_calcu.sh [options] END_COUNT READ_LENGTH GENOME_SIZE [FASTQ ...]
```

首次使用时赋予脚本执行权限：

```bash
chmod +x read_depth_calcu.sh
```

### 单端测序

以下示例计算两个 150 bp 单端文库相对于 dm6 基因组的理论深度：

```bash
./read_depth_calcu.sh 1 150 146638899 sample1.fastq.gz sample2.fastq.gz
```

### 双端测序

双端测序只传入每个样本的 R1（或只传 R2）：

```bash
./read_depth_calcu.sh 2 150 146638899 sample1_R1.fastq.gz sample2_R1.fastq.gz
```

### 自动发现输入文件

省略输入文件时，脚本会扫描当前目录中的 `*.fastq.gz`、`*.fq.gz`、`*.fastq` 和 `*.fq`：

```bash
./read_depth_calcu.sh 1 150 146638899
```

双端数据建议显式列出 R1 文件，避免 R1、R2 分别产生重复的样本结果。

### 指定输出位置

默认结果写入 `read_depth.tsv`。使用 `-o` 修改路径，或用 `-o -` 输出到标准输出：

```bash
./read_depth_calcu.sh -o results/depth.tsv 1 150 146638899 sample.fastq.gz
./read_depth_calcu.sh -o - 1 150 146638899 sample.fastq.gz
```

查看完整帮助：

```bash
./read_depth_calcu.sh --help
```

## 输出格式

输出为制表符分隔的 TSV 文件，每个输入文件对应一行：

| 列名              | 含义                        |
| ----------------- | --------------------------- |
| `file`            | 输入 FASTQ 路径             |
| `read_count`      | FASTQ 记录数                |
| `end_count`       | 单端为 1，双端为 2          |
| `read_length_bp`  | 每个 end 的读长             |
| `sequenced_bases` | 估算的测序碱基总数          |
| `genome_size_bp`  | 单倍体基因组大小            |
| `depth_x`         | 理论测序深度，保留 6 位小数 |

示例：

```text
file    read_count  end_count  read_length_bp  sequenced_bases  genome_size_bp  depth_x
sample.fastq.gz  1000000  1  150  150000000  146638899  1.022920
```

## 输入校验与注意事项

- 参数必须为正整数，`END_COUNT` 只能为 `1` 或 `2`。
- FASTQ 总行数必须能被 4 整除；损坏或不完整的文件会直接报错。
- 脚本只统计 reads 数量，不执行接头去除、质量过滤或比对。
- 变长 reads 数据可将质控后的平均读长作为 `READ_LENGTH`，但结果仍是估算值。
- 输出路径不能与任一输入 FASTQ 相同，脚本会拒绝覆盖原始数据。
- 每个结果先写入临时文件，全部输入处理成功后才替换目标文件，避免失败时留下不完整结果。

## 测试

```bash
bash tests/test_read_depth_calcu.sh
```

## License

[MIT](LICENSE)
