[![English](https://img.shields.io/badge/Language-English-2563eb)](README.md)
[![中文](https://img.shields.io/badge/语言-中文-0f766e)](README_zh-CN.md)

# Flora 私有开发仓库

本仓库保存 Flora 私有源码、测试和打包工具，必须保持 private。
公开仓库 `brilliantlee2/Flora` 保存用户文档，二进制压缩包通过 GitHub Releases 发布。

## 克隆私有仓库

```bash
git clone --branch codex/1028-ont-v330-performance --single-branch \
  https://github.com/brilliantlee2/flora_private.git
cd flora_private
```

需要私有仓库访问权限，使用 GitHub Token 或已授权的 SSH key。
SSH 地址为 `git@github.com:brilliantlee2/flora_private.git`。
修改过的依赖已包含在 `vendor/` 中，不需要配置 submodule。

## 构建环境

在 Linux x86_64 的仓库根目录执行：

```bash
conda env create -n flora-ont330 -f environment.yml
conda activate flora-ont330

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONUTF8=1
```

已有环境时，用 `conda env update -n flora-ont330 -f environment.yml` 代替创建命令，
然后重新激活。必须使用 Python 3.11.x。
`environment.yml` 包含构建工具和运行依赖；`environment.runtime.yml` 用于
已打包的可执行文件；`requirements.txt` 只列出 Python 包。本分支还需要 StringTie
2.2.2、`gffread` 和 `gffcompare`，两个环境文件中已包含这些依赖。
建议将 `flora-ont330` 与 main 或其他 ONT 版本使用的环境分开。

编译前为 `rust-htslib`/`bindgen` 设置 libclang 路径。
这仅是构建依赖，运行已打包的 Flora 不需要单独设置。

```bash
LIBCLANG_SO="$(find -L "$CONDA_PREFIX/lib" -maxdepth 1 -type f \
  \( -name 'libclang.so' -o -name 'libclang.so.*' \) -print -quit)"
if [[ -n "$LIBCLANG_SO" ]]; then
  export LIBCLANG_PATH="$(dirname "$LIBCLANG_SO")"
  echo "LIBCLANG_PATH=$LIBCLANG_PATH"
else
  unset LIBCLANG_PATH
  echo "未找到 libclang，请更新 Conda 环境后再编译。" >&2
fi
```

找不到时，执行 `conda install -c conda-forge libclang`，然后重新查找；成功后再编译。
上述示例不会退出当前 SSH 会话。新开终端后需重新设置，也可在查找成功后执行
`conda env config vars set LIBCLANG_PATH="$LIBCLANG_PATH"`，再重新激活 `flora-ont330`。

## 开发构建与测试

每个准备发版的源码版本执行一次以下检查即可：

```bash
cargo test --locked
bash -n run_all.sh
bash -n run_all_mixed_species.sh
bash -n packaging/build_binary_release.sh
bash -n packaging/refresh_binary_release_metadata.sh
python -m unittest discover -s tests -v
```

每条命令均成功后再继续。`bash -n` 需要逐个检查脚本，不能在后面列出多个路径
来检查所有文件。这些测试不能代替代表性数据的全流程验证。

只生成开发用可执行文件、不打包时：

```bash
cargo build --release --locked
./target/release/flora --version
./target/release/flora --help
./target/release/flora mixed --help
./target/release/flora glycine --help
./target/release/flora analyze --help
```

`flora` 运行单物种全流程，`flora mixed` 运行混合物种全流程；
`flora glycine` 只做全长识别，`flora analyze` 只做 barcode/UMI 提取与细胞分配。
全流程支持 `--fastq` 后输入多个文件，或用 `--fastq-dir` 非递归读取目录。
Glycine 默认 10 个并发任务，共享 64 个总线程。

## 构建 Linux x86_64 发行包

完成环境准备和上述测试后执行：

```bash
PYTHON_BIN="$CONDA_PREFIX/bin/python" \
  bash packaging/build_binary_release.sh
```

打包脚本会自行执行锁定依赖的 release 编译，无需先单独执行 `cargo build`。
它会检查 Linux x86_64 可执行文件和 CLI，复制运行环境与公开文档，
并自动生成、验证 SHA256。Python 资源嵌入在 `flora` 内，不以独立源码或字节码文件分发。

输出如下，版本号自动读取 `Cargo.toml`：

```text
dist/Flora-<version>-linux-x86_64/
dist/Flora-<version>-linux-x86_64.tar.gz
dist/Flora-<version>-linux-x86_64.tar.gz.sha256
```

发行包内的 `BUILD_INFO.txt` 记录所需 glibc 版本。请在兼容的 Linux 主机上构建，
此脚本不能在 macOS 上生成 Linux 发行包。可选前缀 `FLORA_MAX_GLIBC=2.34`
会拒绝依赖更高版本符号的构建，但不会降低二进制的 glibc 要求。
较旧系统需使用兼容的构建环境或 [Singularity](docs/singularity-build_zh-CN.md)。

发布前，将压缩包解压到新目录，激活运行环境，独立检查：

```bash
./flora --version
./flora --help
./flora mixed --help
./flora glycine --help
./flora analyze --help
cat BUILD_INFO.txt
```

发版前还应运行代表性数据；帮助命令只能验证 CLI，不能证明分析结果正确。
压缩包应当只有 `flora` 一个可执行文件，不包含私有源码。

## ONT 对照后端

本分支默认使用 `--feature-backend ont`，参考
[wf-single-cell v3.3.0 历史源码快照](https://github.com/epi2me-labs/wf-single-cell/tree/d46f92121c79c28a02105cf2444029e91dc4b609)
独立实现。ONT 已撤回该版本；本分支固定源码提交和 StringTie 2.2.2 用于对照，
不是直接运行 ONT Nextflow，也不承诺数值结果完全一致。

后端按染色体构建参考引导的样本转录组，重新比对 reads，再按比对分数和覆盖度
选择 gene/transcript。Gene 默认要求基因组 MAPQ >=30；transcript 使用独立的
转录组比对证据。新发现的转录本会导出，但只统计能关联到参考注释的转录本。
Barcode/UMI、细胞判定、QC 和报告沿用 Flora 定义，不套用旧 ONT 的 10x 固定 UMI 长度过滤。

常用 Flora 命令无需增加参数；以下选项按需使用：

| 选项 | 作用 |
| --- | --- |
| `--feature-backend legacy` | 使用原有 GTF-overlap 后端对照，另设输出目录。 |
| `--ont-temp-dir /local_scratch/flora_tmp` | 将临时文件放到可写本地磁盘；默认使用系统临时目录/`TMPDIR`。 |
| `--ont-jobs 2` | 并行处理染色体，默认 1。 |
| `--stringtie-opts '-c 2'` | 设置 StringTie 重构选项。 |
| `--cell-gene-max-reads 0` | ONT 默认不按每个 cell+gene 的 read 数量截断。 |

实际染色体并发数不超过 `min(ont-jobs, max(1, cluster-threads / 8))`，
各任务共享线程预算，工具主线程和辅助 I/O 线程可能额外存在。
提高并发会增加内存和临时盘需求；大样本不要使用空间不足或内存文件系统上的 `/tmp`。
`--transcript-assign-mapq` 和 `--gene-tagging-backend` 仅影响 legacy；
`--skip-isoform` 跳过转录本矩阵，但 ONT 的 gene 判定仍需转录组重构。

除常规输出外，还包括：

- `transcriptome/<sample>.transcriptome.fa.gz` 和 `.transcriptome.gtf.gz`：重构序列和注释。
- `transcriptome/tool_versions.json` 和 `assignment_summary.json`：实际工具版本和分配状态计数。
- `matrix/` 下 gene/transcript 两个层级的 raw/processed MEX 矩阵，
  以及记录处理详情的 `<sample>.matrix_processing.json`。

Raw MEX 为整数 UMI 计数。Processed MEX 要求每细胞至少 200 个特征、每特征至少
3 个细胞、线粒体比例不超过 20%，再按每细胞 10,000 计数归一化并执行
`log10(1+x)`；同时识别 `MT-` 和 `mt-`。
报告仍使用 Flora 原始细胞判定，不是 processed MEX 额外过滤后的细胞集合。
转录组临时比对 BAM 仅在 `--save-intermediate` 或 `--full-output` 时保留。

详细设计和实测：[ONT 后端](docs/ont-v3.3.0-alignment.md)、
[性能优化记录](docs/ont-v330-performance.md)。

## 流程输出策略

默认轻量输出保留最终 tagged BAM/索引、表达矩阵、cell/barcode 对应关系、
RNA 聚类、QC、`qc/metrics_summary.xlsx`、日志和 HTML 报告。
成功结束后删除可重建的 read-level 表、加标签前的 aligned BAM、assignment TSV；
如果 `alignment/` 变空，也会移除该目录。

| 选项 | 作用 |
| --- | --- |
| `--save-intermediate` | 保留 read-level assignment 等调试中间文件。 |
| `--full-output` | 保留上游 FASTQ 输出和流程中间文件。 |
| `--remove-final-bam` | 成功分析后删除 `alignment/` 与 `matrix/` 下全部 BAM/索引，包括最终 tagged BAM。 |

这些参数控制文件保留，不改变分析阈值。BAM 删除选项与保留选项同时使用时，
仍会删除 BAM。全流程不会开启 Glycine 其他分类 FASTQ 的输出；需要这些文件时，
单独运行 `flora glycine` 并添加 `--keep-all-outputs`。

## 发布新版本

1. 修改 `Cargo.toml`，然后执行不带 `--locked` 的 `cargo check`，
   更新 `Cargo.lock` 中 Flora 自身的版本号。
2. 同步更新公开文档中的版本、release notes 和 `packaging/singularity/Flora.def`。
3. 按上述步骤测试、打包，并验证发行包。
4. 将 `docs/repository-templates/public/` 中审核后的文件同步到公开仓库。
   `.tar.gz` 和 `.sha256` 作为 Release 附件上传；SIF 和校验文件为可选附件。

Tag 与标题使用同一版本，例如 `v0.1.2` 和 `Flora v0.1.2`。
不要将压缩包、SIF、Docker tar、Miniforge 安装器或 `dist/` 提交到 Git。
不要向公开仓库复制私有的 `src/`、`scripts/`、`vendor/`、Cargo 文件、测试和打包代码。
编译后的可执行文件并非加密，仍可能被逆向。

## 源码目录

| 路径 | 用途 |
| --- | --- |
| `src/`、`vendor/`、`Cargo.toml`、`Cargo.lock` | Rust 源码与锁定依赖。 |
| `scripts/` | Python 分析、QC 与报告资源。 |
| `run_all.sh`、`run_all_mixed_species.sh` | 单物种与混合物种流程脚本。 |
| `tests/` | 回归与发行检查。 |
| `packaging/` | 二进制与容器打包工具。 |
| `docs/repository-templates/` | 公开/私有仓库文档模板。 |

容器详细步骤：[构建 Singularity 镜像](docs/singularity-build_zh-CN.md)。
