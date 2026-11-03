[English](singularity-build.md) | [中文](singularity-build_zh-CN.md)

# 构建 Singularity 镜像

本教程于 2026-10-27 从私有 README 移出，以精简主文档。
它将已有 Linux 发行包封装为容器，不改变 Flora 的分析逻辑。

## 准备条件

- 构建机为 Linux x86_64，已安装 Singularity，且可以使用 `--fakeroot`。
- Flora 发行压缩包，以及版本对应的 `packaging/singularity/Flora.def`。
- 能联网下载 Ubuntu 软件包和 Conda 依赖。
- 仅从 Docker 镜像准备基础镜像时需要 Docker；也可以直接复制已有的基础 SIF。

以下示例对应当前定义文件的 0.1.2 版本。其他版本需要同步修改定义文件中的
版本标签、压缩包名称和下面的命令。基础系统为 Ubuntu 22.04（glibc 2.35）；
请检查发行包 `BUILD_INFO.txt` 中的 `required_glibc` 不高于容器内版本。
构建 SIF 时不需要激活 Flora Conda 环境。

## 1. 准备构建目录

生成发行压缩包后，在私有仓库根目录执行：

```bash
mkdir -p singularity_build
cp dist/Flora-0.1.2-linux-x86_64.tar.gz singularity_build/
cp packaging/singularity/Flora.def singularity_build/
cd singularity_build
```

## 2. 准备基础镜像

已有 `flora-base-ubuntu22.04.sif` 时，将它复制到本目录，跳过此步骤。
否则，检查并导出 Docker 镜像：

```bash
docker run --rm quay.io/nf-core/ubuntu:22.04 \
  bash -c 'cat /etc/os-release; getconf GNU_LIBC_VERSION'

docker save quay.io/nf-core/ubuntu:22.04 \
  -o flora-base-ubuntu22.04.tar

singularity build --fakeroot \
  flora-base-ubuntu22.04.sif \
  docker-archive://flora-base-ubuntu22.04.tar
```

检查结果应为 Ubuntu 22.04、glibc 2.35。本地未缓存 Docker 镜像时，需要访问镜像仓库。
不要通过将 Docker socket 设置为全局可写来解决权限问题；请使用有权限的账号或
联系管理员。导入已有 Docker archive 不需要访问 Docker daemon。

## 3. 下载本地 Miniforge 安装器

`Flora.def` 从宿主机复制安装器，不在构建中下载它。
可以复用可信的已有安装器，或先下载一次：

```bash
curl -L --fail --retry 5 --retry-delay 5 \
  -o Miniforge3-Linux-x86_64.sh \
  https://mirror.nju.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh
```

此时构建目录中必须有四个文件：

```text
Flora.def
Flora-0.1.2-linux-x86_64.tar.gz
Miniforge3-Linux-x86_64.sh
flora-base-ubuntu22.04.sif
```

## 4. 构建与验证

```bash
env -u LD_LIBRARY_PATH \
singularity build --fakeroot \
  Flora-0.1.2-linux-x86_64.sif \
  Flora.def

singularity run --cleanenv Flora-0.1.2-linux-x86_64.sif --version
singularity run --cleanenv Flora-0.1.2-linux-x86_64.sif --help

sha256sum Flora-0.1.2-linux-x86_64.sif \
  > Flora-0.1.2-linux-x86_64.sif.sha256
sha256sum -c Flora-0.1.2-linux-x86_64.sif.sha256
```

每条命令成功后再继续。定义文件中的 `%test` 会检查 Python 包导入和必要的
可执行工具。最终镜像大小取决于实际解析的依赖；缓存目录默认配置在 `/tmp` 下。

若 `--fakeroot` 提示缺少 subuid 映射，需要管理员配置构建机，或换到满足条件的
主机构建，再复制最终 SIF。若运行时提示 user namespace 被禁用，则目标服务器
需要管理员配置可用的 Singularity 运行环境；修改 Flora 文件不能绕过该限制。

## 5. 运行与分发

用户只需要最终 SIF 和可用的 Singularity，不需要 Docker、基础镜像或单独的
Flora Conda 环境。请绑定流程所需的数据、参考和输出父目录，例如：

```bash
singularity run --cleanenv \
  --bind /data:/data \
  Flora-0.1.2-linux-x86_64.sif \
  --fastq-dir /data/sample_fastq_dir \
  --barcode-list-10bp /data/BC_1536.txt \
  --ref-dir /data/GRCh38_flora \
  --out-dir /data/sample_output \
  --sample-id sample
```

最终 SIF 和校验文件可作为 GitHub Release 的可选附件。
不要把镜像、Docker archive 或 Miniforge 安装器提交到 Git。

[返回私有 README](../README_zh-CN.md)
