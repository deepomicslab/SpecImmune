# 使用 Ubuntu 作为基础镜像
FROM ubuntu:latest

# 安装必要的软件包和依赖
RUN apt-get update && apt-get install -y \
    wget \
    software-properties-common \
    curl \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 设置工作目录
WORKDIR /SpecImmune

# 复制代码和脚本到容器中
COPY . /SpecImmune

# 下载并安装 Miniconda
RUN wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

# 将 Conda 添加到 PATH
ENV PATH="/opt/conda/bin:$PATH"

# 更新 Conda
RUN conda update -n base -c defaults conda -y

# 安装 mamba
#RUN conda install -n base -c conda-forge mamba -y

# 创建 Conda 环境
RUN conda env create -n SpecImmune -f /SpecImmune/environment.yml

# 初始化 Conda 并激活环境
RUN echo "source /opt/conda/bin/activate SpecImmune" >> ~/.bashrc

# 删除不需要的目录, 也可以由 .dockerignore 文件来直接忽略
#RUN rm -rf /SpecImmune/test \
#   /SpecImmune/pharmvar-6.2 \
#   /SpecImmune/simulation \
#   /SpecImmune/evaluation \
#   /SpecImmune/SpecImmune

# 设置 ENTRYPOINT
ENTRYPOINT ["/bin/bash", "-c", "echo \"Passed arguments: $@\"; source /opt/conda/bin/activate SpecImmune && exec python /SpecImmune/scripts/main.py \"$@\"", "--"]
