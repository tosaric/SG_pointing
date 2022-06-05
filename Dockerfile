# The build-stage image:
# FROM python:3.9-slim-buster
# FROM debian:buster-slim
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN groupadd --gid 11010 ccddev && useradd --uid 11010 --gid ccddev --shell /bin/bash --create-home ccddev
#RUN useradd --shell /bin/bash --create-home ccddev
RUN echo "root:root" | chpasswd && echo "ccddev:ccddev" | chpasswd

ENV SG_dir /home/ccddev/Starguider
RUN mkdir -p $SG_dir && chown -R ccddev:ccddev $SG_dir

WORKDIR /tmp

#RUN ./install_astrometry.sh

# basic linux install/build libs, python, and other support
# RUN cd $WORKDIR && \
# apt -y update && \
# apt install -y apt-utils make build-essential python3 python3-pip python3-venv python3-wheel python3-setuptools netpbm libnetpbm10-dev zlib1g-dev libcairo2-dev libjpeg-dev libcfitsio-dev libbz2-dev wget swig wcslib-dev \
# && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN cd $WORKDIR && \
apt-get update -y && apt install -y apt-utils && \
    apt install -y --no-install-recommends \
    build-essential \
    make \
    vim \
    gcc \
    git \
    file \
    pkg-config \
    wget \
    curl \
    swig \
    netpbm \
    wcslib-dev \
    wcslib-tools \
    zlib1g-dev \
    libbz2-dev \
    libcairo2-dev \
    libcfitsio-dev \
    libcfitsio-bin \
    libgsl-dev \
    libjpeg-dev \
    libnetpbm10-dev \
    libpng-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-pil \
    python3-tk \
    python3-setuptools \
    python3-wheel \
    # python3-numpy \
    # python3-scipy \
    python3-venv \
    # python3-matplotlib \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install dependencies:
SHELL ["/bin/bash", "-c"]

USER ccddev

WORKDIR $SG_dir

RUN mkdir -p $SG_dir/venv

ENV VIRTUAL_ENV=$SG_dir/venv

RUN python3 -m venv $VIRTUAL_ENV

ENV PATH="$VIRTUAL_ENV/bin:$PATH"

COPY requirements.txt .
RUN pip3 install -r requirements.txt

USER root

WORKDIR /src
COPY app .

RUN tar xvzf astrometry.net-latest.tar.gz \
&& cd */. \
&& make \
&& make py \
&& make extra \
&& make install

COPY index/index*.fits /usr/local/astrometry/data/

USER ccddev

#RUN mkdir -p /project/SG_astrometry/Data && mkdir -p /project/SG_astrometry/Code
RUN mkdir -p $SG_dir/Data && mkdir -p $SG_dir/Code
COPY Code $SG_dir/Code

ENV VIRTUAL_ENV=$SG_dir/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# ENV PYTHONPATH=/usr/local/astrometry/lib/python
# RUN pip3 install astropy==4.2
#SHELL ["/bin/bash", "-c"]
#CMD ["/bin/bash"]
#CMD ["python3", "/project/SG_astrometry/Code/SG_solve_1.1.py"]
####CMD ["/bin/sh", "-c", "python3 -u /project/SG_astrometry/Code/SG_solve_1.1.py 2>&1 | tee -a /project/SG_astrometry/Data/SG_output.log"]
CMD ["/bin/sh", "-c", "python3 -u $SG_dir/Code/SG_solve_1.2.py 2>&1 | tee -a $SG_dir/Data/SG_output.log"]
#CMD ["./project/SG_astrometry/Code/run_SG_solve.sh"]
#CMD ["/bin/sh", "-c", "python /project/SG_astrometry/Code/SG_solve_1.1.py > /local/home/ccddev/SG_plots/log_output.log 2>&1"]
#CMD ["python", "/project/SG_astrometry/Code/SG_solve_1.1.py 2>&1 | tee /local/home/ccddev/SG_plots/log_output.txt"]
#CMD ["/bin/bash", "-c", "python /project/SG_astrometry/Code/SG_solve_1.1.py 2>&1 | tee -a /local/home/ccddev/SG_plots/log_output.txt"]
#CMD ["/bin/bash", "-c", "python /project/SG_astrometry/Code/SG_solve_1.1.py >> /local/home/ccddev/SG_plots/log_output.txt 2>&1"]
#RUN -v /home/toni/projects/astrometry/Test_2.4/SG_images/:/project/SG_astrometry/Data/
#COPY Data /project/SG_astrometry/Data/

# ENV VIRTUAL_ENV=/project/venv
# ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# SHELL ["/bin/bash", "-c"]


# CMD ["python3", "/project/SG_astrometry/Data/Test_1_1.py"]
