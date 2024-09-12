# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

FROM ubuntu:20.04 AS base
ENV DEBIAN_FRONTEND=noninteractive

COPY DockerImport/OPERA2.9_CL_mcr.tar.gz OPERA.tar.gz
RUN tar -xzf OPERA.tar.gz && rm OPERA.tar.gz
RUN apt-get update && apt-get install -y sudo && apt-get install -y unzip 
RUN apt-get install -y libxt6
RUN sudo OPERA2_CL_mcr/OPERA2.9_mcr_Installer.install -mode silent -agreeToLicense yes
ENV LD_LIBRARY_PATH=/usr/local/MATLAB/MATLAB_Runtime/v912/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v912/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v912/sys/os/glnxa64:MR/v912/sys/opengl/lib/glnxa64
ENV XAPPLRESDIR=/usr/local/MATLAB/MATLAB_Runtime/v912/X11/app-defaults
RUN /usr/OPERA/application/OPERA
RUN echo "/usr/OPERA/application" > /root/.mcrCache9.12/OPERA_installdir.txt

FROM python:3.9

RUN apt-get update
RUN apt-get install -y vim 
RUN apt-get install -y nano

WORKDIR /app

COPY DockerImport/volumes/etcd /etcd_data
COPY DockerImport/volumes/milvus /milvus_data
COPY DockerImport/volumes/minio /minio_data

# Define the volume mount points
VOLUME ["/etcd_data", "/milvus_data", "/minio_data"]


RUN mkdir javadir
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "amd64" ]; then \
        curl -k https://download.oracle.com/java/20/archive/jdk-20.0.2_linux-x64_bin.tar.gz | tar -xvz -C javadir; \
    elif [ "$TARGETARCH" = "arm64" ]; then \
        curl -k https://download.oracle.com/java/20/archive/jdk-20.0.2_linux-aarch64_bin.tar.gz | tar -xvz -C javadir; \
    else \
        echo "Unsupported architecture: $TARGETARCH"; exit 1; \
    fi

ENV JAVA_HOME=/app/javadir/jdk-20.0.2
ENV PATH=${PATH}:/app/javadir/jdk-20.0.2/bin
ENV FLASK_APP=App/App.py
ENV LD_LIBRARY_PATH=/usr/local/MATLAB/MATLAB_Runtime/v912/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v912/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v912/sys/os/glnxa64:MR/v912/sys/opengl/lib/glnxa64
ENV XAPPLRESDIR=/usr/local/MATLAB/MATLAB_Runtime/v912/X11/app-defaults
EXPOSE 5005

COPY --from=base OPERA2_CL_mcr OPERA2_CL_mcr
COPY --from=base /root/.mcrCache9.12 /root/.mcrCache9.12
COPY --from=base /usr/OPERA /usr/OPERA
COPY --from=base /usr/local/MATLAB/MATLAB_Runtime/v912 /usr/local/MATLAB/MATLAB_Runtime/v912

COPY src/Docker/requirements.txt requirements.txt
RUN pip --trusted-host pypi.org --trusted-host files.pythonhosted.org install -r requirements.txt

COPY src/Comparison src/Comparison
COPY src/App src/App
COPY src/main.py src/main.py
COPY test.txt test.txt 
# temp TEST ABOVE ^^ NEED TO REMOVE

RUN rm -r src/Comparison/LocalIO/*
RUN rm -r src/App/static/LocalIO/*
RUN mkdir output
RUN mkdir logs

LABEL org.opencontainers.image.source https://github.com/sandialabs/chemical-recommender-system

CMD ["/bin/sh", "-c", "python /app/src/main.py -w && tail -f /dev/null"]