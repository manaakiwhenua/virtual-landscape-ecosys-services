FROM python:3.8-buster
LABEL maintainer="lawr@landcareresearch.co.nz"
# InVEST code version
ARG VERSION=3.8.9
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get dist-upgrade --assume-yes
RUN apt-get install --assume-yes --no-install-recommends \
  libgdal20 \
  libgdal-dev \
  libspatialindex-c5 \
  python3-rtree

# numpy must be installed before GDAL to avoid problems with GDAL ReadAsArray
# https://gis.stackexchange.com/questions/83138/cannot-import-gdal-array
RUN pip3 install "numpy>=1.11.0,!=1.16.0"
ADD requirements.txt /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements.txt
RUN pip3 install natcap.invest==$VERSION
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal
RUN apt-get purge --assume-yes \
  libgdal-dev
RUN apt-get autoremove --assume-yes
RUN apt-get --assume-yes clean

ENTRYPOINT [ "bin/bash" ]
