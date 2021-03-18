FROM rocker/r-base:latest
RUN apt-get update -qq && apt-get -y install libgdal-dev libproj-dev

COPY ./install_packages.R ./install_packages.R
RUN Rscript ./install_packages.R
