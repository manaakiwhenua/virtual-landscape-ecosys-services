library(landscapemetrics)
library(raster)
library(dplyr)

landscape <- raster(snakemake@params[['landscape']])
check_landscape(landscape)

landscape_metrics = dplyr::bind_rows(
  lsm_l_enn_mn(landscape),
  lsm_l_shdi(landscape),
  lsm_l_contag(landscape)
)

write.csv(landscape_metrics, file=snakemake@output[[1]], fileEncoding='UTF-8', na='', row.names=FALSE)
