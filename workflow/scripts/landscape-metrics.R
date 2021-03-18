library(landscapemetrics)
library(raster)

landscape <- raster(snakemake@params[["landscape"]])
check_landscape(landscape)

landscape_metrics = dplyr::bind_rows(
  lsm_l_shdi(landscape),
  lsm_l_contag(landscape)
)

print(snakemake@output[[1]])
write.csv(landscape_metrics, file=snakemake@output[[1]], fileEncoding='UTF-8')
