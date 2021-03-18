rule landscape_precipitation:
    input:
        f"{config['basepath']}"+"/{landscape}/{landscape}-dem_{topography}.tif"
    output:
        "results/precipitation/{landscape}/{topography}.precipitation.tif"
    params:
        dem="/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif",
        dem_offset="800",
        datatype="Float64"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/precipitation/{landscape}/{topography}.precipitation.log"
    shell:
        'gdal_calc.py -A {params.dem} --outfile={output} --calc="A+{params.dem_offset}" --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'

rule landscape_metrics:
    message: "Calculating landscape metrics (Shannon's diversity index, and contagion index) for {input}"
    input:
        f"{config['basepath']}"+"/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        "results/landscape/{landscape}/{topography}.{landclass}.{proportion}.metrics.csv"
    params:
        landscape="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    container:
        "docker://richardlaw/r-base-landscapemetrics:v0.1"
    log:
        "logs/landscape/{landscape}/{topography}.{landclass}.{proportion}.SHDI.log"
    script:
        "../scripts/landscape-metrics.R"
