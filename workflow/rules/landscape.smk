rule landscape_precipitation:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
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
