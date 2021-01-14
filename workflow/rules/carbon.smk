rule carbon_stock:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        "results/carbon_stocks/{landscape}/{topography}.{landclass}.{proportion}.LandscapeCarbonStock.tif"
    params:
        datatype="UInt16",
        landscape="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/results/carbon_stocks/{landscape}/{topography}.{landclass}.{proportion}.LandscapeCarbonStock.log"
    shell:
        'gdal_calc.py -A {params.landscape} --outfile={output} --calc="(A==0)*5 + (A==1)*10 + (A==2)*18 + (A==3)*30 + (A==4)*140 + (A==5)*200" --type={params.datatype} --overwrite --co TILED=YES  > {log} 2>&1'
