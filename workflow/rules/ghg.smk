rule greenhouse_gas_emission:
    input:
        f"{config['basepath']}"+"/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        "results/greenhouse_gas_emissions/{landscape}/{topography}.{landclass}.{proportion}.LandscapeGreenhouseGasEmissions.tif"
    params:
        datatype="UInt16",
        landscape="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/greenhouse_gas_emissions/{landscape}/{topography}.{landclass}.{proportion}.LandscapeGreenhouseGasEmissions.log"
    shell:
        'gdal_calc.py -A {params.landscape} --outfile={output} --calc="(A==0)*1400 + (A==1)*8800 + (A==2)*1000 + (A==3)*2 + (A==4)*5 + (A==5)*0" --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'
