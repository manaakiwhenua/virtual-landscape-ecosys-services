rule pollinator_requirement:
    input:
        f"{config['basepath']}"+"/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        temp("results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorRequirement.tif")
    params:
        datatype="UInt16",
        landscape="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorRequirement.log"
    shell:
        'gdal_calc.py -A {params.landscape} --outfile={output} --calc="numpy.isin(A,[0,3])*11 + numpy.isin(A,[1,2,4,5])*0 + (A==6)*999" --NoDataValue=999 --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'

rule pollinator_quality:
    input:
        f"{config['basepath']}"+"/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        temp("results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality.tif")
    params:
        datatype="UInt16",
        landscape="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality.log"
    shell:
        'gdal_calc.py -A {params.landscape} --outfile={output} --calc="numpy.isin(A,[0,1,6])*0 + numpy.isin(A,[2,4,5])*5 + numpy.isin(A,[3])*10" --NoDataValue=0 --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'

rule grassdata_landclass:
    input:
        f"{config['basepath']}"+"/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        directory("results/grass/data/region/pollination/{landscape}/{topography}/{landclass}/{proportion}")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        landscape="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
    log:
        "logs/grass/data/region/pollination/{landscape}/{topography}/{landclass}/{proportion}.grass.log"
    shell:
        "grass -c {params.landscape} -e {output} > {log} 2>&1"

rule pollinator_quality_grow:
    input:
        "results/grass/data/region/pollination/{landscape}/{topography}/{landclass}/{proportion}",
        q="results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality.tif"
    output:
        d=temp("results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_distance.tif"),
        v=temp("results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_value.tif")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        quality_layer="{landscape}_{topography}_{landclass}_{proportion}_LandscapePollinatorQuality_grow",
        metric="euclidean",
        grass="grass -c -f /results/grass/data/region/pollination/{landscape}/{topography}/{landclass}/{proportion}/PERMANENT --text --exec"
    log:
        "logs/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_value.log"
    shell: # Import to GRASS, "grow", export the distance raster and the value raster, and use them in gdal_calc classification
        '({params.grass} r.in.gdal -r input={input.q} output={params.quality_layer} --overwrite && '
        '{params.grass} r.grow.distance -m input={params.quality_layer} distance={params.quality_layer}_distance value={params.quality_layer}_value metric={params.metric} --overwrite && '
        '{params.grass} r.out.gdal -c -m input={params.quality_layer}_distance output={output.d} type=Float64 format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite && '
        '{params.grass} r.out.gdal -c -m input={params.quality_layer}_value output={output.v} type=Float64 -f format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite '
        ') > {log} 2>&1'

rule pollinator_quality_expand:
    input:
        q="results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality.tif",
        d="results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_distance.tif",
        v="results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_value.tif"
    output:
        temp("results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_expand.tif")
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    params:
        datatype="UInt16",
        distance="500" # 500 metres
    log:
        "logs/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_expand.log"
    shell: # TODO --hideNoData ? is output correct...?
        'gdal_calc.py -A {input.d} -B {input.v} -C {input.q} --outfile={output} --calc="(C==5)*C + (C==10)*C + logical_and(A>0,A<={params.distance})*B" --type={params.datatype} --overwrite --co TILED=YES --hideNoData > {log} 2>&1'

rule pollinator_overlap:
    input:
        q="results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorQuality_expand.tif",
        r="results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorRequirement.tif"
    output:
        "results/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorOverlap.tif"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    params:
        datatype="UInt16"
    log:
        "logs/pollination/{landscape}/{topography}.{landclass}.{proportion}.LandscapePollinatorOverlap.log"
    shell:
        'gdal_calc.py -A {input.q} -B {input.r} --outfile={output} --calc="A+B" --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'
