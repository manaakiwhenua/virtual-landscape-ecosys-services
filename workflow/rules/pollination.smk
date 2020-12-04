rule pollinator_requirement:
    input:
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif"
    output:
        temp("results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorRequirement.tif")
    params:
        datatype="UInt16",
        landscape="/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    shell:
        'gdal_calc.py -A {params.landscape} --outfile={output} --calc="numpy.isin(A,[0,3])*11 + numpy.isin(A,[1,2,4,5])*0 + (A==6)*999" --NoDataValue=999 --type={params.datatype} --overwrite --co TILED=YES'

rule pollinator_quality:
    input:
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif"
    output:
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality.tif"
    params:
        datatype="UInt16",
        landscape="/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    shell:
        'gdal_calc.py -A {params.landscape} --outfile={output} --calc="numpy.isin(A,[0,1,6])*0 + numpy.isin(A,[2,4,5])*5 + numpy.isin(A,[3])*10" --NoDataValue=0 --type={params.datatype} --overwrite --co TILED=YES'

rule grassdata_landclass:
    input:
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif"
    output:
        directory("results/grass/data/region/{landscape}/{topography}/{landclass}")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        landscape="/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif",
    shell:
        "grass -c {params.landscape} -e {output}"

rule pollinator_quality_grow:
    input:
        "results/grass/data/region/{landscape}/{topography}/{landclass}",
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality.tif"
    output:
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality_distance.tif",
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality_value.tif"
    container:
        "docker://neteler/grassgis7:latest"
    params:
        quality_layer="{landscape}_{topography}_{landclass}_LandscapePollinatorQuality_grow",
        metric="euclidean",
        grass="grass /results/grass/data/region/{landscape}/{topography}/{landclass}/PERMANENT --text --exec"
    shell: # Import to GRASS, "grow", export the distance raster and the value raster, and use them in gdal_calc classification
        '{params.grass} r.in.gdal -r input={input[1]} output={params.quality_layer} --overwrite && '
        '{params.grass} r.grow.distance -m input={params.quality_layer} distance={params.quality_layer}_distance value={params.quality_layer}_value metric={params.metric} --overwrite && '
        '{params.grass} r.out.gdal -c -m input={params.quality_layer}_distance output={output[0]} type=Float64 format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite && '
        '{params.grass} r.out.gdal -c -m input={params.quality_layer}_value output={output[1]} type=Float64 -f format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite'

rule pollinator_quality_expand:
    input:
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality.tif",
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality_distance.tif",
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality_value.tif"
    output:
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality_expand.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    params:
        datatype="UInt16",
        distance="500" # 500 metres
    shell: # TODO --hideNoData ? is output correct...?
        'gdal_calc.py -A {input[1]} -B {input[2]} -C {input[0]} --outfile={output} --calc="(C==5)*C + (C==10)*C + logical_and(A>0,A<={params.distance})*B" --type={params.datatype} --overwrite --co TILED=YES --hideNoData'

rule pollinator_overlap:
    input:
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorQuality_expand.tif",
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorRequirement.tif"
    output:
        "results/pollination/{landscape}/{topography}.{landclass}.LandscapePollinatorOverlap.tif"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    params:
        datatype="UInt16"
    shell:
        'gdal_calc.py -A {input[0]} -B {input[1]} --outfile={output} --calc="A+B" --type={params.datatype} --overwrite --co TILED=YES'
