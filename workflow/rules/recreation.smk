rule grassdata_dem_recreation:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
    output:
        directory("results/grass/data/region/recreation/{landscape}/{topography}")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        dem="/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
    log:
        "logs/grass/data/region/recreation/{landscape}/{topography}.grass.log"
    shell:
        "grass -c {params.dem} -e {output} > {log} 2>&1"

rule grassdata_dem_recreation_landclass:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        directory("results/grass/data/region/recreation/landclass/{landscape}/{topography}/{landclass}/{proportion}")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        landclass="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
    log:
        "logs/grass/data/region/recreation/landclass/{landscape}/{topography}/{landclass}/{proportion}.grass.log"
    shell:
        "grass -c {params.landclass} -e {output} > {log} 2>&1"

RADIUS = 200 # metres
CELL_SIZE = 25 # assumes 25 metre cells
NEIGHBOURHOOD_SIZE = int((RADIUS*2)/CELL_SIZE)
rule landscape_scenic:
    input:
        "results/grass/data/region/recreation/{landscape}/{topography}",
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
    output:
        temp("results/recreation/{landscape}/{topography}.LandscapeScenicValueNormalised.tif")
    params:
        dem="/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif",
        dem_grass="{landscape}_{topography}_dem",
        focal="{landscape}_{topography}_focal",
        tpi="{landscape}_{topography}_tpi",
        tpi_norm="{landscape}_{topography}_tpi_norm",
        grass="grass -f results/grass/data/region/recreation/{landscape}/{topography}/PERMANENT --text --exec",
        datatype="Float64",
        method="average",
        threshold=1.5,
        size=(str(NEIGHBOURHOOD_SIZE if (NEIGHBOURHOOD_SIZE % 2 != 0) else (NEIGHBOURHOOD_SIZE + 1)))
    container:
        "docker://neteler/grassgis7:latest"
    log:
        "logs/recreation/{landscape}/{topography}.LandscapeScenicValueNormalised.log"
    shell:
        '({params.grass} r.in.gdal -r input={params.dem} output={params.dem_grass} --overwrite && '
        '{params.grass} r.neighbors -c input={params.dem_grass} output={params.focal} method={params.method} size={params.size} --overwrite && '
        '{params.grass} r.mapcalc expression="{params.tpi} = {params.dem_grass} - {params.focal}" --overwrite && '
        '{params.grass} r.mapcalc expression="{params.tpi_norm} = if({params.tpi}<{params.threshold},0,1)" --overwrite &&'
        '{params.grass} r.out.gdal -c -m input={params.tpi_norm} output={output} type={params.datatype} format=GTiff createopt="COMPRESS=DEFLATE" --overwrite '
        ') > {log} 2>&1'

rule landscape_natural_value:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        temp("results/recreation/{landscape}/{topography}.{landclass}.{proportion}.LandscapeNaturalValueNormalised.tif")
    params:
        datatype="Float64",
        landscape="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/recreation/{landscape}/{topography}.{landclass}.{proportion}.LandscapeNaturalValueNormalised.log"
    shell:
        'gdal_calc.py -A {params.landscape} --outfile={output} --calc="(A==0)*0 + (A==1)*0.2 + (A==2)*0.7 + (A==3)*0.6 + (A==4)*0.5 + (A==5)*1.0 + (A==6)*999" --NoDataValue=999 --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'

rule landscape_areal_proportion:
    input:
        "results/grass/data/region/recreation/landclass/{landscape}/{topography}/{landclass}/{proportion}",
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        temp("results/recreation/{landscape}/{topography}.{landclass}.{proportion}.Landscape_Landcover_Area_Proportion.tif")
    params:
        grass="grass -c -f results/grass/data/region/recreation/landclass/{landscape}/{topography}/{landclass}/{proportion}/PERMANENT --text --exec",
        landclass="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
        landclass_grass="{landscape}_{topography}_{landclass}_{proportion}_landclass",
        total_area=config['RESOLUTION']^2*config['HEIGHT']*config['WIDTH']
    container:
        "docker://neteler/grassgis7:latest"
    log:
        "logs/recreation/{landscape}/{topography}.{landclass}.{proportion}.Landscape_Landcover_Area_Proportion.log"
    shell:
        '({params.grass} r.in.gdal -r input={params.landclass} output={params.landclass_grass} --overwrite && '
        '{params.grass} r.clump -d input={params.landclass_grass} output={params.landclass_grass}_clump --overwrite && '
        '{params.grass} r.to.vect input={params.landclass_grass}_clump output={params.landclass_grass}_patch type=area --overwrite && '
        '{params.grass} v.db.addcolumn map={params.landclass_grass}_patch columns="Area DOUBLE PRECISION" && '
        '{params.grass} v.to.db map={params.landclass_grass}_patch option=area units=meters columns=Area --overwrite && '
        '{params.grass} v.to.rast input={params.landclass_grass}_patch output={params.landclass_grass}_patch_r use=attr attribute_column=Area --overwrite && '
        '{params.grass} r.mapcalc expression="{params.landclass_grass}_patch_r_proportion=({params.landclass_grass}_patch_r/{params.total_area})" --overwrite && '
        '{params.grass} r.out.gdal -c -m input={params.landclass_grass}_patch_r_proportion output={output} type=Float64 format=GTiff createopt="COMPRESS=DEFLATE" --overwrite '
        ') > {log} 2>&1'

rule landscape_natural_value_patch_weighted:
    input:
        lc="results/recreation/{landscape}/{topography}.{landclass}.{proportion}.Landscape_Landcover_Area_Proportion.tif",
        lv="results/recreation/{landscape}/{topography}.{landclass}.{proportion}.LandscapeNaturalValueNormalised.tif"
    output:
        temp("results/recreation/{landscape}/{topography}.{landclass}.{proportion}.Landscape_Natural_Value_Patch.tif")
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    params:
        datatype="Float64"
    log:
        "logs/recreation/{landscape}/{topography}.{landclass}.{proportion}.Landscape_Natural_Value_Patch.log"
    shell:
        'gdal_calc.py -A {input.lc} -B {input.lv} --outfile={output} --calc="A*B" --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'

rule rivers_buffer:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-rivers_{topography}.tif"
    output:
        temp("results/recreation/{landscape}/{topography}.{landclass}.LandscapeRiverExpand.tif")
    params:
        datatype="Byte",
        rivers="/virtual-landscapes/{landscape}/{landscape}-rivers_{topography}.tif",
        distance=20,
        units="PIXEL"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/recreation/{landscape}/{topography}.{landclass}.LandscapeRiverExpand.log"
    shell:
        '(gdal_proximity.py {params.rivers} {output} -ot {params.datatype} -distunits {params.units} -maxdist {params.distance} -fixed-buf-val 1 -nodata 0 -of GTIff && '
        'gdal_calc.py -A {output} -B {params.rivers} --outfile={output} --calc="A+B" --type={params.datatype} --overwrite --co TILED=YES'
        ') > {log} 2>&1'

rule recreation_value_additive:
    input:
        r="results/recreation/{landscape}/{topography}.{landclass}.LandscapeRiverExpand.tif",
        v="results/recreation/{landscape}/{topography}.LandscapeScenicValueNormalised.tif",
        vp="results/recreation/{landscape}/{topography}.{landclass}.{proportion}.Landscape_Natural_Value_Patch.tif"
    output:
        "results/recreation/{landscape}/{topography}.{landclass}.{proportion}.RecreationValueAdditive.tif"
    params:
        datatype="Float64"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/recreation/{landscape}/{topography}.{landclass}.{proportion}.RecreationValueAdditive.tif.log"
    shell:
        'gdal_calc.py -A {input.r} -B {input.v} -C {input.vp} --outfile={output} --calc="A*0.2 + B*0.2 + C" --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'
