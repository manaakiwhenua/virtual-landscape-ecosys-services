rule landscape_p_factor:
    input:
        "results/precipitation/{landscape}/{topography}.precipitation.tif"
    output:
        temp("results/erosion/{landscape}/{topography}.Landscape_P_Factor.tif")
    params:
        dem_offset="800",
        datatype="Float64"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/erosion/{landscape}/{topography}.Landscape_P_Factor.log"
    shell:
        'gdal_calc.py -A {input} --outfile={output} --calc="(A**2)*0.0012" --type={params.datatype} --overwrite --co TILED=YES > {log} 2>&1'

rule grassdata_dem_erosion:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
    output:
        directory("results/grass/data/region/erosion/{landscape}/{topography}")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        dem="/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
    log:
        "logs/erosion/{landscape}/{topography}.grass.log"
    shell:
        "grass -c {params.dem} -e {output} > {log} 2>&1"

rule slope:
    input:
        "results/grass/data/region/erosion/{landscape}/{topography}",
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
    output:
        temp("results/erosion/{landscape}/{topography}.Landscape_Z_Factor.tif")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        dem="/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif",
        dem_grass="{landscape}_{topography}_dem",
        slope="{landscape}_{topography}_slope",
        z_factor="{landscape}_{topography}_z_factor",
        grass="grass -f results/grass/data/region/erosion/{landscape}/{topography}/PERMANENT --text --exec"
    log:
        "logs/erosion/{landscape}/{topography}.Landscape_Z_Factor.log"
    shell:
        '({params.grass} r.in.gdal -r input="{params.dem}" output={params.dem_grass} --overwrite &&'
        '{params.grass} r.slope.aspect -e elevation={params.dem_grass} slope={params.slope} format="degrees" precision="FCELL" --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.z_factor} = 0.065 + ( 4.56 * ( {params.slope} / 100 ) ) + ( 65.41 * ({params.slope} / 100 ) ^ 2 )" --overwrite &&'
        '{params.grass} r.out.gdal -c -m input={params.z_factor} output={output} type=Float64 format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite '
        ') > {log} 2>&1'

rule slope_length_factor:
    input:
        "results/grass/data/region/erosion/{landscape}/{topography}",
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif"
    output:
        temp("results/erosion/{landscape}/{topography}.LandscapeSlopeLengthFactor.tif")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        dem="/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif",
        dem_grass="{landscape}_{topography}_dem",
        slope_length="{landscape}_{topography}_slope_length",
        slope_length_bounded="{landscape}_{topography}_slope_length_bounded",
        slope_length_factor="{landscape}_{topography}_slope_length_factor",
        flow_acc="{landscape}_{topography}_fa",
        pi="3.141592653589793",
        grass="grass -f results/grass/data/region/erosion/{landscape}/{topography}/PERMANENT --text --exec"
    log:
        "logs/erosion/{landscape}/{topography}.LandscapeSlopeLengthFactor.log"
    shell:
        '({params.grass} r.in.gdal -r input={params.dem} output={params.dem_grass} --overwrite &&'
        '{params.grass} r.flow elevation={params.dem_grass} flowaccumulation={params.flow_acc} --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.slope_length} = sqrt({params.flow_acc} * 625 / {params.pi})" --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.slope_length_bounded} = if({params.slope_length}>=350,350,{params.slope_length})" --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.slope_length_factor} = sqrt({params.slope_length} / 22)" --overwrite &&'
        '{params.grass} r.out.gdal -c -m input={params.slope_length_factor} output={output} type=Float64 format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite '
        ') > {log} 2>&1'

rule landscape_u_values:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    output:
        temp("results/erosion/{landscape}/{topography}.{landclass}.{proportion}.Landscape_U_Values.tif")
    params:
        datatype="Float64",
        input="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/erosion/{landscape}/{topography}.{landclass}.{proportion}.Landscape_U_Values.log"
    shell:
        'gdal_calc.py -A {params.input} --outfile={output} --calc="((A==0)*500 + (A==1)*100 + (A==2)*10 + (A==3)*5 + (A==4)*7 + (A==5)*5 + (A==6)*0)/1000.0" --type={params.datatype} --NoDataValue=0 --overwrite --co TILED=YES  > {log} 2>&1'

rule erosion:
    input:
        p="results/erosion/{landscape}/{topography}.Landscape_P_Factor.tif",
        z="results/erosion/{landscape}/{topography}.Landscape_Z_Factor.tif",
        slp="results/erosion/{landscape}/{topography}.LandscapeSlopeLengthFactor.tif",
        u="results/erosion/{landscape}/{topography}.{landclass}.{proportion}.Landscape_U_Values.tif"
    output:
        "results/erosion/{landscape}/{topography}.{landclass}.{proportion}.LandscapeErosion_t_km_yr.tif"
    params:
        datatype="Float64",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    log:
        "logs/erosion/{landscape}/{topography}.{landclass}.{proportion}.LandscapeErosion_t_km_yr.log"
    shell:
        'gdal_calc.py -A {input.p} -B {input.z} -C {input.slp} -D {input.u} --outfile={output} --calc="A * 0.2 * C * B * D" --overwrite --co TILED=YES > {log} 2>&1'

rule erosion_stats:
    input:
        f"{config['basepath']}"+"/virtual-landscapes/static/VirtualDomain/VirtualDomain.geojson",
        "results/erosion/{landscape}/{topography}.{landclass}.{proportion}.LandscapeErosion_t_km_yr.tif"
    output:
        "results/erosion/{landscape}/{topography}.{landclass}.{proportion}.LandscapeErosion_t_km_yr.json"
    conda:
        "../envs/rasterio.yml"
    params:
        stats='mean sum',
        all_touched=True
    log:
        "logs/erosion/{landscape}/{topography}.{landclass}.{proportion}.LandscapeErosion_t_km_yr.log"
    script:
        "../scripts/rasterstats.py > {log} 2>&1"
