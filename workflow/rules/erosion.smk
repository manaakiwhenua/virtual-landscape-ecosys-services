rule landscape_p_factor:
    input:
        "results/precipitation/{landscape}/{topography}.precipitation.tif"
    output:
        "results/erosion/{landscape}/{topography}.Landscape_P_Factor.tif"
    params:
        dem_offset="800",
        datatype="Float64"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    shell:
        'gdal_calc.py -A {params.dem} --outfile={output} --calc="({input}**2)*0.0012" --type={params.datatype} --overwrite --co TILED=YES'

rule grassdata_dem_erosion:
    input:
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/dem_{topography}.tif"
    output:
        directory("results/grass/data/region/erosion/{landscape}/{topography}")
    container:
        "docker://neteler/grassgis7:latest"
    params:
        dem="/virtual-landscapes/{landscape}/dem_{topography}.tif",
    shell:
        "grass -c {params.dem} -e {output}"

rule slope:
    input:
        "results/grass/data/region/erosion/{landscape}/{topography}",
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/dem_{topography}.tif"
    output:
        "results/erosion/{landscape}/{topography}.Landscape_Z_Factor.tif"
    container:
        "docker://neteler/grassgis7:latest"
    params:
        dem="/virtual-landscapes/{landscape}/dem_{topography}.tif",
        dem_grass="{landscape}_{topography}_dem",
        slope="{landscape}_{topography}_slope",
        z_factor="{landscape}_{topography}_z_factor",
        grass="grass results/grass/data/region/erosion/{landscape}/{topography}/PERMANENT --text --exec"
    shell:
        '{params.grass} r.in.gdal -r input="{params.dem}" output={params.dem_grass} --overwrite &&'
        '{params.grass} r.slope.aspect -e elevation={params.dem_grass} slope={params.slope} format="degrees" precision="FCELL" --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.z_factor} = 0.065 + ( 4.56 * ( {params.slope} / 100 ) ) + ( 65.41 * ({params.slope} / 100 ) ^ 2 )" --overwrite &&'
        '{params.grass} r.out.gdal -c -m input={params.z_factor} output={output} type=Float64 format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite'

rule slope_length_factor:
    input:
        "results/grass/data/region/erosion/{landscape}/{topography}",
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/dem_{topography}.tif"
    output:
        "results/erosion/{landscape}/{topography}.LandscapeSlopeLengthFactor.tif"
    container:
        "docker://neteler/grassgis7:latest"
    params:
        dem="/virtual-landscapes/{landscape}/dem_{topography}.tif",
        dem_grass="{landscape}_{topography}_dem",
        slope_length="{landscape}_{topography}_slope_length",
        slope_length_bounded="{landscape}_{topography}_slope_length_bounded",
        slope_length_factor="{landscape}_{topography}_slope_length_factor",
        flow_acc="{landscape}_{topography}_fa",
        pi="3.141592653589793",
        grass="grass results/grass/data/region/erosion/{landscape}/{topography}/PERMANENT --text --exec"
    shell:
        '{params.grass} r.in.gdal -r input={params.dem} output={params.dem_grass} --overwrite &&'
        '{params.grass} r.flow elevation={params.dem_grass} flowaccumulation={params.flow_acc} --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.slope_length} = sqrt({params.flow_acc} * 625 / {params.pi})" --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.slope_length_bounded} = if({params.slope_length}>=350,350,{params.slope_length})" --overwrite &&'
        '{params.grass} r.mapcalc expression="{params.slope_length_factor} = sqrt({params.slope_length} / 22)" --overwrite &&'
        '{params.grass} r.out.gdal -c -m input={params.slope_length_factor} output={output} type=Float64 format=GTiff createopt="TFW=YES,COMPRESS=DEFLATE" --overwrite'

rule landscape_u_values:
    input:
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif"
    output:
        "results/erosion/{landscape}/{topography}.{landclass}.Landscape_U_Values.tif"
    params:
        datatype="Float64",
        input="/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif"
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    shell:
        'gdal_calc.py -A {params.input} --outfile={output} --calc="((A==0)*500 + (A==1)*100 + (A==2)*10 + (A==3)*5 + (A==4)*7 + (A==5)*5 + (A==6)*0)/1000.0" --type={params.datatype} --NoDataValue=0 --overwrite --co TILED=YES'

rule erosion:
    input:
        "results/precipitation/{landscape}/{topography}.Landscape_P_Factor.tif",
        "results/erosion/{landscape}/{topography}.Landscape_Z_Factor.tif",
        "results/erosion/{landscape}/{topography}.LandscapeSlopeLengthFactor.tif",
        "results/erosion/{landscape}/{topography}.{landclass}.Landscape_U_Values.tif"
    output:
        "results/erosion/{landscape}/{topography}.{landclass}.LandscapeErosion_t_km_yr.tif"
    params:
        datatype="Float64",
    container:
        "docker://osgeo/gdal:ubuntu-small-latest"
    shell:
        'gdal_calc.py -A {input[0]} -B {input[1]} -C {input[2]} -D {input[3]} --outfile={output} --calc="A * 0.2 * C * B * D" --overwrite --co TILED=YES'

rule erosion_stats:
    input:
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/static/VirtualDomain/VirtualDomain.geojson",
        "results/erosion/{landscape}/{topography}.{landclass}.LandscapeErosion_t_km_yr.tif"
    output:
        "results/erosion/{landscape}/{topography}.{landclass}.LandscapeErosion_t_km_yr.json"
    conda:
        "../envs/rasterio.yml"
    params:
        stats='mean sum',
        all_touched=True
    script:
        "../scripts/rasterstats.py"
