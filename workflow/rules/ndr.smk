rule nutrient_delivery_ratio_datastack:
    input:
        "results/precipitation/{landscape}/{topography}.precipitation.tif",
        f"{config['basepath']}"+"/static/invest/Biophysical_Table.csv",
        f"{config['basepath']}"+"/{landscape}/{landscape}-dem_{topography}.tif",
        f"{config['basepath']}"+"/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
        f"{config['basepath']}"+"/static/VirtualDomain/VirtualDomain.v2.geojson",
    output:
        "results/ndr/{landscape}/{topography}.{landclass}.{proportion}.ndr.invs.json"
    params:
        workspace=lambda w, output: str(output).rstrip(".ndr.invs.json"),
        watersheds="/virtual-landscapes/static/VirtualDomain/VirtualDomain.v2.shp",
        lulc="/virtual-landscapes/{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}.tif",
        dem="/virtual-landscapes/{landscape}/{landscape}-dem_{topography}.tif",
        biophysical_table="/virtual-landscapes/static/invest/Biophysical_Table.csv",
        runoff_proxy=lambda w, input: input[0]
    container:
        "docker://richardlaw/natcap-invest:latest"
    script:
        "../scripts/invest-ndr.py"

rule nutrient_delivery_ratio:
    input:
        "results/ndr/{landscape}/{topography}.{landclass}.{proportion}.ndr.invs.json"
    params:
        workspace=lambda w, input: str(input).rstrip(".ndr.invs.json")
    container:
        "docker://richardlaw/natcap-invest:latest"
    output:
        "results/ndr/{landscape}/{topography}.{landclass}.{proportion}/n_export.tif",
        "results/ndr/{landscape}/{topography}.{landclass}.{proportion}/watershed_results_ndr.shp"
    log:
        "logs/ndr/{landscape}/{topography}.{landclass}.{proportion}.n_export.log"
    shell:
        'invest run ndr --headless --datastack {input} -w {params.workspace} > {log} 2>&1'
