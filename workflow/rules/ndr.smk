rule nutrient_delivery_ratio_datastack:
    input:
        "results/precipitation/{landscape}/{topography}.precipitation.tif",
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/static/invest/Biophysical_Table.csv",
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/dem_{topography}.tif",
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif",
        "/cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/static/VirtualDomain/VirtualDomain.geojson",
    output:
        "results/ndr/{landscape}/{topography}.{landclass}.ndr.invs.json"
    params:
        workspace="results/ndr/{landscape}/{topography}.{landclass}",
        watersheds="/virtual-landscapes/static/VirtualDomain/VirtualDomain.shp",
        lulc="/virtual-landscapes/{landscape}/landclass_{topography}_{landclass}.tif",
        dem="/virtual-landscapes/{landscape}/dem_{topography}.tif",
        biophysical_table="/virtual-landscapes/static/invest/Biophysical_Table.csv",
        runoff_proxy="results/precipitation/{landscape}/{topography}.precipitation.tif"
    container:
        "docker://richardlaw/natcap-invest:latest"
    script:
        "../scripts/invest-ndr.py"

rule nutrient_delivery_ratio:
    input:
        "results/ndr/{landscape}/{topography}.{landclass}.ndr.invs.json"
    params:
        workspace="results/ndr/{landscape}/{topography}.{landclass}",
    container:
        "docker://richardlaw/natcap-invest:latest"
    output:
        "results/ndr/test.{landscape}.{topography}.{landclass}.txt"
    shell:
        'invest run ndr --headless --datastack {input} -w {params.workspace}'
