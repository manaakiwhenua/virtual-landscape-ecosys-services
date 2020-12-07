import json
import os
import sys

datastack_dict = {
    "invest_version": "3.8.9",
    "args": {
        "biophysical_table_path": snakemake.params["biophysical_table"],
        "calc_n": True,
        "calc_p": False,
        "dem_path": snakemake.params["dem"],
        "k_param": "2",
        "lulc_path": snakemake.params["lulc"],
        "results_suffix": "",
        "runoff_proxy_path": snakemake.params["runoff_proxy"],
        "subsurface_critical_length_n": "200",
        "subsurface_critical_length_p": "",
        "subsurface_eff_n": "80",
        "subsurface_eff_p": "",
        "threshold_flow_accumulation": "4500",
        "watersheds_path": snakemake.params["watersheds"],
        "workspace_dir": snakemake.params["workspace"]
    },
    "model_name": "natcap.invest.ndr.ndr"
}
with open(str(snakemake.output), 'w') as datastack:
    print(json.dumps(datastack_dict, indent=4, sort_keys=True), file=datastack)

sys.exit(0)
