import json
import sys

from rasterstats import zonal_stats

stats = zonal_stats(snakemake.input[0], snakemake.input[1], stats=snakemake.params['stats'],  all_touched=snakemake.params['all_touched'])
with open(str(snakemake.output), 'w') as outfile:
    json.dump(stats, outfile)

sys.exit(0)
