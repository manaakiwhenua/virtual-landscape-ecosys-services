Run this command from the top-level directory

`snakemake --cores=$(nproc) -j --jobs $(nproc)`

To generate the DAG (requires `GraphViz`):

`snakemake --dag | dot -Tsvg > dag.svg`
