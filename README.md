Run this command from the top-level directory

`snakemake --cores=$(nproc) -p --jobs $(nproc)`

To generate the DAG (requires `GraphViz`):

`snakemake --dag | dot -Tsvg > dag.svg`

Run a linting tool over the Snakefile:

`snakemake --lint`

---

## Singularity/Docker

Run this instead:

```bash
snakemake --use-singularity -p -j1 --singularity-args "-B /cifs/Tlinc/Projects\ K-O/MCSL/virtual-landscapes:/virtual-landscapes -B $PWD/results:/results -B $PWD/logs:/logs" --singularity-prefix /home/users/$USER/.singularity --use-conda
```

Ideally the second volume binding (`-B $PWD/results:/results`) would be `-B /cifs/Tlinc/Projects\ K-O/MCSL/virtual-landscapes/results:/results`, but for some reason Snakemake and/or Singularity are having trouble writing there, or something... the error is not clear, it seems to end up mounting the wrong thing. The whitespace in the path is not helpful!! Error:

```
snakemake --use-singularity -p -j1 --singularity-args "-B /cifs/Tlinc/Projects\\ K-O/MCSL/virtual-landscapes:/virtual-landscapes -B /cifs/Tlinc/Projects\\ K-O/MCSL/workflow/results:/results" --singularity-prefix /home/users/$USER/.singularity
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	4	carbon_stock
	4	greenhouse_gas_emission
	9

[Wed Dec  2 12:43:28 2020]
rule carbon_stock:
    input: /cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes/hills/landclass_t1_c10.tif
    output: results/carbon_stocks/hills/t1.c10.LandscapeCarbonStock.tif
    jobid: 4
    wildcards: landscape=hills, topography=t1, landclass=c10

gdal_calc.py -A "/virtual-landscapes/hills/landclass_t1_c10.tif" --outfile=/results/carbon_stocks/hills/t1.c10.LandscapeCarbonStock.tif --calc="(A==0)*5 + (A==1)*10 + (A==2)*18 + (A==3)*30 + (A==4)*140 + (A==5)*200" --type=UInt16 --overwrite --co TILED=YES
Activating singularity image /home/users/lawr/.singularity/d656e5e84caaafdc61745ab1ef68999a.simg
FATAL:   could not open image /home/users/lawr/Network/Tlinc/Projects K-O/MCSL/workflow/K-O/MCSL/workflow: failed to retrieve path for /home/users/lawr/Network/Tlinc/Projects K-O/MCSL/workflow/K-O/MCSL/workflow: lstat /cifs/Tlinc/Projects K-O/MCSL/workflow/K-O: no such file or directory
[Wed Dec  2 12:43:28 2020]
Error in rule carbon_stock:
    jobid: 4
    output: results/carbon_stocks/hills/t1.c10.LandscapeCarbonStock.tif
    shell:
        gdal_calc.py -A "/virtual-landscapes/hills/landclass_t1_c10.tif" --outfile=/results/carbon_stocks/hills/t1.c10.LandscapeCarbonStock.tif --calc="(A==0)*5 + (A==1)*10 + (A==2)*18 + (A==3)*30 + (A==4)*140 + (A==5)*200" --type=UInt16 --overwrite --co TILED=YES
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: /cifs/Tlinc/Projects K-O/MCSL/workflow/.snakemake/log/2020-12-02T124327.441032.snakemake.log
```

Most relevant part:

```
FATAL:   could not open image /home/users/lawr/Network/Tlinc/Projects K-O/MCSL/workflow/K-O/MCSL/workflow: failed to retrieve path for /home/users/lawr/Network/Tlinc/Projects K-O/MCSL/workflow/K-O/MCSL/workflow: lstat /cifs/Tlinc/Projects K-O/MCSL/workflow/K-O: no such file or directory
```

In the meantime, the workflow can be run from a host machine, writing output to a local drive, and then moved later. This obviously has the user's host drive as an upper limit.

---

## Configuration

There is a configuration file in `config/config.yaml`. This defines all of the possible values of `landclass`, `landscape` and `topography`, as well as the project `"basepath"` (assumed `/cifs/Tlinc/Projects K-O/MCSL`, **but this is only my local setup and is not universally applicable**).
