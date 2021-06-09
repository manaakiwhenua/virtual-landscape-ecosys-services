# Virtual Landscapes Ecosystem Services Workflow

This is a workflow produced for PRJ1968. It takes "virtual landscapes" produced by Tom Etherington and determines ecosystem service indicators for them. Each landscape varies some condition of topography, fragmentation or similar metric. This workflow is used to produce a table of outputs for further statistical analysis that ascertains which input properties are influential for particular ecosystem services.

This could be used for real landscapes, but at present would require some modification.

Contact [Richard Law](lawr@landcareresearch.co.nz) for assistance with this code.

![Snakemake logo](docs/snakemake.png)

## Snakemake

This project relies on the workflow software [Snakemake](https://snakemake.readthedocs.io/en/stable/). Snakemake determines the correct order of execution for the workflow, and can increase processing efficiency by running jobs in parallel where this is possible.

Snakemake also works with the containerisation software [Singularity](https://sylabs.io/singularity/). [Docker](https://www.docker.com/) is a competing containerisation software that is (arguably) easier to set up. Snakemake works with Singularity, not Docker, but Snakemake can convert a Docker image into a Singularity image on the fly; so this workflow specifies Docker images rather than Singularity images.

### Useful Snakemake commands

Command line reference: https://snakemake.readthedocs.io/en/stable/executing/cli.html

In general, to run the whole workflow, run this command from the top-level directory. `all` is a dummy target that sits at the bottom of the workflow, and will cause all (missing) intermediate and ultimate files to be generated or updated.

Unless you just happen to have all the required software installed on your machine, you will want to run the worflow with Singularity images and [Conda](https://docs.conda.io/en/latest/) environments.

```bash
snakemake \
	--jobs $(nproc) \
	--snakefile ./workflow/Snakefile \
	--use-singularity \
	--singularity-args "-B /cifs/Tlinc/Projects\ K-O/MCSL/virtual-landscapes-50km:/virtual-landscapes -B $PWD/results:/results -B $PWD/logs:/logs" \
	--singularity-prefix /home/users/$USER/.singularity \
	--use-conda \
	--conda-frontend conda \
	--rerun-incomplete \
	--keep-going \
	all
```

**This is a highly intertwined workflow where multiple rules rely on the same input files. This can cause issues with concurrent file access, where jobs will fail because they are unable to read an input that does exist. This means that it is likely that you will need to run this command a few times until all outputs are generated, since there is currently no system of file locking that will cause Snakemake to avoid scheduling jobs to run at the same time when they depend on the same input file/s.** Alternatively, copy the input directory to a local directory and use that instead of the networked one.

Alternatively, use `--jobs 1` to avoid concurrent jobs (in exchange for a much longer processing time).

The above command mounts the shared directory of inputs landscapes (`/cifs/Tlinc/Projects\ K-O/MCSL/virtual-landscapes-50km`) to the directory `/virtual-landscapes` inside the container. Because the former is probably different on your machine, you will need to change this to reflect your network drive mapping. It refers to the "T drive" CIFS mount known as "T Lincoln".

Ideally the second volume binding (`-B $PWD/results:/results`) would be (the host-machine-specific equivalent of) `-B /cifs/Tlinc/Projects\ K-O/MCSL/virtual-landscapes-50km/results:/results`, but I have had issues getting Snakemake and/or Singularity to write to that shared directory, so I get it to write locally instead.


<details>
	<summary>Details of this error</summary>
	The error is unclear, but it seems to end up mounting the wrong thing.
	```
	snakemake --use-singularity -p -j1 --singularity-args "-B /cifs/Tlinc/Projects\\ K-O/MCSL/virtual-landscapes-50km:/virtual-landscapes -B /cifs/Tlinc/Projects\\ K-O/MCSL/workflow/results:/results" --singularity-prefix /home/users/$USER/.singularity
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
	    input: /cifs/Tlinc/Projects K-O/MCSL/virtual-landscapes-50km/hills/landclass_t1_c10.tif
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

</details>
</br>
To work around this, I have been writing outputs to a directory on my host machine, and then moving it later to the correct network drive for sharing and posterity.

---

#### Other handy commands

To visualise the workflow as a directed acyclic graph (DAG), i.e. an image that represents the flow of files and processes (this requires the `GraphViz` software):

`snakemake --dag --snakefile ./workflow/Snakefile | dot -Tsvg > dag.svg`

You can also run a linting tool over the Snakefile to check for potential problems while developing new components:

`snakemake --lint --snakefile ./workflow/Snakefile`

You can "clean" the output directory in order to re-run the analysis from a clean slate. This is particularly useful while developing new outputs and making changes.

`snakemake --snakefile ./workflow/Snakefile -j1 clean`

## Configuration

There is a configuration file in `config/config.yaml`. This defines all of the possible values of `landclass`, `landscape` and `topography`, as well as the project `"basepath"` (probably your version of `/cifs/Tlinc/Projects K-O/MCSL`.

The configuration also has options for specifying the pixel height and width of all inputs (assumed to be constant) and the resolution. Together these are used in one part of the workflow related to expressing landscape patch sizes as proportions of the landscape. To adapt this for a non-virtual landscape, the workflow would need to determine these dynamically from inputs, and also know how to correctly deal with null/mask data (e.g. ocean) which aren't valid parts of a landscape. These are not issues that exist with the virtual landscapes.
