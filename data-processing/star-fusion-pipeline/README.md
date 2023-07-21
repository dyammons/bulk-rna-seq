## Installation instructions for Alpine

First, navigate to a safe place to store the STAR-fusion singularity image (`/projects/$USER/software/star-fusion/` or similar location). Use `mkdir /projects/$USER/software/star-fusion/` if the directory does not already exist.

<br>

Then migrate to your scratch space to complete the build. Create a working directory `mkdir /scratch/alpine/dyammons@colostate.edu/star-fusion/` then cd into the directory and run:
```sh
wget https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.12.0.simg
```

> Note: for other versions, browse the ftp site directly [here](https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/).

<br>

Next, clone the GitHub repo that contains some essential scripts into the scratch working directory (maybe put a copy in `projects` so you know exactly what code was used).
```sh
git clone https://github.com/NCIP/ctat-genome-lib-builder.git
```


<br>

If you followed the above file structure you can test the install with:
```sh
singularity exec -e -H /projects/$USER/software/star-fusion/ /projects/$USER/software/star-fusion/star-fusion.v1.12.0.simg /projects/$USER/software/star-fusion/ctat-genome-lib-builder/prep_genome_lib.pl -h
```

When running the above command you should see a help menu. If you do, then the software should be installed and ready for use.

<br>

There are a few more steps, but the next we provide a template to build a fusion library to make predictions using the canine genome. This is a big job that can take time, so reference the associated `cute-star-fusion.sbatch` to get a job submitted. Check all paths match up then run
```sh
sbatch cute-star-fusion.sbatch
```

<br>

This page is underdevelopment... more to come.

<br>

##### Extra

By default singularity points to your home directory, you may run into trouble pulling singularity images. To address this potential problem, the code below will reroute the install path.

Create directories to point to:
```sh
mkdir /scratch/alpine/$USER/cache/ /scratch/alpine/$USER/tmp/
```

Then export the paths so singularity can find them:
```sh
export APPTAINER_CACHEDIR=/scratch/alpine/$USER/cache/
export APPTAINER_TMPDIR=/scratch/alpine/$USER/tmp/
export SINGULARITY_CACHEDIR=/scratch/alpine/$USER/cache/
export SINGULARITY_TMPDIR=/scratch/alpine/$USER/tmp/
```

> pro-tip: if you do not want to export the paths every time you use `singularity pull`, you can add the export statments to your `.bashrc` file in `$HOME`.

