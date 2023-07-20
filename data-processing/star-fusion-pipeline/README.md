## Installation instructions for Alpine

First, navigate to a safe place to install the software (`/projects/$USER/software/star-fusion/` or similar location). Use `mkdir /projects/$USER/software/star-fusion/` if the directory does not already exist.

<br>

In the desired working directory run:
```sh
wget https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/star-fusion.v1.12.0.simg
```

> Note: for other versions, browse the ftp site directly [here](https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/).

<br>

Next, clone the GitHub repo that contains some essential scripts
```sh
git clone https://github.com/NCIP/ctat-genome-lib-builder.git
```

<br>

From here, all of the required external software should be avalible locally. Now it is time to set up singulaity, as the container will be used to run the code. By default singularity points to your home directory, the code below will reroute the install path.

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

pro-tip: if you do not want to export the paths every time you use singularity, you can add the export statments to your `.bashrc` file in $HOME.

<br>

If you followed the above file structure you can test the install with:
```sh
singularity exec -e -H /projects/$USER/software/star-fusion/ /projects/$USER/software/star-fusion/star-fusion.v1.12.0.simg /projects/$USER/software/star-fusion/ctat-genome-lib-builder/prep_genome_lib.pl -h
```

When running the above command you should see a help menu. If you do, then the software should be installed and ready for use.

<br>

There are a few more steps, but the next will provide a template to build a fusion library to make predictions using the canine genome. This is a big job that can take time, so reference the associated `cute-star-fusion.sbatch` and `star-fusion.sh` scripts to get a job submitted.
```sh
singularity exec -e -H /projects/$USER/software/star-fusion/ /projects/$USER/software/star-fusion/star-fusion.v1.12.0.simg /projects/$USER/software/STAR-Fusion/star-fusion/ctat-genome-lib-builder/prep_genome_lib.pl \
--CPU 8 --genome_fa /absolute/path/to/Canis_familiaris.CanFam3.1.dna.toplevel.fa \
--gtf /absolute/path/to/Canis_familiaris.CanFam3.1.98.gtf \
--pfam_db /absolute/path/to/Pfam-A.hmm.gz \
--dfam_db /absolute/path/to/Dfam.hmm.gz
```
