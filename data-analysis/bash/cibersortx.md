## CIBERSORTx installation and basic useage

CIBERSORTx is a deconvolution algorithm that attempts to predict the proportions of cell population in bulk RNA seq datasets using a single cell RNA seq (or sorted bulk seq of cell populations) as a reference. The developers offer a webpage tool to run, but the server is not ideal for the evaluataion of ideal reference parameters and has limited functionality if running many samples. Therefore, the instructions provided here details the steps required to obtain a cibersort `docker` container to run the program on Alpine using `Singularity`.

Below are the steps to get the program running on Alpine:  

---

<br>

#### 1) Create a CIBERSORTx account to get access to the software.

Visit the [CIBERSORTx website](https://cibersortx.stanford.edu/), create an account and request a token so you are able to run thier software. 

FYI: it usually takes a week or so to get approval from the cibersort team.

<br>

#### 2) Get the cibersortx software container.  

The cibersortx container is made with docker and the instructions provided by the cibersort team will not work on Alpine because we must run the software through Singularity.  

To set up our computing space we will install the software in a `projects` space directory and use `singularity` to pull down the `docker` file for `CIBERSORTxFractions`. Before running the command, you may want to check out the [page on singularity - currently broken link](addlink) if this is your first time using the container software, as you will likely want to modify the default cache paths.

<br>

```sh
#make directory
mkdir -p /projects/$USER/software/cibersort/

#pull down the container
singularity pull --dir /projects/$USER/software/cibersort/ docker://cibersortx/fractions
```

<br>

#### 3) Test installation by loading the help menu.

Execute the `.sif` file.
```sh
singularity exec /projects/$USER/software/cibersort/fractions_latest.sif  /src/CIBERSORTxFractions
```

<details>
<summary>Output</summary>

<br>

NOTE: the bind mount instructions are not 100% accurate when using `singularity` instead of `docker` to run the code.

```sh
CIBERSORTx - enumerate cell composition in bulk genomic profiles

For instructions and terms of use, see cibersortx.stanford.edu

Usage:
docker run <bind_mounts> cibersortxfractions [Options] 

Bind Mounting:
> 2 directories must be bind mounted to be accessed within the docker container: 
    1) Input file dir 
        > Format: -v {dir_path}:/src/data 
    2) Outdir 
        > Format: -v {dir_path}:/src/outdir 
> Note: Absolute paths required

Authorization Parameters:
--username      <string>  Email used for login to cibersortx.stanford.edu
--token         <string>  Token associated with current IP address (generated on website)

Primary Options:
--mixture       <file_name>  Mixture matrix [required for running CIBERSORTx, optional for creating a custom signature matrix only]
--sigmatrix     <file_name>  Signature matrix [required: use preexisting matrix or create one]
--perm          <int>   No. of permutations for p-value calculation [default: 0]
--label         <char>  Sample label [default: none]
--rmbatchBmode  <bool>  Run B-mode batch correction [default: FALSE]
--rmbatchSmode  <bool>  Run S-mode batch correction [default: FALSE]
--sourceGEPs    <file_name>  Signature matrix GEPs for batch correction [default: sigmatrix]
--QN            <bool>  Run quantile normalization [default: FALSE]
--absolute      <bool>  Run absolute mode [default: FALSE]
--abs_method    <char>  Pick absolute method ['sig.score' (default) or 'no.sumto1']
--verbose       <bool>  Print verbose output to terminal [default: FALSE]

Options for creating a custom signature matrix:
--refsample     <file_name>  Reference profiles (w/replicates) or labeled scRNA-Seq data [required]
--phenoclasses  <file_name>  Cell type classes [required, if single_cell = FALSE]
--single_cell   <bool>  Create matrix from scRNA-Seq data [default: FALSE]
--G.min         <int>   Minimum number of genes per cell type in sig. matrix [default: 50, if single_cell = TRUE: 300]
--G.max         <int>   Maximum number of genes per cell type in sig. matrix [default: 150, if single_cell = TRUE: 500] 
--q.value       <int>   Q-value threshold for differential expression [default: 0.3, if single_cell = TRUE: 0.01] 
--filter        <bool>  Remove non-hematopoietic genes [default: FALSE] 
--k.max         <int>   Maximum condition number [default: 999] 
--remake        <bool>  Remake signature gene matrix [default: False] 
--replicates    <int>   Number of replicates to use for building scRNAseq reference file [default: 5] 
--sampling      <float> Fraction of available single cell GEPs selected using random sampling [default: 0.5] 
--fraction      <float> Fraction of cells of same identity showing evidence of expression [default: 0.75] 
```

</details>

<br>

#### 4) Set up bind mounts.

Next, step is to add in paths to directories that Singularity can mount to give access to the files we need to load in (and save as output) when running `CIBERSORTx`.  

For ease of keeping the code universal, we will run this example in out `scratch` space. If applying the code else where you will need to modify the paths.

<br>

```sh
#make directories for example
mkdir -p /scratch/alpine/$USER/proj_cibersort/src/data
mkdir -p /scratch/alpine/$USER/proj_cibersort/src/outDir
```

<br>

Run command with bind mounts (`-B`) and this should produce the help menu again.
```sh
singularity exec -B /scratch/alpine/$USER/proj_cibersort/src/data:/src/data -B /scratch/alpine/$USER/proj_cibersort/src/outDir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions
```

<br>


