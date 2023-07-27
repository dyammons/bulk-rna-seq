## CIBERSORTx installation and basic useage

CIBERSORTx is a deconvolution algo that attempts to prediect the percentage of cell types in a bulk RNA seq sample. The developers offer a webpage tool to run, but the server is not condusive to running/transfering large amounts of data.

---

Here are the steps to get the program running on Alpine:  
1. Visit the [CIBERSORTx website](https://cibersortx.stanford.edu/), create an account and request a token so you are able to run thier software.

<br>

2. We will install in a projects space and use singularity to pull down the docker file for `CIBERSORTxFractions`. Before running the command, please check out the [page on singularity](addlink) if this is your first time using the container software.

```sh
mkdir /projects/$USER/software/cibersort/
singularity pull --dir /projects/$USER/software/cibersort/ docker://cibersortx/fractions
```

<br>

3. Test installation by loading the help menu.
```sh
singularity exec /projects/$USER/software/cibersort/fractions_latest.sif  /src/CIBERSORTxFractions
```

This is just the basics, more to come...
