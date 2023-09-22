## kallisto pipeline installation and useage


### Softare installation


### Index generation

Retrieve the cdna (transcriptome from ensembl or other source)
```sh
curl -O ftp://ftp.ensembl.org/pub/release-104/fasta/canis_lupus_familiaris/cdna/Canis_lupus_familiaris.CanFam3.1.cdna.all.fa.gz
```

### index generation
kallisto index -i canFam3.1.fa.idx Canis_lupus_familiaris.CanFam3.1.cdna.all.fa
