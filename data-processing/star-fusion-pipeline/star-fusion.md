In the directory with your samples run:
```sh
ls | grep ".fq" | cut -d"_" -f1 | sort -u > samples.tmp
```
This will generate a .tmp file that can the be read and used to generate a metadata file. In the same input directory run (be sure to modify the `ext` variable if needed):

```sh
ext=".fastq.gz" #set this to match you file extensions

while read line
do
     echo $line"_1"$ext","$line"_2"$ext","$line >> metadata.csv
done < samples.tmp
```


Now you can use the metadata file to loop through you samples.  
You can run the following in command line or put it in an executable bash script and submit as a job.
```sh
while read line
do

  sample1=( $(echo $line | cut -f1 -d',' --output-delimiter=' ') )
  sample2=( $(echo $line | cut -f2 -d',' --output-delimiter=' ') )
  name=( $(echo $line | cut -f3 -d',' --output-delimiter=' ') )

  singularity exec -e -B `pwd` -B /path/to/ctat_genome_lib_build_dir \ ###set paths
  star-fusion.v1.12.0.simg \ ###check version
  STAR-Fusion \
  --left_fq $sample1 \ ###this is pulled from metadata
  --right_fq $sample2 \ ###this is pulled from metadata
  --genome_lib_dir /path/to/ctat_genome_lib_build_dir \ # set paths
  -O $name \ #this is pulled from metadata
  --FusionInspector validate \
  --examine_coding_effect \
  --denovo_reconstruct

done < /pwd/to/metadata.csv ###set to point to metadata file
```

        
