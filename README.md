# dsp_cnv_exp1

##

Log into `Freya` through PuTTy. Port 22. Default settings. Enter `Freya`


## Data

```
/home/qumulo/NextSeq2000/201121_VH00121_143_AAAH23NM5
```

- Made directory in `NAS_data`

`mkdir /home/NAS_data/zlewis`


## Using `DNDC`

```
# cd to directory with fastq files
scp -r -i ~/.ssh/prod-key-efs-dndc.pem DSP*.fastq.gz ec2-user@44.241.187.245:/mnt/efs

# copy over .ini file
scp -r -i ~/.ssh/prod-key-efs-dndc.pem *.ini ec2-user@44.241.187.245:/mnt/efs

```


## Generating deduplicated counts from the DCC files

- Edit the hard paths in [dnd_summary.R](DSP_CNV_e1/dnd_summary.R)
- Run the `R` script
