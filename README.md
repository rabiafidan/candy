Automation of the workflow originally developped by Irene Lobon

# Candy lab analysis workflow

1) Mutect2 variant calling
2) Quality filtering
3) Multi sample merging and rescuing
4) .. TO BE ADDED: annotation

### Install
```
git clone --branch pipeline https://github.com/rabiafidan/SRM.git
```

### configure
edit `config.yaml` as you need

### run
check `snakemake.sh` for the pipeline in cluster mode.
run it locally:

```
snakemake --latency-wait 100 --rerun-incomplete -p --cluster --keep-going --use-conda --use-envmodules
```

## Output files
1) logs/... - output and error files of SLURM (they will be empth if ran locally)
2) Mutect2  - unfiltered calls
3) Mutect2/single - single tumour patients vcf
4) Mutect2/multi - multi tumour patients vcf
