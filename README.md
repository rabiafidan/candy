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
and change directory to the repo


### configure
edit `config.yaml` as you need

### run
check `snakemake.sh` for usage. Comment/uncomment as you wish.

```
sbatch snakemake.sh
```

## Output files
1) logs/... - output and error files of SLURM (they will be empth if ran locally)
2) Mutect2  - unfiltered calls
3) Mutect2/single - single tumour patients vcf
4) Mutect2/multi - multi tumour patients vcf
