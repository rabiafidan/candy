# Candy lab analysis workflow

1) Mutect2 variant calling
2) Quality filtering
3) Multi sample merging and rescuing
4) annotation to be added

### Install
```
git clone https://github.com/rabiafidan/SRM.git
```

### configure
edit `config.yaml` as you need

### run

## Output files
1) logs/... - output and error files of SLURM (they will be empth if ran locally)
2) Mutect2  - unfiltered calls
3) Mutect2/single - single tumour patients vcf
4) Mutect2/multi - multi tumour patients vcf
