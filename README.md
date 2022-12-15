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

## Input
You need bam files and recalibrated.tsv/csv from SAREK. You will enter this in `config.yaml`. 

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
3) Mutect2/single - single tumour patients vcfs
4) Mutect2/multi - multi tumour patients vcfs

## Notes
You can add new lines to the recalibrated.tsv/csv it should only run the new ones (95% sure but if it does not you might need to add `--rerun-trigger mtime` to the running script `snakemake.sh`) (Please let me know how it goes if you test this so that I can update here)

If you add new samples of an existing patient, the patient VCF should now appear under Mutect2/multi with all listed (old and new) samples of the patient but it would not delete the old Mutect2/single/patient.vcf file, you should delete it manually. (If this is a big inconvinience, let me know)
