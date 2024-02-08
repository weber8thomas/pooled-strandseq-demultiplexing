#!/bin/bash

elements=(
    "2022-11-15-H33JMAFX5" # HGSVCpool1xulOPxEcho, HGSVCpool1xulOPxmanual
    "2022-11-25-H37MNAFX5" # HGSVCpool1quadrant2KAPA
    # "2023-02-08-HCN3VAFX5" # HGSVCpool2
    # "2023-03-08-HCNGHAFX5" # HGSVCpool2inWell2ul, HGSVCpool2inWell5ul, HGSVCpool2OPSfromFrozen2ul
    "2023-04-21-HGF2CAFX5" # LanexHGSVCpool2500nlEcho
    # "2023-04-26-HCMMNAFX5" # HGSVCpool2iinWell2ulLS, HGSVCpool2OPS500nl
    # "2023-06-23-HGFLGAFX5" # HGSVCpool3UVled
)

# Iterate over each element in the array
for element in "${elements[@]}"; do
    # Split the element into run and plate using '/' as the delimiter
    IFS='/' read -r run plate <<<"$element"

    # Output or use the run and plate variables
    echo "Run: $run"
    snakemake --config bam_folder="/scratch/tweber/DATA/MC_DATA/STOCKS/$run" --profile ../snakemake_profiles/HPC/dev/slurm_legacy_conda/

    # Example of filling in config for snakemake, assuming you want to set variables or create files
    # For demonstration purposes, showing how you might use these in a command:
    # echo "snakemake --config run=$run plate=$plate"
    # Actual implementation would depend on your specific snakemake setup
done
