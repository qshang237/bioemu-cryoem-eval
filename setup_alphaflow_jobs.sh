#!/bin/bash
# Run this from ~/alphaflow on the HPC login node:
#   bash setup_alphaflow_jobs.sh
set -e
cd "$(dirname "$0")"

# 1. Split alphaflow_input.csv into one CSV per remaining protein
for p in GPR4 GltPh SLC37A4 SPNS2; do
    echo "name,seqres" > "alphaflow_input_${p}.csv"
    grep "^${p}," alphaflow_input.csv >> "alphaflow_input_${p}.csv"
done

# 2. Generate one sbatch script per protein
for p in GPR4 GltPh SLC37A4 SPNS2; do
cat > "run_alphaflow_${p}.sh" << EOF
#!/bin/bash
#SBATCH -A pilot_sae_gpu
#SBATCH -p sae
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie:1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=alphaflow_${p}
#SBATCH --output=/data/home/bty189/alphaflow/alphaflow_${p}_%j.out

module load miniforge
module load cuda/11.8.0-gcc-12.2.0
mamba activate alphaflow

cd /data/home/bty189/alphaflow

python predict.py \\
    --mode alphafold \\
    --input_csv alphaflow_input_${p}.csv \\
    --msa_dir msa_dir/ \\
    --weights params/alphaflow_md_base_202402.pt \\
    --samples 500 \\
    --outpdb alphaflow_results/

echo "Done"
EOF
done

# 3. Submit all 4
for p in GPR4 GltPh SLC37A4 SPNS2; do
    sbatch "run_alphaflow_${p}.sh"
done

echo "---"
squeue -u "$USER"
