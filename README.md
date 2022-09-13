# CRON shell script for running methylation study on Oxford Nanopore data (ONT)

This repo contains some shell scripts which will perform a few analysis on different machine. The main aim is to offload basecalling from the sequencer machine, Promethion, and perform the compute intensive basecalling on available DGX machine with 8x V100 GPU graphics card available on-premise. 
## Analysis

Tabe below summarizes the analysis performed by which machine and which script.

| Analysis | Machine | Script |
| ------ | ------ | ------ |
| Raw data generation (fast5) | Promethion CPU  | MiniKNOW |
| Sync "fast" basecalling fast5 and fastq to network | Promethion CPU  | Promethion_machine/syncFile2hexagon.sh |
| Sync "fast" basecalling Completed file to network | Promethion CPU  | Promethion_machine/syncCompleted2hexagon.sh |
| Guppy basecalling | DGX with 8x V100 - 32GB RAM GPU | DGX_machine/basecalling/cron_Guppy_v6_mongoDB.sh |
| Methylation analysis  | DGX with 8x V100 - 32GB RAM GPU | DGX_machine/mapping/run_methylatedBAM_minimap2.sh |

## Installation
#### Promethion

Create folder eg. CRON_HEXAGON and put the script there. 
```sh
cd /data
mkdir CRON_HEXAGON
cd CRON_HEXAGON 
cp Promethion_machine/syncFile2hexagon.sh /data/CRON_HEXAGON
cp Promethion_machine/syncCompleted2hexagon.sh /data/CRON_HEXAGON
```
Create a cron
```sh
crontab -e
# Include below lines into last part 
# Adjust frequency of performing action accordingly to needed frequency 
# 0 */2 * * * /data/CRON_HEXAGON/syncFile2hexagon.sh
# 0 */3 * * * /data/CRON_HEXAGON/syncCompleted2hexagon.sh
```
#### DGX
Needs following installed
- Miniconda3 [https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html]
If host machine is not setup with CUDA and Docker properly linked
- Docker + CUDA [https://grids-bioinformatics.atlassian.net/wiki/spaces/GRIDSBIOIN/pages/355336193/Ronnin+with+Guppy+gpu+Megalodon]
```sh
bash Miniconda3-latest-Linux-x86_64.sh
# Install mamba to base conda environment
conda install mamba -n base -c conda-forge
# Create environment and install modbam2bed
mamba create -n modbam2bed -c bioconda -c conda-forge -c epi2melabs modbam2bed
# activate 
conda activate modbam2bed
# Install latest samtools=1.15.1
conda install -c bioconda samtools=1.15.1
# Install mosdepth
conda install -c bioconda mosdepth
# Install minimap2
conda install -c bioconda minimap2
```
Create folder eg. methyl_test and put the script there. 
```sh
mkdir -p /raid/scratch/user/methyl_test/mapping
cd /raid/scratch/user/methyl_test
cp DGX_machine/basecalling/cron_Guppy_v6_mongoDB.sh /raid/scratch/user/methyl_test
cp DGX_machine/mapping/run_methylatedBAM_minimap2.sh /raid/scratch/user/methyl_test/mapping
```
 Create a cron
```sh
crontab -e
# Include below lines into last part 
# Adjust frequency of performing action accordingly to needed frequency 
# 0 */2 * * * /raid/scratch/user/methyl_test/cron_Guppy_v6_mongoDB.sh
# 0 */3 * * * /raid/scratch/user/methyl_test/mapping/run_methylatedBAM_minimap2.sh
```

# Description of script
#### Promethion_machine/syncFile2hexagon.sh 
rsync the generated fast5 and fastq from the Promethion machine "fast" mode run to network drive

#### Promethion_machine/syncCompleted2hexagon.sh
- Checks the local folder of Proemthion machine for completed flowcell run signal - final_summary.txt
- rsync the Completed flowcell folder from the Promethion machine to network drive
- Update MongoDB with the completed run. Status -> "To_Basecall"

#### DGX_machine/basecalling/cron_Guppy_v6_mongoDB.sh
- Checks mongoDB for Status:"To_Basecall"
- Get full path of the flowcell directory
- Download fast5 data from network drive to DGX scratch
- Run Guppy 6 basecalling (HAC) on GPU 
- Rsync result back to network
- Set Status: "DONE_BASECALL"

#### DGX_machine/mapping/run_methylatedBAM_minimap2.sh
Requires created environment modbam2bed
- Check mongoDB for Status:"DONE_BASECALL"
- samtools merge bam files from Guppy 6.2.1 with methylation status of reads (Output: bam)
- pipe 3 process = samtools fastq | minimap2 | samtools sort
    - Convert bam to fastq 
    - Minimap2 map to GRCh38
    - Create sorted bam 


## License

MIT

**Free Software, Hell Yeah!**
