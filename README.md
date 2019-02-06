# HPGAP
Helminth Population Genomic Analysis Pipeline

This repo contains the implementation for `HPGAP` , along with its associated support scripts, libraries and example population genomic data. You can also check the installed programs and packages with the DockerFile.

The workflow for analysis using `HPGAP` consists of five basic parts. 1) Variant calling; 2) Sample filtering; 3) Inferring genetic relationships; 4) Intra-population analysis; 5) identifying loci under nature selection. 

## Dependencies

* miniconda
* udocker

## Installation

To install the pipeline, just need to clone the git directory and install udocker and miniconda with the following instructions.

~~~bash
git clone https://github.com/vetscience/HPGAP.git

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
echo -e "\nyes\n\nyes\n" | bash Miniconda3-latest-Linux-x86_64.sh
export PATH=~/miniconda3/bin:$PATH
conda update -y -n base conda
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

echo "########################"
echo "# Installing udocker ..."
conda install -y -p ~/miniconda3 udocker==1.1.1
rm -f udocker.py
wget https://raw.githubusercontent.com/indigo-dc/udocker/7f6975c19c63c3d65ec6256c7cf5b2369d5c115d/udocker.py
sed 's/proot_killonexit = True/proot_killonexit = False/1' udocker.py > ~/miniconda3/bin/udocker.py
rm -f udocker.py
chmod uog+x ~/miniconda3/bin/udocker.py
sed 's/#!\/bin\/bash/#!\/bin\/bash\nexport CONDA_PREFIX=~\/miniconda3/1' ~/miniconda3/bin/udocker > udocker
chmod uog+x udocker
mv udocker ~/miniconda3/bin
wget https://github.com/proot-me/proot-static-build/raw/master/static/proot-x86_64 -P ~/miniconda3/bin
wget https://github.com/proot-me/proot-static-build/raw/master/static/proot-x86 -P ~/miniconda3/bin
chmod uog+x ~/miniconda3/bin/proot*

echo "###########################"
echo "# Installing containers ..."
udocker pull walkbay3000/hpgap:latest
udocker create --name=HPGAP_c1 walkbay3000/hpgap:latest
~~~



## Run the pipeline

To run the pipeline, simply follow this command:

~~~bash
udocker run -v <the path to your git cloned directory>:<the path to your git cloned directory> -v <the path to your working directory>:<the path to your working directory> --env="PATH=<the path to your git cloned directory>" HPGAP_c1 /bin/bash -c 'HPGAP.pl --config <the path to your configuration file>'
~~~

This command will generate an HPGAP.main.sh in your working directory. You should first try to run the HPGAP.main.sh step by step. 



## Usage

~~~
Usage
	HPGAP.pl --samplelist sample.list --reference reference.fa [-options]

	--run <String> use this option choose one of steps below (the option should not be used with --step at the same time)
		step0_indexing
		step1_read_filtering
		step1_read_mapping
		step1_recalibration
		step1_variant_calling
		step1_combine_calling
		step1_variant_filtering
		step2_relatedness
		step3_phylogeny
		step3_admixture
		step5_homozygosity
		step5_roh
		step5_ld
		step5_slidingwindow
		step5_sfs

	--config path to the .yml config file
	--step <String>	specified steps , separated by commas (e.g., "0:A;1:A,B,C,E,F,G").
		0:indexing;
		1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;
		3:phylogeny,admixture;
		5:homozygosity,roh,ld,slidingwindow,sfs
	--skipsh use this to skip running bash inside the pipeline
	--help

Note
	1. No symbolic link to the files outside the mounted volume, which means all the data files themselves should be located within the mounted volume.
	2. For each pair of fastq data, the second colomn (library or flowcell code) should always be unique.
	3. All the paths need to be absolute path
~~~



## Setting the configuration file

Here is an example configuration file,  you can change the path within the file and run the example analysis.

~~~yaml
---
args:
  container: HPGAP_c1
  env: PATH=/root/admixture_linux-1.3.0/:/root/gatk:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin/:/usr/sbin:/usr/bin:/sbin:/bin:/sbin:/bin:<your path to the pipeline>:<your path to the pipeline>/lib:<your path to the pipeline>/Tools
  mount:
    -
      dockerpath: <your working directory in the container>
      hostpath: <your working directory>
  outdir: <your working directory>
  ploidy: '2'
fqdata:
  1-1:
    rawdata:
      CL1-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/1-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/1-1_r2.fastq.gz
  1-2:
    rawdata:
      CL1-2:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/1-2_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/1-2_r2.fastq.gz
  14-1:
    rawdata:
      CL14-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/14-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/14-1_r2.fastq.gz
  17-1:
    rawdata:
      CL17-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/17-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/17-1_r2.fastq.gz
  26-1:
    rawdata:
      CL26-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/26-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/26-1_r2.fastq.gz
  32-1:
    rawdata:
      CL32-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/32-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/32-1_r2.fastq.gz
  35-1:
    rawdata:
      CL35-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/35-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/35-1_r2.fastq.gz
  44-1:
    rawdata:
      CL44-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/44-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/44-1_r2.fastq.gz
  49-1:
    rawdata:
      CL49-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/49-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/49-1_r2.fastq.gz
  56-1:
    rawdata:
      CL56-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/56-1_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/56-1_r2.fastq.gz
  56-2:
    rawdata:
      CL56-2:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: /home/darcy/PopGen_WorkFlow/Example/Input/Data/56-2_r1.fastq.gz
        fq2: /home/darcy/PopGen_WorkFlow/Example/Input/Data/56-2_r2.fastq.gz
population:
  1-1:
    'presumed population': South
  1-2:
    'presumed population': South
  14-1:
    'presumed population': South
  17-1:
    'presumed population': South
  26-1:
    'presumed population': South
  32-1:
    'presumed population': South
  35-1:
    'presumed population': South
  44-1:
    'presumed population': South
  49-1:
    'presumed population': South
  56-1:
    'presumed population': North
  56-2:
    'presumed population': North
ref:
  choose: Cs-c1
  db:
    Cs-c1:
      name: Cs-c1
      path: /home/darcy/PopGen_WorkFlow/Example//00.INDEXING//Cs-c1.example.fa
    Cs-k2:
      name: Cs-k2
      path: /home/darcy/PopGen_WorkFlow/Example//00.INDEXING//Cs-k2.example.fa
step1:
  variant_filtering:
    indel: 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
    ldcutoff: '0.3'
    ldwindowsize: '50'
    ldwindowstep: '10'
    scaffold_length_cutoff: '0'
    scaffold_number_limit: '2'
    snp: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
step3:
  admixture: ~
step4:
  discoal:
    hard_simulation_times: '100'
    neut_simulation_times: '100'
    soft_simulation_times: '100'
  slidingwindow:
    gff: /home/darcy/PopGen_WorkFlow/Example/Input/Data/clonorchis_sinensis.example.gff
    scaffold_length_cutoff: '5000'
    scaffold_number_limit: '10000'
    snpeff_species: Clonorchis_sinensis_henan
    windowsize: '5000'
~~~

