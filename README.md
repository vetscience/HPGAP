# HPGAP
Helminth Population Genomic Analysis Pipeline

This repo contains the implementation for `HPGAP` , along with its associated support scripts, libraries and example population genomic data. You can also check the installed programs and packages with the DockerFile.

The workflow for analysis using `HPGAP` consists of five basic parts. 1) Variant calling; 2) Sample filtering; 3) Inferring genetic relationships; 4) Intra-population analysis; 5) identifying loci under nature selection. 

## Dependencies

* git
* wget
* miniconda3
* udocker v1.1.1

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

To run the pipeline, first, you need to prepare the configuration file (see "setting the configuration file"). After that, you can use the following command to generate a "HPGAP.main.sh", which contains the specific command lines for running each step of the analysis.

~~~bash
udocker run -v <the path to your git cloned directory>:<the path to your git cloned directory> -v <the path to your working directory>:<the path to your working directory> --env="PATH=<the path to your git cloned directory>" HPGAP_c1 /bin/bash -c 'HPGAP.pl --config <the path to your configuration file>'
~~~

This command will generate an HPGAP.main.sh in your working directory. You should first try to run the HPGAP.main.sh step by step. 



## Usage

~~~
Usage
	HPGAP.pl --config data.yml --run <String> [-options]
	
	--run <String> use this option to choose one of steps below (this option should not be used with --step at the same time)
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

Here is an example configuration file,  you can edit the paths and values within the file and use that for the example analysis.

~~~yaml
---
args:
  container: HPGAP_c1 # Specify the name of the hpgap container
  env: PATH=/root/admixture_linux-1.3.0/:/root/gatk:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin/:/usr/sbin:/usr/bin:/sbin:/bin:/sbin:/bin:<your path to the pipeline>:<your path to the pipeline>/lib:<your path to the pipeline>/Tools # Set up the $PATH in the container. The paths to git cloned directory and the lib and Tools directories as well.
  mount: # Set up the mapping to mount your host directories to the directories of the container 
    -
      dockerpath: <your container pipeline directory>
      hostpath: <your host pipeline directory >
    -
      dockerpath: <your container working directory> 
      hostpath: <your host working directory >
  outdir: <your working directory>
  ploidy: '2' # set 1 for haploid; set 2 for diploid.
fqdata: # Set up the sample fqdata
  1-1: # Sample id
    rawdata:
      CL1-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  1-2: # Sample id
    rawdata:
      CL1-2: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  14-1: # Sample id
    rawdata:
      CL14-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  17-1: # Sample id
    rawdata:
      CL17-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  26-1: # Sample id
    rawdata:
      CL26-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  32-1: # Sample id
    rawdata:
      CL32-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  35-1: # Sample id
    rawdata:
      CL35-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  44-1: # Sample id
    rawdata:
      CL44-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  49-1: # Sample id
    rawdata:
      CL49-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  56-1: # Sample id
    rawdata:
      CL56-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
  56-2: # Sample id
    rawdata:
      CL56-2:
        Flag: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform 
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <your path to read 1 file >
        fq2: <your path to read 2 file >
population: 
  1-1: # Sample id
    'presumed population': South # Here, you can assign population labels for your samples (eg. South and North Populations)
  1-2: # Sample id
    'presumed population': South
  14-1: # Sample id
    'presumed population': South
  17-1: # Sample id
    'presumed population': South
  26-1: # Sample id
    'presumed population': South
  32-1: # Sample id
    'presumed population': South
  35-1: # Sample id
    'presumed population': South
  44-1: # Sample id
    'presumed population': South
  49-1: # Sample id
    'presumed population': South
  56-1: # Sample id
    'presumed population': North
  56-2: # Sample id
    'presumed population': North
ref:
  choose: Cs-c1 # You can choose the reference you prefer for the following analysis
  db:
    Cs-c1: # label of reference 1
      name: Cs-c1 # label of reference 1
      path: /home/darcy/PopGen_WorkFlow/Example//00.INDEXING//Cs-c1.example.fa
    Cs-k2: # label of reference 2
      name: Cs-k2 # label of reference 2
      path: /home/darcy/PopGen_WorkFlow/Example//00.INDEXING//Cs-k2.example.fa
step1: # parameter settings for step1
  variant_filtering:
    indel: 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' # the settins for indel filtering
    ldcutoff: '0.3' # cut off for filtering SNVs with high LD
    ldwindowsize: '50' # window size for for filtering SNVs with high LD
    ldwindowstep: '10' # window step for filtering SNVs with high LD
    scaffold_length_cutoff: '0' # only analyse the scaffolds with length larger than the cutoff
    scaffold_number_limit: '2' # the maximum number of scaffolds
    snp: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' # the settins for snv filtering
step3: # parameter settings for step2
  admixture: ~
step4: # parameter settings for step4
  discoal: 
    hard_simulation_times: '100' # the number of simulation for generating simulated hard sweep data
    neut_simulation_times: '100' # the number of simulation for generating simulated neutral data
    soft_simulation_times: '100' # the number of simulation for generating simulated soft sweep data
  slidingwindow:
    gff: /home/darcy/PopGen_WorkFlow/Example/Input/Data/clonorchis_sinensis.example.gff # the gff file for the species
    scaffold_length_cutoff: '5000' # only analyse the scaffolds with length larger than the cutoff
    scaffold_number_limit: '10000' # the maximum number of scaffolds
    snpeff_species: Clonorchis_sinensis_henan # the species name in the snpeff database 
    windowsize: '5000' # window size
~~~

