# HPGAP
Helminth Population Genomic Analysis Pipeline

This repo contains the implementation for `HPGAP` , along with its associated support scripts, libraries and example population genomic data. You can also check the installed programs and packages with the DockerFile.

The workflow for analysis using `HPGAP` consists of five basic parts. 1) Variant calling; 2) Sample filtering; 3) Inferring genetic relationships; 4) Intra-population analysis; 5) identifying loci under nature selection. 

## Dependencies

* git v2.20.1
* wget
* miniconda3
* udocker v1.1.1

## Installation

To install the pipeline, just need to clone the git directory and install udocker and miniconda with the following instructions.

~~~bash

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
echo -e "\nyes\n\nyes\n" | bash Miniconda3-latest-Linux-x86_64.sh
export PATH=~/miniconda3/bin:$PATH
conda update -y -n base conda
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

echo "########################"
echo "# Installing git ..."
conda install -c conda-forge git=2.20.1
git clone https://github.com/vetscience/HPGAP.git # We will refer to the path that you clone the git repository as HPGAP path (eg. /home/foobar/HPGAP)

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
udocker create --name=HPGAP_c1 walkboy3000/hpgap:latest
~~~



## Run the pipeline

To run the pipeline, first, you need to prepare the configuration file (see "setting the configuration file"). After that, you can use the following command to generate a "HPGAP.main.sh", which contains the specific command lines for running each step of the analysis.

~~~bash
udocker run -v <HPGAP path>:<HPGAP path> -v <your host working directory>:<your host working directory> --env="PATH=<HPGAP path>:/usr/local/bin/:/usr/bin:/bin/" HPGAP_c1 /bin/bash -c 'HPGAP.pl --config <the path to your configuration file>'
~~~

This command will generate an HPGAP.main.sh in your working directory. You should first try to run the HPGAP.main.sh step by step. 



## Usage

~~~
Usage
	HPGAP.pl --config <path to the .yml config file> --run <String> [-options]
	
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
		step4_homozygosity
		step4_roh
		step4_ld
		step4_slidingwindow
		step4_sfs

	--config path to the .yml config file
	
	--step <String>	specified steps separated by semicolon(;). The names of analyses in each step are separated by comma (,);
		(e.g. "0:indexing;1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;3:phylogeny,admixture;4:homozygosity,roh,ld,slidingwindow,sfs").

		All the avaliable analyses in each step: 
			0:indexing;
			1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;
			3:phylogeny,admixture;
			4:homozygosity,roh,ld,slidingwindow,sfs
			6:mktest

	--skipsh use this to skip running bash inside the pipeline
	
	--help

Note 
	1. No symbolic link to the files outside the mounted volume, which means all the data files themselves should be located within the mounted volume.
	2. For each pair of fastq data, the second colomn (library or flowcell code) should always be unique.
	3. All the paths need to be absolute path

Example 
	To be done
USAGE
~~~



## Setting the configuration file

Here is an example configuration file,  you can edit the paths and values within the file and use that for the example analysis.

~~~yaml
---
args:
  container: HPGAP_c1 # Specify the name of the hpgap container
  env: PATH=/root/admixture_linux-1.3.0/:/root/gatk:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin/:/usr/sbin:/usr/bin:/sbin:/bin:/sbin:/bin:<HPGAP path>:<HPGAP path>/lib:<HPGAP path>/Tools # Set up the $PATH in the container. The paths to git cloned directory and the lib and Tools directories as well.
  mount: # Set up the mapping to mount your host directories to the directories of the container. We recommend you to create a local tmp directory (e.g. /home/foobar/tmp) and mount it to the container.
    -
      host_tmp: <your tmp directory >
    -
      host_path: <HPGAP path>
    -
      host_path: <your working directory>
  outdir: <your output directory> # This directory should be located within the working directory.
  ploidy: '2' # set 1 for haploid; set 2 for diploid.
fqdata: # Set up the sample fqdata
  1-1: # Sample id
    rawdata:
      CL1-1: #library id
        Flag: PE # PE (paired end) or SE (single end)
        PL: BGISEQ500 # Sequencing platform
        Phred: '33' # Phred scoring system of the fastq data
        fq1: <HPGAP path>/Example/Input/Data/1-1_r1.fastq.gz # You need to edit the path here and make it points to the fastq file (e.g. /home/foobar/Example/Input/Data/1-1_r1.fastq.gz)
        fq2: <HPGAP path>/Example/Input/Data/1-1_r2.fastq.gz
  1-2:
    rawdata:
      CL1-2:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/1-2_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/1-2_r2.fastq.gz
  14-1:
    rawdata:
      CL14-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/14-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/14-1_r2.fastq.gz
  17-1:
    rawdata:
      CL17-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/17-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/17-1_r2.fastq.gz
  26-1:
    rawdata:
      CL26-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/26-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/26-1_r2.fastq.gz
  32-1:
    rawdata:
      CL32-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/32-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/32-1_r2.fastq.gz
  35-1:
    rawdata:
      CL35-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/35-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/35-1_r2.fastq.gz
  44-1:
    rawdata:
      CL44-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/44-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/44-1_r2.fastq.gz
  49-1:
    rawdata:
      CL49-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/49-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/49-1_r2.fastq.gz
  56-1:
    rawdata:
      CL56-1:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/56-1_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/56-1_r2.fastq.gz
  56-2:
    rawdata:
      CL56-2:
        Flag: PE
        PL: BGISEQ500
        Phred: '33'
        fq1: <HPGAP path>/Example/Input/Data/56-2_r1.fastq.gz
        fq2: <HPGAP path>/Example/Input/Data/56-2_r2.fastq.gz
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
      path: <HPGAP path>/Example/Input/Cs-c1.example.fa
    Cs-k2: # label of reference 2
      name: Cs-k2 # label of reference 2
      path: <HPGAP path>/Example/Input/Cs-k2.example.fa
step1: # parameter settings for step1
  variant_filtering:
    indel: 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' # The settings for indel filtering. To modify that, please check the manual of GATK4
    ldcutoff: '0.3' # cut off for filtering SNVs with high LD
    ldwindowsize: '50' # window size for for filtering SNVs with high LD
    ldwindowstep: '10' # window step for filtering SNVs with high LD
    scaffold_length_cutoff: '0' # only analyse the scaffolds with length larger than the cutoff
    scaffold_number_limit: '2' # the maximum number of scaffolds
    snp: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' # The settings for snv filtering. To modify that, please check the manual of GATK4
step4: # parameter settings for step4
  discoal: 
    hard_simulation_times: '100' # the number of simulation for generating simulated hard sweep data
    neut_simulation_times: '100' # the number of simulation for generating simulated neutral data
    soft_simulation_times: '100' # the number of simulation for generating simulated soft sweep data
  slidingwindow:
    gff: <HPGAP path>/Example/Input/clonorchis_sinensis.example.gff # the gff file for the species
    scaffold_length_cutoff: '5000' # only analyse the scaffolds with length larger than the cutoff
    scaffold_number_limit: '10000' # the maximum number of scaffolds
    snpeff_species: Clonorchis_sinensis_henan # the species name in the snpeff database 
    windowsize: '5000' # window size for scanning the whole genome in a sliding window manner
~~~

