# HPGAP
Helminth Population Genomic Analysis Pipeline (`HPGAP`)  overcomes the difficulties in variant calling and is able to explore key aspects centering the population genetics of helminths, including defining population structures, estimating reproductive modes and identifying loci under local adaptation. 

`HPGAP` automates six key bioinformatic components in population genomic analysis using customised R, Perl, Python and Unix shell scripts. These scripts and associated programs were then packaged into a platform independent image to enable setting up an identical system environment.

This repo contains the implementation for `HPGAP`, along with its associated support scripts, libraries and example population genomic data. You can also check the installed programs and packages with the DockerFile.

The workflow for analysis using `HPGAP` consists of five basic parts. 1) Variant calling; 2) Sample filtering; 3) Inferring genetic relationships; 4) Intra-population analysis; 5) Identifying loci under natural selection.



## Dependencies

* git v2.20.1
* wget
* miniconda3
* udocker v1.1.1

## Installation

To install the pipeline, just need to clone the git directory and install udocker and miniconda with the following instructions.

### Installing conda and git

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
~~~

### Installing udocker

Udocker requires either **pycurl** or the **curl** to download both the binaries and/or pull containers from repositories and **tar** to unpackage binaries and libraries. Please check https://github.com/indigo-dc/udocker/blob/master/doc/installation_manual.md for more details of installing udocker.

~~~bash
echo "########################"
echo "# Installing udocker ..."
curl https://download.ncg.ingrid.pt/webdav/udocker/udocker-1.1.2.tar.gz > udocker-tarball.tgz
export UDOCKER_TARBALL=$(pwd)/udocker-tarball.tgz
tar xzvf $UDOCKER_TARBALL udocker
./udocker install
mv ./udocker $HOME/bin/   # move the executable to your preferred location for binaries

wget https://github.com/proot-me/proot-static-build/raw/master/static/proot-x86_64 -P $HOME/bin
wget https://github.com/proot-me/proot-static-build/raw/master/static/proot-x86 -P $HOME/bin
chmod uog+x $HOME/bin/proot*
~~~

### Creating a udocker container

~~~bash
echo "###########################"
echo "# Setting up container ..."
udocker pull walkboy3000/hpgap:latest
udocker create --name=HPGAP_c1 walkboy3000/hpgap:latest
~~~

### Installing SNPEff database


After you create a udocker container and decide the preferred genome, you can run these commands to install the database for SNPEFF into the container and also revise the path of gff file for the following steps. For more details, please check <http://snpeff.sourceforge.net/SnpEff_manual.html#databases>.

* First, use any text editor to create a bash script (e.g. run.sh), then copy, paste and edit the following commands into this file. Pleae make sure the source files (fasta and gff files) are not **soft linked**.

```bash
cd /root/miniconda3/share/snpeff-4.3k-0
#dbname is the preferred name of the SNPEff database.
dbname="Clonorchis_sinensis_henan" # example
mkdir -p "data/$dbname"
cp <the path of genome fasta file> data/$dbname/sequences.fa
cp <the path of gff file> data/$dbname/genes.gff
echo "$dbname.genome : $dbname" >> snpEff.config
snpEff build -gff3 -v $dbname
cd -
```

* Second, run the bash script (e.g. run.sh) with udocker to create the database on the container . Please move your source files (fasta and gff files) to your working directory. If the source files (fasta and gff files are not located in your current working directory, please mount the directory as well). 

```bash
udocker run -v <your host working directory>:<your host working directory> -v <the directory of the source files>:<the directory of the source files> --env="PATH=/usr/local/bin/:/usr/bin:/bin/" HPGAP_c1 /bin/bash -c 'sh run.sh'
```



## Run the pipeline

To run the pipeline, first, you need to prepare the configuration file (see "setting the configuration file"). After that, you can use the following command to generate a "HPGAP.main.sh", which contains the specific command lines for running each step of the analysis. 

**Warning**: HPGAP.main.sh and allcfg.yml should be regenerated if you change the configuration file.

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
		step1_read_report
		step1_recalibration
		step1_variant_calling
		step1_combine_calling
		step1_variant_filtering
		step2_relatedness
		step3_phylogeny
		step3_admixture
		step4_homozygosity
		step4_roh (run of homozygosity)
		step4_ld (linkage disequilibrium)
		step4_slidingwindow
		step4_sfs (snp frequency spectrum)

	--config path to the .yml config file
	
	--step <String>	specified steps separated by semicolon(;). The names of analyses in each step are separated by comma (,);
		(e.g. "0:indexing;1:read_filtering,read_mapping,read_report,recalibration,variant_calling,combine_calling,variant_filtering;3:phylogeny,admixture;4:homozygosity,roh,ld,slidingwindow,sfs").

		All the avaliable analyses in each step: 
			0:indexing;
				  1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;
			2:relatedness
			3:phylogeny,admixture;
			4:homozygosity,roh,ld,slidingwindow,sfs

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



## Required Input

### step1 (from read_filtering to variant calling):
* FASTQ files of each sample
* FASTA file of all the reference sequences. All reference FASTA sequences should be ended with .fa. 

### step2 (relatedness among samples)

* HPGAP generated or custom VCF file (including variant information of all the samples)

### step3 (phylogeny and admixture)

* HPGAP generated or custom VCF file (including variant information of all the samples)

### step4 (homozygosity, roh, ld and sliding window analysis within a population)

* HPGAP generated or custom VCF file (including variant information of all the samples)
* GFF file of the reference genome



## Setting the configuration file

Here is an example configuration file,  you can edit the paths and values within the file and use that for the example analysis.

Warning: YAML file is sensitive to indentis and colon of keys should be followed by one space. Please follow the indents and spaces of the example configuration file carefully.

~~~yaml
---
args:
  container: HPGAP_c1 # Specify the name of the hpgap container
  env: PATH=/root/admixture_linux-1.3.0/:/root/gatk:/root/miniconda3/bin:/usr/local/sbin:/usr/local/bin/:/usr/sbin:/usr/bin:/sbin:/bin:/sbin:/bin:<HPGAP path>:<HPGAP path>/lib:<HPGAP path>/Tools # Set up the $PATH in the container. The paths to git cloned directory and the lib and Tools directories as well.
  mount: # Set up the mapping to mount your host directories to the directories of the container. The directory including the rawdata should be mounted as well. We recommend you to create a local tmp directory (e.g. /home/foobar/tmp) and mount it to the container. 
    -
      host_tmp: <your tmp directory >
    -
      host_path: <HPGAP path>
    -
      host_path: <your working directory>
  outdir: <your output directory> # This directory should be located within the working directory.
  ploidy: '2' # set 1 for haploid; set 2 for diploid.
  threads: 40 # set the number of threads
fqdata: # Set up the sample fqdata
  1-1: # Sample id
    rawdata:
      CL1-1: #library id. When you have more than one library, plese put them together under the same sample ID. 
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
  variant_calling_mode: 'fast' # if fast mode is selected, base recalibration of gatk will be skipped to spead up the variant calling, otherwise, two rounds of base recalibration will be applied.
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



## Structure of output directories

```bash
00.INDEXING
  ├── [Reference genome ID].dict #Dictionary file of the genome
  ├── [Reference genome ID].fa #FASTA file of the genome
  ├── [Reference genome ID].fa.amb #text file recording the "N" in the reference genome 
  ├── [Reference genome ID].fa.ann #text file recording reference sequences, name, length
  ├── [Reference genome ID].fa.bwt #binary file of the Burrows-Wheeler transformed sequence
  ├── [Reference genome ID].fa.fai #fasta index file
  ├── [Reference genome ID].fa.pac #binary file of packaged sequence 
  └── [Reference genome ID].fa.sa #binary file of suffix array index.
01.QualityControl
  ├── Combined
  │   ├── high_confidence_prunned.vcf.gz #VCF file of SNPs with low LD
  │   ├── high_confidence.vcf.gz #VCF file of the SNPs with complete infomation on all the samples
  │   ├── PASS.SNP.DP.vcf.gz #VCF file of filtered SNPs with enough read depth
  │   ├── singletons.list #List of SNP singletons
  │   └── snv.summary.txt #Summary of variants
  ├── read_filtering
  │   ├── [Sample ID]
  │   │   ├── [Library ID + Sample ID]_1.filt.fq.gz # filtered read 1
  │   │   ├── [Library ID + Sample ID]_2.filt.fq.gz # filtered read 2
  ├── read_mapping.[of reference genome ID]
      ├── [Sample ID]
          ├── [Sample ID].sorted.markdup.bam #Duplication marked read mapping BAM file
          └── bam.stats.txt #Summary of the BAM file

03.GeneticRelationships
  ├── Admixture
  │   ├── CV.error.txt # Cross validation error of each K value
  │   ├── high_confidence_prunned.vcf.gz # VCF file used for this analysis
  │   ├── K[K value].png #Admixture graph
  └── Phylogeny
      ├── plink_data.fasta # FASTA file of concatenated SNPs
      └── result2.tre # The output phylogeny

05.IntraPopulation
├── LD  #Linkage disequilibrium
│   ├── South.GLD.png # LD decay graph 
│   ├── South.scf00001.GLD_window_1M.summary #Average LD on each distance
│   └── South.SNP.vcf.gz #VCF file of the SNPs used for this analysis
├── ROH	#Run of Homozygosity
│   ├── South.hom.indiv #Average length of homozygous regions of each individual 
│
├── SFS
│   ├── [Population name]
│   │   ├── EstimatedNe.list #Estimated population size of each run
│   │   ├── [nonsyn_sfs|syn_sfs|total_sfs] 
│   │   │   └── fastsimcoal2
│   │   │       └── South_MSFS.obs #SNP frequency spectrum
│   │   ├── observedFVFiles
│   │   │   └── [Scaffold ID].[Population ID].SNP.preds # predicted loci under natural selection
│   │   ├── trainingSets
│   │   │   ├── hard.fvec #simulated data of hard sweep loci
│   │   │   ├── linkedHard.fvec #simulated data of regions linking to hard sweep loci
│   │   │   ├── linkedSoft.fvec #simulated data of regions linking to soft sweep loci
│   │   │   ├── neut.fvec #simulated data of neutral loci
│   │   │   └── soft.fvec #simulated data of soft sweep loci
└── Slidingwindow
    ├── final.[population ID].allstat.txt # Statistics of each sliding window
    ├── snpEff.nonsyn.vcf.gz #VCF file of nonsynnonymous SNPs
    ├── snpEff.syn.vcf.gz #VCF file of synnonymous SNPs
    ├── snpEff.vcf.gz #functional annotated VCF file
    ├── [population ID].snpEff.nonsyn.pi.list #polymorphism of synnonymous SNPs
    ├── [population ID].snpEff.syn.pi.list #polymorphismfile of nonsynnonymous SNPs


```
