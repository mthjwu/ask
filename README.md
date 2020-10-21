# ask #
--------
Seek amplicons from high throughput sequencing data 

### How to download? ###
--------
Download the package using git clone or using the "download" link.

     git clone https://github.com/mthjwu/ask

### Software dependencies ###
--------
The software has been tested in MacOSX and Linux system.
The software does not depend on any other softwares except some basic python packages.
Prequired python packages: pysam, pandas, numpy, statsmodels, matplotlib, seaborn

### How to install python ###
--------
* install Miniconda by following https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
        bash Miniconda3-latest-Linux-x86_64.sh
* create a environment with the prequired python packages
        conda create -n ask --no-channel-priority pysam, pandas, numpy, statsmodels, matplotlib, seaborn
* activate environment
        conda activate ask
Now, you are ready to run ask

### How to run from bam file? ###
--------
* run ask from sorted, deduplicated bam files with index file in the same folder (see below for how to generate such bam file)
        <ask_dir>/ask_cmd.py -i test.bam -o test_ask/test_ask

### How to prepare bam file? ###
--------
* Map fastq file to the genome
        # paired end
        bwa mem -t 5 <bwa_index> test_R1.fastq.gz test_R2.fastq.gz | samtools view -Shb - > test_unsorted.bam
        # single end
        bwa mem -t 5 <bwa_index> test.fastq.gz | samtools view -Shb - > test_unsorted.bam
* sort and mark duplicates
        samtools fixmate --threads 5 -m test_unsorted.bam - \
        |samtools sort --threads 5 -T ./ - \
        |samtools markdup --threads 5 -T ./ -S -s - test.bam
* make index
        samtools index test.bam

