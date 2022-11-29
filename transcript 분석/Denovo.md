# **0. Raw data 및 환경설정**

File format: FASTQ.gz  
File:        Forward Sequence, Reverse Sequence pair  

파일이 압축되어 있는 경우 압축 해제

### 커맨드

tar.gz 
    tar -xf rawdata.fastq_1.tar.gz
gz 
    gzip -r rawdata.fastq_1.gz

### 결과
rawdata_1.fastq  
fawdata_2.fastq  

# **1. Preprocessing(Quality Control)**

File format: FASTQ  
File:        Forward Sequence, Reverse Sequence pair  
Tool: Trimmomatic  

### 커맨드 (.sh 파일로 만들어 실행)

    java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq input_reverse.fq output_forward_paired.fq output_forward_unpaired.fq output_reverse_paired.fq output_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### 옵션

### 결과
output_forward_paired.fq  
output_reverse_paired.fq  
output_forward_unpaired.fq  
output_reverse_unpaired.fq  

# **2. Concatenation of reads (reference genome 역할을 하는 cDNA 제작)**

De novo RNA-seq 분석을 하게 된다면 일반적으로 여러 samples(주로 multiple tissues)을 가지고 시작하게 됩니다.  
또한 genome 이 없다보니 genome 과 유사한 역할을 할 수 있는 transcriptome assembly를 만들어야 합니다.  

> 2-1. 이때 comprehensive 한 transcriptome assembly 를 만들기 위해 여러 samples 의 RNA-seq data를 concatenation 해야합니다. 
forward끼리, reverse끼리 합쳐줌  

### 커맨드

    cat Brain_1.fastq, Liver_1.fastq, Testis_1.fastq >> Merged_tissues_1.fastq 
    cat Brain_2.fastq, Liver_2.fastq, Testis_2.fastq >> Merged_tissues_2.fastq


> 2-2. De novo assembly
Tool: Trinity  

### 커맨드
Trinity --seqType fq --left Merged_tissues_1.fastq --right Merged_tissues_2.fastq --output trinity_out --max_memory 100G --CPU 8

### 옵션
--seqType:  reads format 을 지정합니다. (fq: fastq, fa: fasta)  
--left, right: foward, reverse reads 를 지정합니다.  
--output: 생성될 output directory 입니다. (trinity 단어가 들어가야합니다.)  
--max_memory:  assembly 과정중 할당할 최대 memory  입니다.  
--CPU:  사용할 cpu 수 입니다.  

### 결과
Trinity.fasta  
TrinityStats.pl  

완료되면 여러 파일들이 생성되는데 이때  Trinity.fasta 라는 파일이 assembly 된 transcriptome 입니다.  
이 과정 이후에 statistics 를 구하고 싶으면 Trinity tool 의 util directory 내의 TrinityStats.pl 을 통해 할 수 있다.  

> 2-3. Gene prediction