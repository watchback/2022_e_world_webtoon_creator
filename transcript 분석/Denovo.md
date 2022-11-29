# **0. Raw data 및 환경설정**
### 필요한 sortware

- Trinity  

- Bowtie2  

- Samtools  

- TransDecoder  

- BLAST  

- CD-hit  

- RSEM  

- R  

- R packages (edgeR, ctc, Biobase, ape, gplots 등)  

PATH 설정 필요  

---

File format: **FASTQ.gz**  
File:        **Forward Sequence, Reverse Sequence pair**  

파일이 압축되어 있는 경우 압축 해제

### 커맨드

tar.gz 

    tar -zxvf rawdata.fastq_1.tar.gz

gz 

    gzip -d rawdata.fastq_1.gz

### 결과
**rawdata_1.fastq**  
**rawdata_2.fastq**  

---

FASTQ 파일에는 sequence 서열 뿐만 아니라 quality 정보도 포함하고 있는데  

첫번째 줄은 “@” 로 시작하며 Sequence ID  

두번째 줄은 Sequence 서열  

세번째 줄은 “+” 하나만 있거나, 또는 그 뒤에 첫번째줄의 Sequence ID 의 반복  

네번째 줄은 각 서열의 quality 를 나타내는 기호 (ASCII code)로 구성 

 

# **1. Preprocessing(Quality Control)**

File format: **FASTQ**  
File: **Forward Sequence, Reverse Sequence pair**  
Tool: **Trimmomatic**  

### 커맨드 (.sh 파일로 만들어 실행)

    java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq input_reverse.fq output_forward_paired.fq output_forward_unpaired.fq output_reverse_paired.fq output_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### 옵션
-threads: thread 수를 지정   
ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>  
<fastaWithAdaptersEtc> == adaptor 파일, Trimmomatic 프로그램 내에 adaptor 디렉토리에 존재  
<seed mismatches> == 최초 16 bases를 seed로 놓고 이를 full match가 허용하는만큼 확장한다. full match에 허용할 mismatch의 최소값(2)  
<palindrome clip threshold> == paried-ended data의 경우 score 값(30) /약 50 bases  
<simple clip threshold> == single-ended data의 경우 score 값(10) /약 17 bases  

SLIDINGWINDOW (SLIDINGWINDOW:[num1]:[num2]) : 주어진 window(similar to mer) 상수값[num1]만큼 sequence를 sliding하며 window falls [num1]내의 average quality가 주어진 값[num2]보다 낮을 경우 제거한다.  
LEADING (LEADING:3): a read의 앞쪽이 주어진 threshold quality보다 낮은 경우 제거한다.  
TRAILING (TRAILING:3): a read의 뒤쪽이 주어진 threshold quality보다 낮은 경우 제거한다.  
CROP : 명시된 길이만큼 read의 뒤쪽부분을 제거한다.  
HEADCROP : 명시된 길이만큼 read의 앞쪽부분을 제거한다.  
MINLEN:(num) read trimming 과정중에 (num)bp 미만의 read 는 버리라는 명령인데, sequencing data 가 101bp 인 경우에는 MINLEN:36, 151bp 인 경우에는 MINLEN:50을 줌  
TOPHRED33 : quality scores를 Phread-33으로 변경한다.  
TOPHRED64 : quality scores를 Phread-64로 변경한다.  

### 결과
**output_forward_paired.fq**  
**output_reverse_paired.fq**  
**output_forward_unpaired.fq**  
**output_reverse_unpaired.fq**  

# **2. Concatenation of reads (reference genome 역할을 하는 cDNA 제작)**

De novo RNA-seq 분석을 하게 된다면 일반적으로 여러 samples(주로 multiple tissues)을 가지고 시작
genome 이 없다보니 genome 과 유사한 역할을 할 수 있는 transcriptome assembly 제작 필요  

## 2-1. 이때 comprehensive 한 transcriptome assembly 를 만들기 위해 여러 samples 의 RNA-seq data를 concatenation 해야합니다. 

forward끼리, reverse끼리 합쳐줌  

### 커맨드

    cat Brain_1.fastq, Liver_1.fastq, Testis_1.fastq >> Merged_tissues_1.fastq 
    cat Brain_2.fastq, Liver_2.fastq, Testis_2.fastq >> Merged_tissues_2.fastq

### 결과
**Merged_tissues_1.fastq**  
**Merged_tissues_2.fastq**  

## 2-2. *De novo* assembly
Trinity를 활용하여 assembly 진행

### 커맨드

    Trinity --seqType fq --left Merged_tissues_1.fastq --right Merged_tissues_2.fastq --output trinity_out --max_memory 100G --CPU 8

### 옵션
--seqType:  reads format 을 지정 (fq: fastq, fa: fasta)  
--left, right: foward, reverse reads 를 지정  
--output: 생성될 output directory (trinity 단어가 들어가야함)  
--max_memory:  assembly 과정중 할당할 최대 memory  
--CPU:  사용할 cpu 수  

### 결과
**Trinity.fasta**  
**TrinityStats.pl**  

---

완료되면 여러 파일들이 생성되는데 이때  Trinity.fasta 라는 파일이 assembly 된 transcriptome  
이 과정 이후에 statistics 를 구하고 싶으면 Trinity tool 의 util directory 내의 TrinityStats.pl을 이용  

## 2-3. Gene prediction