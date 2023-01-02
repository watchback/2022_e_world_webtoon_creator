# **0. Raw data 및 환경설정**
#### 필요한 software  
  
PATH 설정 필요  

- Trimmomatic  

- Trinity  

- Bowtie2  

- Samtools  

- TransDecoder  

- BLAST  

- CD-hit  

- RSEM  

- R  

- R packages (edgeR, ctc, Biobase, ape, gplots 등)  



---
#### RAW DATA
  
FASTQ 파일에는 sequence 서열 뿐만 아니라 quality 정보도 포함하고 있는데  

첫번째 줄은 “@” 로 시작하며 Sequence ID  

두번째 줄은 Sequence 서열  

세번째 줄은 “+” 하나만 있거나, 또는 그 뒤에 첫번째줄의 Sequence ID 의 반복  

네번째 줄은 각 서열의 quality 를 나타내는 기호 (ASCII code)로 구성  
  

파일이 압축되어 있는 경우 압축 해제  

File:**Forward Sequence.gz, Reverse Sequence.gz**  
Tool:**Linux command**  

#### 커맨드

tar.gz 

    tar -zxvf rawdata.fastq_1.tar.gz

gz 

    gzip -d rawdata.fastq_1.gz

#### 결과
**rawdata_1.fastq**  
**rawdata_2.fastq**  


# **1. Preprocessing(Quality Control)**

File: **Forward Sequence.fastq, Reverse Sequence.fastq**  
Tool: **Trimmomatic**  

#### 커맨드 (.sh 파일로 만들어 실행)

    java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq input_reverse.fq output_forward_paired.fq output_forward_unpaired.fq output_reverse_paired.fq output_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#### 옵션
- -threads: thread 수를 지정   
  
- ILLUMINACLIP: fastaWithAdaptersEtc : seed mismatches : palindrome clip threshold : simple clip threshold  
  
fastaWithAdaptersEtc == adaptor 파일, Trimmomatic 프로그램 내에 adaptor 디렉토리에 존재  
  
seed mismatches == 최초 16 bases를 seed로 놓고 이를 full match가 허용하는만큼 확장한다. full match에 허용할 mismatch의 최소값(2)  
  
palindrome clip threshold == paried-ended data의 경우 score 값(30) /약 50 bases  
simple clip threshold == single-ended data의 경우 score 값(10) /약 17 bases  
  
- SLIDINGWINDOW:[num1]:[num2]: 주어진 window(similar to mer) 상수값[num1]만큼 sequence를 sliding하며 window falls [num1]내의 average quality가 주어진 값[num2]보다 낮을 경우 제거   
  
- LEADING (LEADING:3): a read의 앞쪽이 주어진 threshold quality보다 낮은 경우 제거  
  
- TRAILING (TRAILING:3): a read의 뒤쪽이 주어진 threshold quality보다 낮은 경우 제거  
  
- CROP : 명시된 길이만큼 read의 뒤쪽부분을 제거  
  
- HEADCROP : 명시된 길이만큼 read의 앞쪽부분을 제거  
  
- MINLEN:(num) read trimming 과정중에 (num)bp 미만의 read 는 버리라는 명령인데, sequencing data 가 101bp 인 경우에는 MINLEN:36, 151bp 인 경우에는 MINLEN:50을 줌  
  
- TOPHRED33 : quality scores를 Phread-33으로 변경  
    
- TOPHRED64 : quality scores를 Phread-64로 변경  
   

#### 결과
**output_forward_paired.fq**  
**output_reverse_paired.fq**  
**output_forward_unpaired.fq**  
**output_reverse_unpaired.fq**  