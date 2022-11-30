# **0. Raw data 및 환경설정**
#### 필요한 software

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

PATH 설정 필요  

---
#### RAW DATA

File format: **FASTQ.gz**  
File:        **Forward Sequence, Reverse Sequence pair**  

파일이 압축되어 있는 경우 압축 해제

#### 커맨드

tar.gz 

    tar -zxvf rawdata.fastq_1.tar.gz

gz 

    gzip -d rawdata.fastq_1.gz

#### 결과
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

#### 커맨드 (.sh 파일로 만들어 실행)

    java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq input_reverse.fq output_forward_paired.fq output_forward_unpaired.fq output_reverse_paired.fq output_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#### 옵션
- -threads: thread 수를 지정   
- ILLUMINACLIP: fastaWithAdaptersEtc : seed mismatches : palindrome clip threshold : simple clip threshold  
fastaWithAdaptersEtc == adaptor 파일, Trimmomatic 프로그램 내에 adaptor 디렉토리에 존재  
seed mismatches == 최초 16 bases를 seed로 놓고 이를 full match가 허용하는만큼 확장한다. full match에 허용할 mismatch의 최소값(2)  
palindrome clip threshold == paried-ended data의 경우 score 값(30) /약 50 bases  
simple clip threshold == single-ended data의 경우 score 값(10) /약 17 bases  
- SLIDINGWINDOW:[num1]:[num2]: 주어진 window(similar to mer) 상수값[num1]만큼 sequence를 sliding하며 window falls [num1]내의 average quality가 주어진 값[num2]보다 낮을 경우 제거한다.  
- LEADING (LEADING:3): a read의 앞쪽이 주어진 threshold quality보다 낮은 경우 제거한다.  
- TRAILING (TRAILING:3): a read의 뒤쪽이 주어진 threshold quality보다 낮은 경우 제거한다.  
- CROP : 명시된 길이만큼 read의 뒤쪽부분을 제거한다.  
- HEADCROP : 명시된 길이만큼 read의 앞쪽부분을 제거한다.  
- MINLEN:(num) read trimming 과정중에 (num)bp 미만의 read 는 버리라는 명령인데, sequencing data 가 101bp 인 경우에는 MINLEN:36, 151bp 인 경우에는 MINLEN:50을 줌  
- TOPHRED33 : quality scores를 Phread-33으로 변경한다.  
- TOPHRED64 : quality scores를 Phread-64로 변경한다.  

#### 결과
**output_forward_paired.fq**  
**output_reverse_paired.fq**  
**output_forward_unpaired.fq**  
**output_reverse_unpaired.fq**  

# **2. Concatenation of reads**

De novo RNA-seq 분석을 하게 된다면 일반적으로 여러 samples(주로 multiple tissues)을 가지고 시작
genome 이 없다보니 genome 과 유사한 역할을 할 수 있는 transcriptome assembly 제작 필요  
이때 comprehensive 한 transcriptome assembly 를 만들기 위해 여러 samples 의 RNA-seq data를 concatenation  
forward끼리, reverse끼리 concat  

#### 커맨드

    cat Brain_1.fastq, Liver_1.fastq, Testis_1.fastq >> Merged_tissues_1.fastq 
    cat Brain_2.fastq, Liver_2.fastq, Testis_2.fastq >> Merged_tissues_2.fastq

#### 결과
**Merged_tissues_1.fastq**  
**Merged_tissues_2.fastq**  

# **3. *De novo* assembly**
Trinity를 활용하여 assembly 진행

#### 커맨드

    Trinity --seqType fq --left Merged_tissues_1.fastq --right Merged_tissues_2.fastq --output trinity_out --max_memory 100G --CPU 8

#### 옵션
--seqType: reads format 을 지정 (fq: fastq, fa: fasta)  
--left, right: foward, reverse reads 를 지정  
--output: 생성될 output directory (trinity 단어가 들어가야함)  
--max_memory:  assembly 과정중 할당할 최대 memory  
--CPU:  사용할 cpu 수  

#### 결과
**Trinity.fasta**  
**TrinityStats.pl**  

---

완료되면 여러 파일들이 생성되는데 이때  Trinity.fasta 라는 파일이 assembly 된 transcriptome  
이 과정 이후에 statistics 를 구하고 싶으면 Trinity tool 의 util directory 내의 TrinityStats.pl을 이용  

# **4. Gene prediction**
Assembled transcripts 의 protein coding genes 을 찾기 위해 TransDecoder 를 이용  

## 4-1. 100 amino acids 이상의 ORF 포함하는 transcripts 추출

#### 커맨드

    TransDecoder.LongOrfs -t Trinity.fasta

#### 옵션
-t: assembly 를 통해 생성된 Trinity.fasta  

#### 결과
**longest_orfs.pep** (이후 분석에 사용)  

**longest_orfs.gff3**  

**longest_orfs.cds** 등 생성  

## 4-2. Homology-based search

기존에 잘알려진 protein database (Swissprot DB) 에 homology search 를 하여 gene prediction 과정에 활용  
NCBI 의 BLAST 사용  

### 4-2-1. swissprot.fasta 파일을 download 후 blast db 형식으로 제작

#### 커맨드

    ncbi-blast-2.x.x+/bin/makeblastdb -in swissprot.fasta -dbtype prot  

#### 옵션
-in:  db 형식으로 만들고자 하는 fasta 파일  
-dbtype:  fasta 파일이 nucleotide 또는 protein sequence 인지에 따라 (nucleotide: nucl,  protein: prot)  

#### 결과

**swissprot.fasta**  

### 4-2-2. db 파일에 transcript assembly 를 homology search

#### 커맨드

    ncbi-blast-2.x.x+/bin/blastp -query logest_orfs.pep -db swissprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6  

#### 옵션
-query: db에 search 할 sequence 파일  

-db: database  

-max_target_seqs: 해당 query가 matching 되는 sequence 중 몇개를 보여줄지에 대한 parameter  

-outfmt: output의 format에 대한 정보 (0~11 까지 있으며 6 은 tabular format)  

-evalue: blast evalue cutoff  

-num_threads: 사용할 cpu 수  

#### 결과
**blastp.outfmt6**  

### 4-2-3. Gene prediction

이제 blast 결과를 토대로 gene prediction 을 진행  
transcripts.transdecoder_dir 가 있는 directory 에서 진행  

#### 커맨드

    TransDecoder.Predict -t Trinity.fasta --retain_blastp_hits blastp.outfmt6 --cpu 8 --single_best_orf

#### 옵션
-t:  Trinity assembly 파일  

--retain_blastp_hits: (4-2-2) 과정에서 얻은 blastp 결과 파일  

--cpu: 사용할 cpu 수  

--single_best_orf: gene prediction 과정 중 하나의 transcripts 에서 6개의 ORFs 를 탐색했을때 가장 best orf 출력  

#### 결과
**Trinity.fasta.transdecoder.pep** 등 생성  

### 4-2-4. Removing redundant transcripts
(4-2-1)~(4-2-3) 과정을 통해 gene prediction 은 완료되었으나, transcriptome 이다 보니 불가피하게 isoforms 이 존재  
Reference genome 의 역할을 하기 위해 representative protein coding genes 만 있는게 좋으므로, redundant 한 transcripts 를 제거하기 위해 CD-hit 사용  

#### 커맨드

    cdhit -i Trinity.fasta.transdecoder.pep -o Trinity.fasta.transdecoder.pep.cdhit -c 0.99 -T 8

#### 옵션
-i:  TransDecoder 를 통해 geneprediction 이 된 fasta 파일  

-o: otuput   

-c: identity cutoff (0.99: 99% 비슷한 transcripts 의 cluster 중 가장 긴 transcript 선택)  

-T:  사용할 cpu 수  

#### 결과
**Trinity.fasta.transdecoder.pep.cdhit**  
**Trinity.fasta.transdecoder.pep.cdhit.clstr**  

---

.cdhit 파일은 sequence,  .cdhit.clstr 파일에는 cluster 가 어떻게 묶였는지에 대한 정보가 담긺  
얻어진 Trinity.fasta.transdecoder.pep.cdhit 파일을 Non-redundant protein coding sequences (NRCDS) 로 명명하고 진행  

### 4-2-5. NRCDS_Trinity.fasta 파일 생성
NRCDS 파일에 대응되는 nucleotide sequence 가 필요, NRCDS 파일의 sequence id 를 *de novo* assembly 파일인 Trinity.fasta 에 matching 하여 nucleotide sequence 를 얻는다.  

#### 커맨드

    awk '/>/ {print $1}' Trinity.fasta.transdecoder.pep.cdhit > output.txt
    sed 's/.p/ /g' output.txt > output2.txt
    awk '{print $1}' output2.txt > output3.txt
    sed 's/>//g' output3.txt > TrinitypepID.txt
    xargs samtools faidx Trinity.fasta <TrinitypepID.txt> NRCDS_sequence.fasta

#### 옵션
xargs: 대량으로 진행  
file1.fasta <[ID.txt]> file2.fasta: 괄호 안의 ID파일과 일치하는 file 1의 시퀀스 파일을 파싱해서 file2로 저장  


#### 결과
**NRCDS_sequence.fasta**  
