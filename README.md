# nanoDoc

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Introduction

RNA modification detection using Nanopore raw reads with Deep One Class classification

This software requires In vitro unmodified RNA raw read and 
Native RNA raw read to process.

GPU enviroment is strongly recommended to run this software.

You can find more infromation about nanoDoc in this preprint:

nanoDoc: RNA modification detection using Nanopore raw reads with Deep One-Class Classification

https://www.biorxiv.org/content/10.1101/2020.09.13.295089v1


## Dependency

### Python
Python (>= 3.6), packages,tensorflow, checked with tensorflow 2.3, cuda/10.1, cudnn/7.6


## Install

```
$ git clone https://github.com/uedaLabR/nanoDoc.git  
$ cd nanoDoc/src
$ python3 -m venv venv3
$ source venv3/bin/activate
(venv3) $ pip install --upgrade pip
(venv3) $ pip install -r requirements.txt
(venv3) $ cd..
(venv3) $ mkdir weight5mer
(venv3) $ mv./weight5mer_1/* ./weight5mer
(venv3) $ mv./weight5mer_2/* ./weight5mer
(venv3) $ rm -r /weight5mer_*

```


## Prepareation

Basecalling and signal alignment is required prior to run this program.

Tombo (https://github.com/nanoporetech/tombo) resiggle command is used for propressing.


## Commands

  #### formatFile:   create uniform bin sized parque file from tombo resiggled signle fast5 files.
  
  python ./nanoDoc.py formatFile -i fast5dir -o outputdir -r fast_genome_reference -t thread
  
  
  fast5dir - directory contains fast5 files (files will be searched under the directory recursively)   
  outputdir - oututdir  
  fast_genome_reference - reference fasta file  
  thread - number of threads (defult 4), this process is slow. a large number of thread (e.g. 10) is recommended if resources allowed.  
  
  e.g.
  ```
  python ./nanoDoc.py formatfile -i /mydir/testIVT/singleFast5 -o /mydir/testIVTout -r /reference/NC000913.fa -t 10
  ```
  
  #### analysis:  analyse modification sites, given IVT and Native raw reads sequence
  
  python ./nanoDoc.py analysis -w 5-mer wight -p parameter_file -rraw dir_to_IVT_raw_parquet_file -traw dir_to_Native_raw_parquet_file -o result_output
                               -chrom chromosome(defult first chromosome in the reference) -s start -e end -strand strand(defult "+")  
  
  5-mer wight - please download prebuild weight for each 5mer  
  parameter_file - please download parameter file for nanoDoc  
  dir_to_IVT_raw_parquet_file - directory containing reference Invitro parquet files, created by formatFile command  
  dir_to_Native_raw_parquet_file - directory containing target Native files, created by formatFile command  
  result_output - outout text format file.  
  chromosome - chromosome, which have to match to reference file  
  start - start position (defult 0)  
  end - end position (defult end of the reference)  
  strand - strand "+" or "- (defult "+")  
  e.g.
  ```
  python ./nanoDoc.py analysis -w /weight5mer/ -p /param20.txt -r /reference/NC000913.fa -rraw /equalbinnedpq/ecrRnaIvt -traw /equalbinnedpq/ecrRnaNative -o /result/result.txt -s 4035631 -e 4037072
```  
please download the preculculate weight from repository weight5mer_1/weight5mer_2 and merge inner directory to /weight5mer
  
  
## Commands for preparetion (calculate weight)

#### make 5mer parquet for training
python ./nanoDocPrep.py make5mer -r /path/to/reference.fa -rraw /path/to/ivt/parquetfile/testIVTout -ssize 12000 -o /path/to/out/fivemerparquet

#### training first cnn
python ./nanoDocPrep.py traincnn -in /path/to/out/fivemerparquet -o /path/to/out/weight

#### DOC training
python ./nanoDocPrep.py traindoc -in /path/to/out/fivemerparquet -wight /path/to/out/weight/xxx.hdf -o /path/to/out/fivemer/weight






