![](https://github.com/sinnamone/InteractomeSeq/blob/master/log.png.png)
# InteractomeSeq CLI pipeline with Prokaryotic dataset

### 1) Download Prokaryotic dataset


```bash
pwd
```

    /home/spuccio


#### Activate conda env


```bash
conda activate interactomeseq
```

    (interactomeseq) 



#### Create Working Folder


```bash
mkdir -p InteractomSeqCLIProkaryotic
```


```bash
cd InteractomSeqCLIProkaryotic
```


```bash
WorkingDir=$(pwd)
```


```bash
processor=10
```


```bash
echo $WorkingDir
```

    /home/spuccio/InteractomSeqCLIProkaryotic


#### Data download with wget 


```bash
time wget http://interactomeseq.ba.itb.cnr.it/services/InteractomeSeq/Prokaryote.Uploading/Prokaryote_Hp26695.zip -P $WorkingDir
```

    --2020-04-09 17:09:54--  http://interactomeseq.ba.itb.cnr.it/services/InteractomeSeq/Prokaryote.Uploading/Prokaryote_Hp26695.zip
    Risoluzione di interactomeseq.ba.itb.cnr.it (interactomeseq.ba.itb.cnr.it)... 150.145.82.8
    Connessione a interactomeseq.ba.itb.cnr.it (interactomeseq.ba.itb.cnr.it)|150.145.82.8|:80... connesso.
    Richiesta HTTP inviata, in attesa di risposta... 200 OK
    Lunghezza: 1200780141 (1,1G) [application/zip]
    Salvataggio in: "/home/spuccio/InteractomSeqCLIProkaryotic/Prokaryote_Hp26695.zip"
    
    Prokaryote_Hp26695. 100%[===================>]   1,12G  5,81MB/s    in 4m 4s   
    
    2020-04-09 17:13:58 (4,69 MB/s) - "/home/spuccio/InteractomSeqCLIProkaryotic/Prokaryote_Hp26695.zip" salvato [1200780141/1200780141]
    
    
    real	4m4,398s
    user	0m0,648s
    sys	0m6,992s



```bash
cd $WorkingDir
ls 
```




```bash
mkdir rawdata
```

Unzip folder


```bash
unzip -d rawdata Prokaryote_Hp26695.zip 
```

    Archive:  Prokaryote_Hp26695.zip
      inflating: rawdata/GCF_000008525.1_ASM852v1_genomic.fna  
      inflating: rawdata/GCF_000008525.1_ASM852v1_genomic.gff  
      inflating: rawdata/Sel_AG_R1.fastq.gz  
      inflating: rawdata/Sel_AG_R2.fastq.gz  
      inflating: rawdata/Sel_HPneg_R1.fastq.gz  
      inflating: rawdata/Sel_HPneg_R2.fastq.gz  
      inflating: rawdata/Sel_HPpos_R1.fastq.gz  
      inflating: rawdata/Sel_HPpos_R2.fastq.gz  
      inflating: rawdata/HP_genomic_26695_R1.fastq.gz  
      inflating: rawdata/HP_genomic_26695_R2.fastq.gz  
      inflating: rawdata/Prokaryote_Hp26695_Tutorial.pdf  



```bash
cd rawdata
```


```bash
gunzip *.fastq.gz 
```

Alternatively they can be downloaded from here: https://interactomeseq.ba.itb.cnr.it/Prokaryote#!/Hp26695/Uploading


```bash
RawData=$(pwd)
```


```bash
ls
```

    GCF_000008525.1_ASM852v1_genomic.fna  Sel_AG_R2.fastq
    GCF_000008525.1_ASM852v1_genomic.gff  Sel_HPneg_R1.fastq
    HP_genomic_26695_R1.fastq             Sel_HPneg_R2.fastq
    HP_genomic_26695_R2.fastq             Sel_HPpos_R1.fastq
    Prokaryote_Hp26695_Tutorial.pdf       Sel_HPpos_R2.fastq
    Sel_AG_R1.fastq


#### Now we are ready to start the analysis with InteractomeSeqCLI


```bash
cd $WorkingDir
```


```bash
ls
```

   



```bash
git clone https://github.com/sinnamone/InteractomeSeq
```

    Cloning into 'InteractomeSeq'...
    remote: Enumerating objects: 34, done.
    remote: Counting objects: 100% (34/34), done.
    remote: Compressing objects: 100% (26/26), done.
    remote: Total 581 (delta 13), reused 22 (delta 8), pack-reused 547[K
    Ricezione degli oggetti: 100% (581/581), 263.09 KiB | 0 bytes/s, done.
    Risoluzione dei delta: 100% (423/423), done.



```bash
ls 
```

    

Reference files (FASTA and GFF) are downloaded from https://www.ncbi.nlm.nih.gov/genome/ , before to launch the script fasta index should be indexed while gff should be parsed. For parsing we provide a custom script call gff_parser_ncbi_newformat.py


```bash
samtools faidx  "$RawData"/GCF_000008525.1_ASM852v1_genomic.fna 
```

    (interactomeseq) 




```bash
python "$WorkingDir"/InteractomeSeq/PyinteraseqProkaryotic/gff_parser_ncbi_newformat.py --gff "$RawData"/GCF_000008525.1_ASM852v1_genomic.gff --outputfolder "$RawData" --outputid GCF_000008525.1_ASM852v1_genomic_parsed 
```

    (interactomeseq) 




```bash
ls -t "$RawData" | head
```

    Sel_AG_R2.fastq
    Sel_AG_R1.fastq
    Sel_HPneg_R2.fastq
    Sel_HPneg_R1.fastq
    HP_genomic_26695_R1.fastq
    HP_genomic_26695_R2.fastq
    Sel_HPpos_R2.fastq
    Sel_HPpos_R1.fastq
    GCF_000008525.1_ASM852v1_genomic_parsed.bed
    GCF_000008525.1_ASM852v1_genomic.fna.fai
    (interactomeseq) 



### 2) Mapping reads against reference genome


```bash
mkdir -p output
```


```bash
cd "$WorkingDir"/InteractomeSeq/PyinteraseqProkaryotic
```

#### During this step Adapters are removed with Cutdapt and cleaned reads assigned to the reference genome with Blastn. Input dataset could be executed both with SingleEnd and PairedEnd modality, using Fasta and Fastq. 
#### The main purpose of this first mapping is to filter the reads with many mismatches that could be sequencing artifacts.
#### Blast output are filtered and only the reads with the best match are considered. 


```bash
time for i in HP_genomic_26695 Sel_AG Sel_HPneg Sel_HPpos; do python pyinteraseq_main_mapping.py --readforward "$RawData"/"$i"_R1.fastq --readreverse "$RawData"/"$i"_R2.fastq  --outputfolder "$WorkingDir"/output --primer5forward GCAGCAAGCGGCGCGCATGCCACTAGTGGGAT --primer3forward CCCAGAGCAA  --primer5reverse GGGATTGGTTTGCCGCTAGCGGAGAT --primer3reverse CCCAGAGCAA --fastasequence "$RawData"/GCF_000008525.1_ASM852v1_genomic.fna --thread $processor --outputid $i ; done
```

    
    real	14m42,317s
    user	44m33,984s
    sys	1m15,116s
    (interactomeseq) 



#### Output files are two:
#### 1) File (.tab) with the aligned sequences; 
#### 2) File (.log) with the log of the execution and metrics;
#### File tab will have mutuated from Blast, outformat7 with this columns '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'


```bash
ls "$WorkingDir"/output/*_mapping.* 
```

    /home/spuccio/InteractomSeqCLIProkaryotic/output/HP_genomic_26695_mapping.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/HP_genomic_26695_mapping.tab
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_mapping.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_mapping.tab
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPneg_mapping.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPneg_mapping.tab
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_mapping.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_mapping.tab
    (interactomeseq) 




```bash
head "$WorkingDir"/output/HP_genomic_26695_mapping.log 
```

    Thu Apr  9 21:22:30 2020
    
    ###Setting Parameters
    
    ###Sequencing input:	Paired-End
    ###Input dataset forward:	HP_genomic_26695_R1.fastq
    ###Input dataset reverse:	HP_genomic_26695_R2.fastq
    ###Primer 5' forward read:	GCAGCAAGCGGCGCGCATGCCACTAGTGGGAT
    ###Primer 3' forward read:	CCCAGAGCAA
    ###Primer 5' reverse read:	GGGATTGGTTTGCCGCTAGCGGAGAT
    (interactomeseq) 



### 3) Domain Definition

#### This is the core of the analysis, in this step reads are clustered using sumaclust algoritm. The representative sequence of each cluster, represented by the most_abundant sequence of the cluster are picked using pick_rep_seq.py from qiime (--rep_set_picking_method most_abundant ). The set of representative sequences is mapped against the Input Fasta Sequence. 


```bash
time  for i in HP_genomic_26695 Sel_AG Sel_HPneg Sel_HPpos; do python pyinteraseq_main_domain_definition.py --mappingoutput "$WorkingDir"/output/"$i"_mapping.tab  --outputfolder "$WorkingDir"/output --outputid $i --fastasequence "$RawData"/GCF_000008525.1_ASM852v1_genomic.fna --annotation "$RawData"/GCF_000008525.1_ASM852v1_genomic_parsed.bed --thread $processor ; done
```

    
    real	5m11,850s
    user	37m41,560s
    sys	0m16,640s
    (interactomeseq) 




```bash
ls "$WorkingDir"/output/*_definition.* 
```

    /home/spuccio/InteractomSeqCLIProkaryotic/output/HP_genomic_26695_definition.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/HP_genomic_26695_definition.txt
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_definition.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_definition.txt
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPneg_definition.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPneg_definition.txt
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_definition.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_definition.txt
    (interactomeseq) 



#### Output files are two:
#### 1) File (.txt) with the predicted genomic domains; 
#### 2) File (.log) with the log of the execution and metrics;

#### Output tabular have the following columns : "Chr     CloneStart      CloneEnd        CloneLength     Start   End     GeneID  Strand  Description     NuclSeq"

### 4) Domain Enrichment

#### This step, like the following ones in points 4,5,6, is optional. After the domain definition, the user obtains the list of domains present in the genomic library or selection.

#### The following points have been designed to allow the user to compare / subtract / intersect and calculate a differential statistic of this predicted genomic domains.

#### This depends on the experimental design and serves to assist biologists in interpretation.


```bash
time for i in Sel_AG Sel_HPneg Sel_HPpos; do python pyinteraseq_enrichment_prokaryotic.py --blastnoutputgenomic "$WorkingDir"/output/HP_genomic_26695_mapping.tab --blastnoutputarget "$WorkingDir"/output/"$i"_mapping.tab --domainstarget "$WorkingDir"/output/"$i"_definition.txt --outputfolder  "$WorkingDir"/output/ --outputid "$i"; done 
```

    
    real	0m27,569s
    user	0m26,140s
    sys	0m5,256s
    (interactomeseq) 




```bash
ls "$WorkingDir"/output/*_enrichment.*
```

    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_enrichment.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_enrichment.txt
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPneg_enrichment.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPneg_enrichment.txt
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_enrichment.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_enrichment.txt
    (interactomeseq) 



#### Output files are two:
#### 1) File (.txr) with the statistically enriched domains;
#### 2) File (.log) with the log of the execution and metrics;

#### Output tabular have the following columns: Chr	CloneStart	CloneEnd	CloneLength	Start	End	GeneID	logFC	PValue	AdjPValue	Strand	Description	NuclSeq

### 5) Domain Subtraction

#### This step allows the user to subtract common domains through the use of bedtools intersect.


```bash
time for i in Sel_HPpos Sel_AG ; do python pyinteraseq_subtraction_prokaryotic.py --enrichedcontrol "$WorkingDir"/output/Sel_HPneg_enrichment.txt --enrichedselection "$WorkingDir"/output/"$i"_enrichment.txt --outputfolder "$WorkingDir"/output --outputid $i; done
```

    
    real	0m0,889s
    user	0m0,712s
    sys	0m1,228s
    (interactomeseq) 




```bash
ls -t "$WorkingDir"/output/*_subtraction.* | head
```

    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_subtraction.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_AG_subtraction.txt
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_subtraction.log
    /home/spuccio/InteractomSeqCLIProkaryotic/output/Sel_HPpos_subtraction.txt
    (interactomeseq) 



#### Output files are two:
#### 1) File (.txt) with the subtrated domains;
#### 2) File (.log) with the log of the execution and metrics;

##### Output tabular have the following columns: Chr	CloneStart	CloneEnd	CloneLength	Start	End	GeneID	logFC	PValue	AdjPValue	Strand	Description	NuclSeq

### 6) Domain Intersection

#### This step allows the user to intersect the lists with the predicted genomic domains. This scripts can be used to intersect at most 3 samples and produce a Venn diagram in png format. 


```bash
python pyinteraseq_intersection_prokaryotic.py --selection "$WorkingDir"/output/Sel_HPpos_subtraction.txt --selection "$WorkingDir"/output/Sel_AG_subtraction.txt --label HpPositiveControl --label AtrophicGastritis --outputfolder "$WorkingDir"/output  --outputid HpPositiveControlAtrophicGastritis
```

    (interactomeseq) 




```bash
ls -t "$WorkingDir"/output | head 
```

    HpPositiveControlAtrophicGastritis_AtrophicGastritis_unique.txt
    HpPositiveControlAtrophicGastritis_HpPositiveControl_AtrophicGastritis_intersection.txt
    HpPositiveControlAtrophicGastritis_HpPositiveControl_unique.txt
    HpPositiveControlAtrophicGastritis.png
    Sel_AG_subtraction.log
    Sel_AG_subtraction.txt
    Sel_HPpos_subtraction.log
    Sel_HPpos_subtraction.txt
    Sel_HPpos_enrichment.log
    Sel_HPpos_enrichment.txt
    (interactomeseq) 



#### Output files are two:
#### 1) File (.txr) with the intersected domains unique and common between samples ;¶
#### 2) File (.log) with the log of the execution and metrics;
#### 3) File (.png) Venn diagram plot 

##### Output tabular have the following columns: Chr	CloneStart	CloneEnd	CloneLength	Start	End	GeneID	logFC	PValue	AdjPValue	Strand	Description	NuclSeq


```bash
echo "Done"
```

    Done
    (interactomeseq) 


