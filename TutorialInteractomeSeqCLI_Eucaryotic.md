![](https://github.com/sinnamone/InteractomeSeq/blob/master/log.png.png)
### 1) Data download 

### Remenber to activate conda interactomeseq env 


```bash
conda activate interactomeseq
```

    (interactomeseq) 




```bash
pwd -P
```


```bash
mkdir -p InteractomeSeqCLIEukaryotic
```

    (interactomeseq) 




```bash
cd InteractomeSeqCLIEukaryotic
```

    (interactomeseq) 




```bash
WorkingDir=$(pwd)
```

    (interactomeseq) 




```bash
processor=10
```

    (interactomeseq) 




```bash
echo $WorkingDir
```

    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic
    (interactomeseq) 




```bash
time wget http://interactomeseq.ba.itb.cnr.it/services/InteractomeSeq/Eukaryote.Uploading/Eukaryote_RnaBindProt.zip -P $WorkingDir
```

    --2020-04-10 14:15:35--  http://interactomeseq.ba.itb.cnr.it/services/InteractomeSeq/Eukaryote.Uploading/Eukaryote_RnaBindProt.zip
    Risoluzione di interactomeseq.ba.itb.cnr.it (interactomeseq.ba.itb.cnr.it)... 150.145.82.8
    Connessione a interactomeseq.ba.itb.cnr.it (interactomeseq.ba.itb.cnr.it)|150.145.82.8|:80... connesso.
    Richiesta HTTP inviata, in attesa di risposta... 200 OK
    Lunghezza: non specificato [text/plain]
    Salvataggio in: "/home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/Eukaryote_RnaBindProt.zip"
    
    Eukaryote_RnaBindPr     [               <=>  ] 105,76M  6,03MB/s    in 17s     
    
    2020-04-10 14:15:52 (6,36 MB/s) - "/home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/Eukaryote_RnaBindProt.zip" salvato [110902666]
    
    
    real	0m16,682s
    user	0m0,092s
    sys	0m0,428s
    (interactomeseq) 




```bash
time wget interactomeseq.ba.itb.cnr.it/services/InteractomeSeq/Eukaryote.Uploading/Homo_sapiens.GRCh38.99.zip -P $WorkingDir
```

    --2020-04-10 14:16:54--  http://interactomeseq.ba.itb.cnr.it/services/InteractomeSeq/Eukaryote.Uploading/Homo_sapiens.GRCh38.99.zip
    Risoluzione di interactomeseq.ba.itb.cnr.it (interactomeseq.ba.itb.cnr.it)... 150.145.82.8
    Connessione a interactomeseq.ba.itb.cnr.it (interactomeseq.ba.itb.cnr.it)|150.145.82.8|:80... connesso.
    Richiesta HTTP inviata, in attesa di risposta... 200 OK
    Lunghezza: non specificato [text/plain]
    Salvataggio in: "/home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/Homo_sapiens.GRCh38.99.zip"
    
    Homo_sapiens.GRCh38     [     <=>            ]   1,90G  3,57MB/s    in 6m 46s  
    
    2020-04-10 14:23:40 (4,81 MB/s) - "/home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/Homo_sapiens.GRCh38.99.zip" salvato [2043656730]
    
    
    real	6m45,648s
    user	0m1,576s
    sys	0m10,396s
    (interactomeseq) 




```bash
cd $WorkingDir
ls 
```

    (interactomeseq) [0m[01;31mEukaryote_RnaBindProt.zip[0m  [01;31mHomo_sapiens.GRCh38.99.zip[0m
    (interactomeseq) 




```bash
mkdir -p rawdata
```

    (interactomeseq) 




```bash
mkdir -p references
```

    (interactomeseq) 




```bash
time unzip -d rawdata Eukaryote_RnaBindProt.zip
```

    Archive:  Eukaryote_RnaBindProt.zip
      inflating: rawdata/Delta5Uchl1ssDNA.fastq.gz  
      inflating: rawdata/Delta5Uchl1tRNA.fastq.gz  
      inflating: rawdata/invSINEB2ssDNA.fastq.gz  
      inflating: rawdata/invSINEB2tRNA.fastq.gz  
      inflating: rawdata/NotSelected.fastq.gz  
      inflating: rawdata/Eukaryote_RnaBindProt_Tutorial.pdf  
    
    real	0m0,834s
    user	0m0,632s
    sys	0m0,080s
    (interactomeseq) 




```bash
time unzip -d references Homo_sapiens.GRCh38.99.zip
```

    Archive:  Homo_sapiens.GRCh38.99.zip
      inflating: references/Homo_sapiens.GRCh38.99.bed  
      inflating: references/Homo_sapiens.GRCh38.99.cdna.all.fa  
      inflating: references/Homo_sapiens.GRCh38.99.sizes.genome  
      inflating: references/Homo_sapiens.GRCh38.99.chr.gtf  
      inflating: references/Homo_sapiens.GRCh38.99.idx  
      inflating: references/Homo_sapiens.GRCh38.99.sizes.transcriptome  
    
    real	0m33,296s
    user	0m30,424s
    sys	0m2,856s
    (interactomeseq) 




```bash
cd rawdata
```

    (interactomeseq) 




```bash
time gunzip *.fastq.gz 
```

    
    real	0m6,333s
    user	0m5,548s
    sys	0m0,776s
    (interactomeseq) 




```bash
ls
```

    Delta5Uchl1ssDNA.fastq  Eukaryote_RnaBindProt_Tutorial.pdf  invSINEB2tRNA.fastq
    Delta5Uchl1tRNA.fastq   invSINEB2ssDNA.fastq                NotSelected.fastq
    (interactomeseq) 




```bash
pwd
```

    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/rawdata
    (interactomeseq) 




```bash
RawData=$(pwd)
```

    (interactomeseq) 




```bash
cd $WorkingDir
```

    (interactomeseq) 




```bash
cd references
```

    (interactomeseq) 




```bash
Reference=$(pwd)
```

    (interactomeseq) 



#### Now we are ready to start the analysis with InteractomeSeqCLI


```bash
cd $WorkingDir
```

    (interactomeseq) 




```bash
git clone https://github.com/sinnamone/InteractomeSeq
```

    Cloning into 'InteractomeSeq'...
    remote: Enumerating objects: 61, done.[K
    remote: Counting objects: 100% (61/61), done.[K
    remote: Compressing objects: 100% (48/48), done.[K
    remote: Total 608 (delta 24), reused 41 (delta 13), pack-reused 547[K
    Ricezione degli oggetti: 100% (608/608), 304.66 KiB | 0 bytes/s, done.
    Risoluzione dei delta: 100% (434/434), done.
    (interactomeseq) 




```bash
samtools faidx  "$Reference"/Homo_sapiens.GRCh38.99.cdna.all.fa
```

    (interactomeseq) 



### 2) Mapping reads against reference Homo Sapiens genome


```bash
cd "$WorkingDir"
```

    (interactomeseq) 




```bash
ls 
```

    [0m[01;31mEukaryote_RnaBindProt.zip[0m   [01;34mInteractomeSeq[0m  [01;34mreferences[0m
    [01;31mHomo_sapiens.GRCh38.99.zip[0m  [01;34mrawdata[0m
    (interactomeseq) 




```bash
mkdir -p output
```

    (interactomeseq) 




```bash
cd "$WorkingDir"/InteractomeSeq/PyinteraseqEukaryotic
```

    (interactomeseq) 




```bash
time for i in NotSelected invSINEB2tRNA invSINEB2ssDNA Delta5Uchl1tRNA  Delta5Uchl1ssDNA; do python pyinteraseq_main_mapping.py --readforward "$RawData"/"$i".fastq --outputfolder  "$WorkingDir"/output/ --primer5forward GCAGCAAGCGGCGCGCATGC --primer3forward GCGCTTCGTCAT --fastasequence  "$Reference"/Homo_sapiens.GRCh38.99.cdna.all.fa --thread $processor --outputid "$i" --kallistoindex "$Reference"/Homo_sapiens.GRCh38.99.idx --genomesize "$Reference"/Homo_sapiens.GRCh38.99.sizes.genome --transcriptsize "$Reference"/Homo_sapiens.GRCh38.99.sizes.transcriptome --gtf "$Reference"/Homo_sapiens.GRCh38.99.chr.gtf ; done
```

    
    real	11m51,442s
    user	10m34,396s
    sys	1m15,256s
    (interactomeseq) 




```bash
ls "$WorkingDir"/output/*_mapping.*
```

    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1ssDNA_mapping.bam
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1ssDNA_mapping.bam.bai
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1ssDNA_mapping.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNA_mapping.bam
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNA_mapping.bam.bai
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNA_mapping.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2ssDNA_mapping.bam
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2ssDNA_mapping.bam.bai
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2ssDNA_mapping.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNA_mapping.bam
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNA_mapping.bam.bai
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNA_mapping.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/NotSelected_mapping.bam
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/NotSelected_mapping.bam.bai
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/NotSelected_mapping.log
    (interactomeseq) 



### 3) Domain Definition


```bash
time for i in Delta5Uchl1ssDNA Delta5Uchl1tRNA invSINEB2ssDNA invSINEB2tRNA NotSelected ; do python pyinteraseq_main_domain_definition.py --mappingoutput "$WorkingDir"/output/"$i"_mapping.bam --outputfolder "$WorkingDir"/output/ --outputid "$i" --fastasequence "$Reference"/Homo_sapiens.GRCh38.99.cdna.all.fa --annotation "$Reference"/Homo_sapiens.GRCh38.99.bed --threshold 30 --genome "$Reference"/Homo_sapiens.GRCh38.99.sizes.transcriptome ; done 
```

    sys:1: DtypeWarning: Columns (0) have mixed types. Specify dtype option on import or set low_memory=False.
    sys:1: DtypeWarning: Columns (0) have mixed types. Specify dtype option on import or set low_memory=False.
    sys:1: DtypeWarning: Columns (0) have mixed types. Specify dtype option on import or set low_memory=False.
    sys:1: DtypeWarning: Columns (0) have mixed types. Specify dtype option on import or set low_memory=False.
    sys:1: DtypeWarning: Columns (0) have mixed types. Specify dtype option on import or set low_memory=False.
    
    real	13m30,204s
    user	12m30,980s
    sys	0m52,764s
    (interactomeseq) 




```bash
ls "$WorkingDir"/output/*_definition.*
```

    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1ssDNA_definition.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1ssDNA_domain_definition.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNA_definition.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNA_domain_definition.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2ssDNA_definition.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2ssDNA_domain_definition.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNA_definition.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNA_domain_definition.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/NotSelected_definition.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/NotSelected_domain_definition.log
    (interactomeseq) 



# Domain Enrichment


```bash
time python pyinteraseq_enrichment_eukaryotic.py --outputcontrol "$WorkingDir"/output/invSINEB2ssDNA_definition.txt  --outputarget "$WorkingDir"/output/invSINEB2tRNA_definition.txt --outputfolder "$WorkingDir"/output/ --outputid invSINEB2tRNAvsinvSINEB2ssDNA 
```

    
    real	0m2,024s
    user	0m1,804s
    sys	0m0,472s
    (interactomeseq) 




```bash
time python pyinteraseq_enrichment_eukaryotic.py --outputcontrol "$WorkingDir"/output/Delta5Uchl1ssDNA_definition.txt  --outputarget "$WorkingDir"/output/Delta5Uchl1tRNA_definition.txt --outputfolder "$WorkingDir"/output/ --outputid Delta5Uchl1tRNAvsDelta5Uchl1ssDNA 
```

    
    real	0m2,418s
    user	0m2,268s
    sys	0m0,620s
    (interactomeseq) 




```bash
ls "$WorkingDir"/output/*_intervals.*
ls "$WorkingDir"/output/*_enrichment.log
```

    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNAvsDelta5Uchl1ssDNA_common_intervals.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNAvsDelta5Uchl1ssDNA_unique_intervals.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNAvsinvSINEB2ssDNA_common_intervals.txt
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNAvsinvSINEB2ssDNA_unique_intervals.txt
    (interactomeseq) /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/Delta5Uchl1tRNAvsDelta5Uchl1ssDNA_enrichment.log
    /home/spuccio/InteractomSeqCLIEukaryotic/InteractomeSeqCLIEukaryotic/output/invSINEB2tRNAvsinvSINEB2ssDNA_enrichment.log
    (interactomeseq) 




```bash
echo "Done"
```

    Done
    (interactomeseq) 




```bash

```
