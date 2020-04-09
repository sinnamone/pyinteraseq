# Testing pipeline with Prokaryotic dataset


```bash
https://www.dropbox.com/s/kx79noqq3kypbw6/Eukaryote_RnaBindProt.zip?dl=0
```


```bash
mkdir -p PhageRnaBinding
```


```bash
unzip -d PhageRnaBinding/ Eukaryote_RnaBindProt.zip
```


```bash
cd PhageRnaBinding/ 
```


```bash
gunzip *.gz
```


```bash
cd /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/
```


```bash
mkdir -p output
```


```bash
thread=30
```

# Mapping 

python pyinteraseq_main_mapping.py --readforward  /home/spuccio/PhageRnaBinding/Delta5Uchl1tRNA.fastq --outputfolder  /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --primer5forward GCAGCAAGCGGCGCGCATGC --primer3forward GCGCTTCGTCAT  --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa  --thread $thread --outputid Delta5Uchl1tRNA  --organism  Homo_sapiens --log /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/Delta5Uchl1tRNA.log

python pyinteraseq_main_mapping.py --readforward  /home/spuccio/PhageRnaBinding/Delta5Uchl1ssDNA.fastq --outputfolder  /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --primer5forward GCAGCAAGCGGCGCGCATGC --primer3forward GCGCTTCGTCAT  --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa  --thread $thread --outputid Delta5Uchl1ssDNA  --organism  Homo_sapiens --log /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/Delta5Uchl1ssDNA.log

python pyinteraseq_main_mapping.py --readforward  /home/spuccio/PhageRnaBinding/invSINEB2ssDNA.fastq --outputfolder  /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --primer5forward GCAGCAAGCGGCGCGCATGC --primer3forward GCGCTTCGTCAT  --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa  --thread 20 --outputid invSINEB2ssDNA --organism  Homo_sapiens --log /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2ssDNA.log

python pyinteraseq_main_mapping.py --readforward  /home/spuccio/PhageRnaBinding/invSINEB2tRNA.fastq --outputfolder  /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --primer5forward GCAGCAAGCGGCGCGCATGC --primer3forward GCGCTTCGTCAT  --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa  --thread 20 --outputid invSINEB2tRNA --organism  Homo_sapiens --log /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2tRNA.log

python pyinteraseq_main_mapping.py --readforward  /home/spuccio/PhageRnaBinding/NotSelected.fastq --outputfolder  /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --primer5forward GCAGCAAGCGGCGCGCATGC --primer3forward GCGCTTCGTCAT  --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa  --thread 20 --outputid NotSelected --organism  Homo_sapiens --log /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/NotSelected.log

# Domain Definition


```bash
python pyinteraseq_main_domain_definition.py --mappingoutput /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2tRNA_mapping.bam --outputfolder /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --outputid invSINEB2tRNA --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa   --annotation /home/spuccio/PhageRnaBinding/Human_annotation.bed --threshold 30 --genome /home/spuccio/PhageRnaBinding/sizes.transciptome 

```


```bash
python pyinteraseq_main_domain_definition.py --mappingoutput /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2ssDNA_mapping.bam --outputfolder /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --outputid invSINEB2ssDNA --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa   --annotation /home/spuccio/PhageRnaBinding/Human_annotation.bed --threshold 30 --genome /home/spuccio/PhageRnaBinding/sizes.transciptome 
```


```bash
python pyinteraseq_main_domain_definition.py --mappingoutput /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/Delta5Uchl1ssDNA_mapping.bam --outputfolder /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --outputid Delta5Uchl1ssDNA --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa   --annotation /home/spuccio/PhageRnaBinding/Human_annotation.bed --threshold 30 --genome /home/spuccio/PhageRnaBinding/sizes.transciptome 
```


```bash
python pyinteraseq_main_domain_definition.py --mappingoutput /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/Delta5Uchl1tRNA_mapping.bam --outputfolder /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --outputid Delta5Uchl1tRNA --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa   --annotation /home/spuccio/PhageRnaBinding/Human_annotation.bed --threshold 30 --genome /home/spuccio/PhageRnaBinding/sizes.transciptome 
```


```bash
python pyinteraseq_main_domain_definition.py --mappingoutput /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/NotSelected_mapping.bam --outputfolder /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --outputid NotSelected --fastasequence /home/spuccio/PhageRnaBinding/Homo_sapiens.GRCh38.cdna.all.fa   --annotation /home/spuccio/PhageRnaBinding/Human_annotation.bed --threshold 30 --genome /home/spuccio/PhageRnaBinding/sizes.transciptome 
```

# Domain Enrichment

python pyinteraseq_enrichment_eukaryotic.py --outputcontrol /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2ssDNA_definition.txt  --outputarget /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2tRNA_definition.txt --outputfolder /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --outputid invSINEB2tRNA

python pyinteraseq_enrichment_eukaryotic.py --outputcontrol /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2ssDNA_definition.txt  --outputarget /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/invSINEB2tRNA_definition.txt --outputfolder /home/spuccio/InteractomeSeq/PyinteraseqEukaryotic/output/ --outputid invSINEB2tRNA
