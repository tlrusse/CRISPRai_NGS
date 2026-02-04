# Notes for CRISPRai_NGS
---
## 20260204
- go over mageck count

## 20260203
- ran fastqc
- demuxed the 4 samples using bowtie2 and samtools
- settled on just using the forward reads with CIGAR `6S30M115S`
  - also compared using paired end reads, but we lose 15% of the reads
  - fastqc report indicates that sequencing from reverse primer has quality issues which makes alignment wrse
  - also considered shifted matches for barcode, but only increases number of reads by \< 1%
- started project from template 🤖
