## Install MetaCurator and Dependencies

MetaCurator (https://github.com/RTRichar/MetaCurator)  
Python 2 or Python 3  
Perl 5.16.0 or compatible  
MAFFT 7.270 or compatible  
VSEARCH v2.8.1 or compatible  
HMMER v3.1 or compatible  

## Curate Data and Train Metaxa2
### Download sequence data from NCBI Nucleotide 
- rbcL, trnL, ITS2 and chloroplast genome sequences downloaded on 07/05/2018
- trnH sequences downloaded on 01/23/2019
- Short COI downloaded on 04/19/2019
- Arthropod whole mitochondrion genomes downloaded on 04/10/2019
#### NCBI Nucleotide search terms used:
*COI sequences:* '(cytochrome c oxidase subunit 1[All Fields] OR COI[All Fields]) AND (Arthropoda[All Fields] OR "Arthropoda"[Organism]) AND ("200"[SLEN] : "3000"[SLEN])'

*mitochondrial genome sequences:* '("Arthropoda"[Organism] OR Arthropoda[All Fields]) AND mitochondrion[All Fields] '

*trnH sequences:* 'trnH[All Fields] AND (plants[filter] AND ("100"[SLEN] : "3000"[SLEN]))'

*trnL sequences:* 'trnL[All Fields] AND plants[filter] AND ("0"[SLEN] : "10000"[SLEN])'

*rbcL sequences:* 'rbcL[All Fields] OR rubisco[All Fields] AND (plants[filter] AND ("0"[SLEN] : "10000"[SLEN]))'

*contiguous 5.8S-ITS2-28S ribosomal sequences:* '5.8S[All Fields] AND 28S[All Fields] AND plants[filter]'

*assorted ribosomal sequences:* 'ITS2[All Fields] OR 5.8S[All Fields] OR 28S[All Fields] AND (plants[filter] AND ("0"[SLEN] : "10000"[SLEN]))'

*chloroplast genome sequences:* 'genome[All Fields] AND (plants[filter] AND chloroplast[filter])'

#### Obtain NCBI accession for each entry and reformat sequences (We just show the trnL example, but this should be applied to all markers)
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < trnL.fasta | sed '/^$/d' > trnL_rmNL.fasta
cat rbcL_rmNL.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' > trnL_rmNL_AccH.fasta
grep '>' trnL.fasta | cut -d'|' -f4 > trnL_acc.csv
```
* Since NCBI changes how the data are formated over time, and depending on whether the data is downloaded through Edirect or from the site manually, these commands won't work for all instances. The desired format of the output is shown below so that users can adjust the re-formating to account for variation in input formatting. 
* Note the awk one-liner used to remove next line symbols within individual sequence entries. MetaCurator will fail if sequences are chunked into multiple lines, as is the case with sequences downloaded from NCBI.
##### After these steps you should have two important files, trnL_rmNL_AccH.fasta and trnL_acc.csv. The files should look like the examples below. No extra lines, no text after the accessions in the fasta header or the accessions text file
Contents of **trnL_rmNL_AccH.fasta**  
>\>AH015566.2  
>ATCTTAGCTATTAACTAGTTCGAAATTTTAAGTTCTACTTATAACTTATACTTATAAAAAAAAATACTAAAACTTCTTACAGATAAAGTTAGCTTGATATGCTTAACTAGAAGATATCTTTAAAAAACATTATATAATTTATTGAACTTTCTTTTTATTTTATTTCTCTAATTCGCAAATGCATTTTTCTATCTTTCTATCATAGAATAGAATTGATTCCAATTTCTATAATGGAACTGGATTTCAAATATTTTCAATTTGATATGGCTCGG  
>\>AH015565.2  
>AATAGTGTAACAAATAGAAATAGGTATAGTATAGGAAATCCGTAAAATCTCAGATCTTAGTTATTAATCTTAGCTATTAACTAGTTCGAAATTTTAAGTTCTACTTATAACTTATACTTATAAATTATAAAAAAAATACTAAAACTTCTTACAGATAAAGTTAGCTTGATATGCTTAACTAGAAGATATCTTTAAAAAACATTATATAATTTATTGAACTTTCTTTTTATTTTATTTCTCTAATTCGCAAATGCATTTTTCTATCATAGAAT  
>\>AH015564.2  
>CGTAAAATCTCAGATCTTAGGTTATTAATCTTAGCTATTAACTAGTTCGAAATTTTAAGTTCTACTTATAACTTATACTTATAAAAAAAAAATACTAAAACTTCTTACAGATAAAGTTAGCTTGATATGCTTAACTAGAAGATATCTTTAAAAAACATTATATAATTTATTGAACTTTCTTTTTATTTTATTTCTTTATTTCTCTAATTCGCAAATGCATTTTTCTATCATAGAATAGAATTGATTCCAATTTCTATAATGGAACTGGATTT

Contents of **trnL_acc.csv**  
>AH015566.2  
>AH015565.2  
>AH015564.2  

#### Upload accessions file to R and use Taxonomizr module to obtain lineage information for each accession with the following commands (In order to query the taxonomies of each accession, you must first download a local sql NCBI taxonomy database: https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html)
```
library(taxonomizr)
Accession <- as.vector(trnL_acc$V1)
TaxIDs <- accessionToTaxa(Accession,'Path/To/accessionTaxa.sql')
Lineages <- getTaxonomy(TaxIDs,'Path/To/accessionTaxa.sql')
Final.df <- cbind(Accession, TaxIDs, Lineages)
write.csv(Final.df, file = "trnL_Lineages.csv")
```
##### After these steps you should have two important files, trnL_rmNL_AccH.fasta and trnL_Lineages.csv. The taxonomy file should look like the example below. 
Contents of **trnL_Lineages.csv**  
>"  49992","Z37472.1","49992","Eukaryota","Streptophyta",NA,"Lamiales","Lamiaceae","Thymus","Thymus vulgaris"  
>"  49992","Z37471.1","49992","Eukaryota","Streptophyta",NA,"Lamiales","Lamiaceae","Thymus","Thymus vulgaris"  
>"  49991","Z37470.1","49991","Eukaryota","Streptophyta",NA,"Lamiales","Lamiaceae","Thymus","Thymus alsinoides"  
## Run data through MetaCurator pipeline 
### The short, automated way
```
# commands used for short NCBI Nucleotide entries
MetaCurator.py -i trnL_rmNL_AccH.fa -r trnLReps.fa -it trnL_Lineages.csv -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 5 -e 0.005 -ct True -tf True -of trnL_Amplicons.fa -ot trnL_Amplicons.tax
MetaCurator.py -i trnH_rmNL_AccH.fa -r trnH_Reps.fa -it trnH_Lineages.csv -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 5 -e 0.00005 -ct True -tf True -of trnH_Amplicons.fa -ot trnH_Amplicons.tax
MetaCurator.py -i rbcL_rmNL_AccH.fa -r rbcLReps.fa -it rbcL_Lineages.csv -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 5 -e 0.005 -ct True -tf True -of rbcL_Amplicons.fa -ot rbcL_Amplicons.tax
MetaCurator.py -i ITS2_rmNL_AccH.fa -r ITS2Reps.fa -it ITS2_lineages.csv -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 5 -e 0.005 -ct True -tf True -of ITS2_Amplicons.fa -ot ITS2_Amplicons.tax
MetaCurator.py -i COI.fa -r COI_Reps.fa -it COI_lineages.csv -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 5 -e 0.005 -ct True -tf True -of COI_Amplicons.fa -ot COI_Amplicons.tax
# commands used for whole chloroplast genomes
MetaCurator.py -i WCP_rmNL_AccH.fa -r rbcLReps.fa -it WCP_lineages.tax -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 10 -e 0.005 -of rbcL_WCP_Amplicons.fa -ot rbcL_WCP_Amplicons.tax
MetaCurator.py -i WCP_rmNL_AccH.fa -r trnLReps.fa -it WCP_lineages.tax -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 10 -e 0.005 -of trnL_WCP_Amplicons.fa -ot trnL_WCP_Amplicons.tax
MetaCurator.py -i WCP_rmNL_AccH.fa -r trnHReps.fa -it WCP_lineages.tax -is 8,3 -cs 1.0,0.98 -t 10 -e 0.00005 -of trnH_WCP_Amplicons.fa -ot trnH_WCP_Amplicons.tax
MetaCurator.py -i WMG_rmNL_AccH.fa -r COI_Reps.fa -it WMG_lineages.csv -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85 -t 5 -e 0.005 -ct True -tf True -of COI_WMG_Amplicons.fa -ot COI_WMG_Amplicons.tax
```
### An example of the long, manual way for trnL
#### Format tax lineages, clean up artifacts, extract marker using HMM-based IterRazor tool and dereplicate
```
Rtaxa2Mtaxa.py -i trnL_Lineages.csv -o trnL.tax
CleantTax.sh trnL.tax
ReviseIntNAs.py trnL_clean.tax
IterRazor.py -i trnL_rmNL_AccH.fasta -r trnLReps.fasta -o trnL_Amplicons.fasta -t 2 -st True -is 6,6,3,3,3 -cs 1.0,0.98,0.95,0.9,0.85
TaxFastaConsensus.py -it trnL_clean.tax -if trnL_Amplicons_rmNNN.fasta -ot trnL_trim.tax -of trnL_trim.fasta
DerepByTaxonomy.py -i trnL_trim.fasta -t trnL_trim.tax -o trnL_Amplicons_NR.fasta -st True
TaxFastaConsensus.py -it trnL_clean.tax -if trnL_Amplicons_NR.fasta -ot trnL_Final.tax -of trnL_Final.fasta
```
