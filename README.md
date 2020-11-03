--This file includes my capstone project's PDF file,Microsoft powerpoint presentation and R codes.--The datasets were taken from the previously published paper of Miranda et al (https://academic.oup.com/hmg/article/25/22/4962/2525924)


#  Identification of Novel Peptides from Human Brain Samples: A Pilot Study


###### Capstone project in Bioinformatics: Hood College                     	          
###### Supervisor: Dr. Miranda Darby, Assistant Professor, Department of Bioinformatics
###### Title: “Novel Peptides Identification from human brain samples”
###### Brief Description: Identified multiple novel peptides expressed in human brain samples. Using R/BIOCONDUCTOR Package, ‘PGA’, I analyzed RNA-Seq data expressed in human brain tissues, converted to peptide sequences encoded by novel exons and finally search for their expression in samples available in public proteomic databases such as PRIDE using cluster computing. I have shown expression of novel repetitive elements in human brain samples and using USCS genome browser determining their genomic location.


# Abstract

Identification of novel peptides is very important to reveal the cause of different diseases, and the
field of Proteogenomics plays an important role in this regard. Here, novel peptides are identified
by searching MS/MS spectra against customized protein sequence databases. These databases
contain both known and predicted novel protein sequences, as well as sequence variants that are
generated based on genomic and transcriptomic sequence information. Previously our lab has
performed a strand-specific RNA sequencing on transcripts from 59 human orbitofrontal cortex
samples, and discovered a large number of transcribed and exonized Repetitive Elements (RE).
However, it is not known if these novel RE-containing exons are translated into humans, which
might have neurological disease relevance. The overall goal of this capstone project was to
determine if any of these putative RE-exons are translated to produce proteins in human cells
using the proteogenomics approach. Here, we used PGA, an R/Bioconductor package, that
enables an automatic process for constructing customized proteomic databases based upon RNA-
Seq data, and subsequently search peptides using MS/MS data from publicly available proteomic
databases. In this study, we used PRIDE, the most widely used proteomic database which has a
very rich deposit of proteomics samples in MS/MS form. After searching against 88 human brain
cortex samples, we have discovered 33 different peptides in 26 samples., among them 16 are
isomers. Based on these data, we predict that transcripts of RE-exons are potentially translated in
the human orbitofrontal cortex. Our study suggests that the PGA based Proteogenomic approach
could be a useful tool to identify novel peptides in RNA seq data.
