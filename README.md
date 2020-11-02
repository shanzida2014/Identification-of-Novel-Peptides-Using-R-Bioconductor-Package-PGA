# Identification-of-Novel-Peptides-Using-R-Bioconductor-Package-PGA
Identified Novel Peptides From RNA seq data of Human Brain Samples by R/ Bioconductor Package PGA
Identification of Novel Peptides from Human Brain Samples: A Pilot Study

Abstract

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
