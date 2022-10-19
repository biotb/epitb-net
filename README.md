![Title](https://github.com/biotb/Epistasis_TB/blob/master/figures/title.jpg)
# An epistatic network describes *oppA* and *glgB* as relevant genes for *Mycobacterium tuberculosis*
## Abstract
*Mycobacterium tuberculosis* is an acid-fast bacteria that causes tuberculosis disease worldwide. The role of epistatic interactions among different sites of the *M. tuberculosis* genome under selection by selective pressure may be crucial for understanding the disease and the molecular basis of antibiotic resistance acquisition. Here, we analyzed polymorphic loci interactions applying a model-free method for epistasis detection, SpydrPick, on a pan-genome-wide alignment created from a set of 254 complete reference genomes available to date. By means of the analysis of an epistatic network created with the detected epistatic interactions, we found that *glgB* (α-1,4-glucan branching enzyme) and *oppA* (oligopeptide-binding protein) are putative targets of co-selection in *Mycobacterium tuberculosis* since they were associated in the network with Tuberculosis genes related to virulence, pathogenesis, transport systems modulators of the immune response, and antibiotic resistance. In addition, our work unveiled potential pharmacological applications for genotypic antibiotic resistance inherent to mutations of *glgB* and *oppA*, as they epistatically interact with *fprA* and *embC*, two genes recently included as antibiotic-resistant genes in the catalogue of the World Health Organization. Our findings showed that this approach allows the identification of relevant epistatic interactions that may lead to a better understanding of *Mycobacterium tuberculosis* and that may be applied to different bacterial populations.

## R scripts used in this project

- [box-plot.R](https://github.com/biotb/epitb-net/blob/master/scripts/box-plot.R): Create box plot of MI scores.
- [circular-plot.R](https://github.com/biotb/epitb-net/blob/master/scripts/circular-plot.R): Create circular plots with GO annotations.
- [david-annotation.R](https://github.com/biotb/epitb-net/blob/master/scripts/david-annotation.R): Gene annotation with DAVID GOs.
- [geneid-annotation.R](https://github.com/biotb/epitb-net/blob/master/scripts/geneid-annotation.R): Gene annotation with ENTREZ id using GFF file.
- [genename-annotation.R](https://github.com/biotb/epitb-net/blob/master/scripts/genename-annotation.R): Gene name annotation using partitions generated by AMAS.
- [genes-in-whocatalogue.R](https://github.com/biotb/epitb-net/blob/master/scripts/genes-in-whocatalogue.R): Look for genes in WHO catalogue and annotation with WHO drugs.
- [phandango.R](https://github.com/biotb/epitb-net/blob/master/scripts/phandango.R): Create table with annotations for phandango plot with allele distribution.
- [snp-sites-outliers.R](https://github.com/biotb/epitb-net/blob/master/scripts/snp-sites-outliers.R): Compare positions from outliers to SNP-sites.

## Citation
### Article
Posada-Reyes, A. B., Balderas-Martínez, Y. I., Ávila-Ríos, S., Vinuesa, P., & Fonseca-Coronado, S. (2022). An epistatic network describes oppA and glgB as relevant genes for Mycobacterium tuberculosis. *Frontiers in Molecular Biosciences*, 9.
- Bibtex
```
@ARTICLE{10.3389/fmolb.2022.856212,
AUTHOR={Posada-Reyes, Ali-Berenice and Balderas-Martínez, Yalbi I. and Ávila-Ríos, Santiago and Vinuesa, Pablo and Fonseca-Coronado, Salvador},   
TITLE={An Epistatic Network Describes oppA and glgB as Relevant Genes for Mycobacterium tuberculosis},      
JOURNAL={Frontiers in Molecular Biosciences},      
VOLUME={9},           
YEAR={2022},      
URL={https://www.frontiersin.org/articles/10.3389/fmolb.2022.856212},       
DOI={10.3389/fmolb.2022.856212},      
ISSN={2296-889X},   
ABSTRACT={Mycobacterium tuberculosis is an acid-fast bacterium that causes tuberculosis worldwide. The role of epistatic interactions among different loci of the M. tuberculosis genome under selective pressure may be crucial for understanding the disease and the molecular basis of antibiotic resistance acquisition. Here, we analyzed polymorphic loci interactions by applying a model-free method for epistasis detection, SpydrPick, on a pan–genome-wide alignment created from a set of 254 complete reference genomes. By means of the analysis of an epistatic network created with the detected epistatic interactions, we found that glgB (α-1,4-glucan branching enzyme) and oppA (oligopeptide-binding protein) are putative targets of co-selection in M. tuberculosis as they were associated in the network with M. tuberculosis genes related to virulence, pathogenesis, transport system modulators of the immune response, and antibiotic resistance. In addition, our work unveiled potential pharmacological applications for genotypic antibiotic resistance inherent to the mutations of glgB and oppA as they epistatically interact with fprA and embC, two genes recently included as antibiotic-resistant genes in the catalog of the World Health Organization. Our findings showed that this approach allows the identification of relevant epistatic interactions that may lead to a better understanding of M. tuberculosis by deciphering the complex interactions of molecules involved in its metabolism, virulence, and pathogenesis and that may be applied to different bacterial populations.}
}
```
### Git citation
Please cite:
- APA

Posada-Reyes, Ali-Berenice. (2022). An epistatic network describes *oppA* and *glgB* as relevant genes for *Mycobacterium tuberculosis*. Retrieved from https://github.com/biotb/epitb-net.

- Vancouver

Posada-Reyes, Ali-Berenice. An epistatic network describes *oppA* and *glgB* as relevant genes for *Mycobacterium tuberculosis* [Internet]. 2022. Available from: https://github.com/biotb/epitb-net.

- Harvard

Posada-Reyes, Ali-Berenice, 2022. An epistatic network describes *oppA* and *glgB* as relevant genes for *Mycobacterium tuberculosis*, Available at: https://github.com/biotb/epitb-net.

- Bibtex
```
@misc{posada2022epitb-net,
  author={Posada-Reyes, Ali-Berenice},
  title={An epistatic network describes \emph{oppA} and \emph{glgB} as relevant genes for \emph{Mycobacterium tuberculosis},
  year={2022},
  url={https://github.com/biotb/epitb-net},
}
```

