# SCIFER

Scripts to analyze single cell LINE-1 mRNA expression at single loci. Obtained from Stow et al. 2022.

## Abstract

### Background
Endogenous expression of L1 mRNA is the first step in an L1-initiated mutagenesis event. However, the contribution of individual cell types to patterns of organ-specific L1 mRNA expression remains poorly understood, especially at single-locus resolution. We introduce a method to quantify expression of mobile elements at the single-locus resolution in scRNA-Seq datasets called Single Cell Implementation to Find Expressed Retrotransposons (SCIFER). SCIFER aligns scRNA-Seq reads uniquely to the genome and extracts alignments from single cells by cell-specific barcodes. In contrast to the alignment performed using default parameters, this alignment strategy increases accuracy of L1 locus identification by retaining only reads that are uniquely mapped to individual L1 loci. L1 loci expressed in single cells are unambiguously identified using a list of L1 loci manually validated to be expressed in bulk RNA-Seq datasets generated from the same cell line or organ.

### Results
Validation of SCIFER using MCF7 cells determined technical parameters needed for optimal detection of L1 expression in single cells. We show that unsupervised analysis of L1 expression in single cells exponentially inflates both the levels of L1 expression and the number of expressed L1 loci. Application of SCIFER to analysis of scRNA-Seq datasets generated from mouse and human testes identified that mouse Round Spermatids and human Spermatogonia, Spermatocytes, and Round Spermatids express the highest levels of L1 mRNA. Our analysis also determined that similar to mice, human testes from unrelated individuals share as much as 80% of expressed L1 loci. Additionally, SCIFER determined that individual mouse cells co-express different L1 sub-families and different families of transposable elements, experimentally validating their co-existence in the same cell.

### Conclusions
SCIFER detects mRNA expression of individual L1 loci in single cells. It is compatible with scRNA-Seq datasets prepared using traditional sequencing methods. Validated using a human cancer cell line, SCIFER analysis of mouse and human testes identified key cell types supporting L1 expression in these species. This will further our understanding of differences and similarities in endogenous L1 mRNA expression patterns in mice and humans.




## Citation
Stow, E.C., Baddoo, M., LaRosa, A.J. et al. SCIFER: approach for analysis of LINE-1 mRNA expression in single cells at a single locus resolution. Mobile DNA 13, 21 (2022). https://doi.org/10.1186/s13100-022-00276-0
