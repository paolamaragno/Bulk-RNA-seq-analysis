# Bulk-RNA-seq analysis

The aim of the analysis is finding the differentially expressed genes across three different tissues - brain, liver and lung.
The data are downloaded from GTEx portal and, for each tissue, they consist of lot of samples coming from different individuals. Only three samples per tissue are selected depending on three quality measures - the RIN, the percentage of reads mapping on rRNA and the percentage of paired end reads both uniquely mapping on the genome.

The three samples are used as replicates of the same tissue: the overall nine replicates are analysed to identify the DE genes between each pair of tissues comparison.
Subsequently the lists of genes up regulated in one tissue with respect to both the other two are extracted and a functional enrichment analysis is performed to check if the enriched term are consistent with the tissue in which we know these genes are over expressed.

The same analysis is repeated filtering out from the nine replicates all the reads mapping on rRNAs, mitochondrial genes, psuedogenes, non canonical chromosomes or on genes of unknown type. 

The work ends with the comparison between the whole analysis starting from non-filtered and filtered data.
