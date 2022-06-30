Bulk RNA-seq analysis on brain, liver and lung replicates
================
Paola Maragno, Alberto Pettenella
17-06-2022

``` r
library(recount3)
library(edgeR)
```

    ## Warning: package 'limma' was built under R version 4.1.3

``` r
library(writexl)
```

# Differential expression analysis without filtering pseudogenes, rRNA genes, mitochondrial genes and genes of unknown type

## Data loading and pre-processing

**Load the data from GTEx portal**

``` r
rse_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)

rse_liver <-recount3::create_rse_manual(
  project = "LIVER",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)

rse_lung <- recount3::create_rse_manual(
  project = "LUNG",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)
```

In the assays of each Large Ranged Summarized Experiment object there is
the raw_counts table that contains the overall read coverage over the
genes’ exons. **To transform it into raw counts with the read counts for
each gene:**

``` r
assays(rse_brain)$counts <- recount3::transform_counts(rse_brain)
assays(rse_liver)$counts <- recount3::transform_counts(rse_liver)
assays(rse_lung)$counts <- recount3::transform_counts(rse_lung)
```

**Add another table to the assays containing the counts converted into
TPM values:**

``` r
assays(rse_brain)$TPM <- recount::getTPM(rse_brain,length_var = 'bp_length')
assays(rse_liver)$TPM <- recount::getTPM(rse_liver,length_var = 'bp_length')
assays(rse_lung)$TPM <- recount::getTPM(rse_lung,length_var = 'bp_length')
```

## Quality control for replicates selection

-   RIN -\> 6 or higher is considered acceptable

``` r
colData(rse_brain)$gtex.smrin[74]
```

    ## [1] 6.8

``` r
colData(rse_brain)$gtex.smrin[75]
```

    ## [1] 6.9

``` r
colData(rse_brain)$gtex.smrin[76]
```

    ## [1] 6.8

``` r
colData(rse_liver)$gtex.smrin[74]
```

    ## [1] 7.2

``` r
colData(rse_liver)$gtex.smrin[75]
```

    ## [1] 6.6

``` r
colData(rse_liver)$gtex.smrin[76]
```

    ## [1] 5.9

``` r
# too low
colData(rse_liver)$gtex.smrin[77]
```

    ## [1] 5.8

``` r
# too low
colData(rse_liver)$gtex.smrin[78]
```

    ## [1] 9.8

``` r
colData(rse_lung)$gtex.smrin[74]
```

    ## [1] 8.1

``` r
colData(rse_lung)$gtex.smrin[75]
```

    ## [1] 8

``` r
colData(rse_lung)$gtex.smrin[76]
```

    ## [1] 6.9

-   Estimated fraction of rRNA -\> never higher than 10%

``` r
colData(rse_brain)$gtex.smrrnart[74]
```

    ## [1] 0.0685402

``` r
colData(rse_brain)$gtex.smrrnart[75]
```

    ## [1] 0.0646821

``` r
colData(rse_brain)$gtex.smrrnart[76]
```

    ## [1] 0.0505885

``` r
colData(rse_liver)$gtex.smrrnart[74]
```

    ## [1] 0.00958842

``` r
colData(rse_liver)$gtex.smrrnart[75]
```

    ## [1] 0.018535

``` r
colData(rse_liver)$gtex.smrrnart[78]
```

    ## [1] 0.00737068

``` r
colData(rse_lung)$gtex.smrrnart[74]
```

    ## [1] 0.00452413

``` r
colData(rse_lung)$gtex.smrrnart[75]
```

    ## [1] 0.00296733

``` r
colData(rse_lung)$gtex.smrrnart[76]
```

    ## [1] 0.00345299

-   The percentage of mapped reads -\> at least 85% of reads uniquely
    mapped

``` r
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[74]
```

    ## [1] 89.1

``` r
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[75]
```

    ## [1] 91.8

``` r
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[76]
```

    ## [1] 91.4

``` r
colData(rse_liver)$"recount_qc.star.uniquely_mapped_reads_%_both"[74]
```

    ## [1] 86

``` r
colData(rse_liver)$"recount_qc.star.uniquely_mapped_reads_%_both"[75]
```

    ## [1] 88.3

``` r
colData(rse_liver)$"recount_qc.star.uniquely_mapped_reads_%_both"[78]
```

    ## [1] 88.6

``` r
colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"[74]
```

    ## [1] 90.3

``` r
colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"[75]
```

    ## [1] 92.1

``` r
colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"[76]
```

    ## [1] 89.7

**Build a rse object containing only the selected three replicates for
each tissue**

``` r
rse_brain_selected <- rse_brain[,c(74,75,76)]
rse_liver_selected <- rse_liver[,c(74,75,78)]
rse_lung_selected <- rse_lung[,c(74,75,76)]
```

**Fish out the count tables since edgeR needs only them to perform DE
analysis**

``` r
counts_brain_selected <- assays(rse_brain_selected)$counts
counts_liver_selected <- assays(rse_liver_selected)$counts
counts_lung_selected <- assays(rse_lung_selected)$counts
```

**Create a single count table containing the count tables of the three
replicates of each tissue,** assign to each gene its official gene
symbol that will be useful for the final enrichment analysis and rename
the columns

``` r
final_count_table <- cbind(counts_brain_selected, counts_lung_selected, counts_liver_selected)

rownames(final_count_table) <- rowData(rse_brain_selected)$gene_name

colnames(final_count_table) <- c("Brain74", "Brain75", "Brain76", "Lung74", "Lung75", "Lung76", "Liver74", "Liver75", "Liver78")
```

Check the library sizes

``` r
size <- colSums(final_count_table)
size
```

    ##  Brain74  Brain75  Brain76   Lung74   Lung75   Lung76  Liver74  Liver75 
    ## 27223813 31066981 29044162 31977213 29950710 29657908 33175303 30141545 
    ##  Liver78 
    ## 35166457

These are the number of reads of each replicate that will be used for
quantification. We can see that the sizes are not so much variable.

## Differential expression analysis

**Creation of “DGEList” object** that by the end will contain all the
information about the replicates as well as the parameters estimated
during the different steps of the DE analysis. This is the object on
which edgeR works.

``` r
y <- DGEList(counts=final_count_table)
y
```

    ## An object of class "DGEList"
    ## $counts
    ##              Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75
    ## SNX18P15           0       0       0      0      0      0       0       0
    ## SNX18P16           0       0       0      0      0      0       0       0
    ## ANKRD20A12P        0       0       0      0      0      0       0       0
    ## ANKRD20A15P        0       0       0      0      0      0       0       0
    ## LOC105379272       0       0       0      0      0      0       0       0
    ##              Liver78
    ## SNX18P15           0
    ## SNX18P16           0
    ## ANKRD20A12P        0
    ## ANKRD20A15P        0
    ## LOC105379272       0
    ## 54037 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors
    ## Brain74     1 27223813            1
    ## Brain75     1 31066981            1
    ## Brain76     1 29044162            1
    ## Lung74      1 31977213            1
    ## Lung75      1 29950710            1
    ## Lung76      1 29657908            1
    ## Liver74     1 33175303            1
    ## Liver75     1 30141545            1
    ## Liver78     1 35166457            1

Assign the group label to each replicate, that is the tissue from which
it derives

``` r
group <- as.factor(c("Brain", "Brain", "Brain", "Lung", "Lung", "Lung", "Liver", "Liver", "Liver"))
y$samples$group <- group
y
```

    ## An object of class "DGEList"
    ## $counts
    ##              Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75
    ## SNX18P15           0       0       0      0      0      0       0       0
    ## SNX18P16           0       0       0      0      0      0       0       0
    ## ANKRD20A12P        0       0       0      0      0      0       0       0
    ## ANKRD20A15P        0       0       0      0      0      0       0       0
    ## LOC105379272       0       0       0      0      0      0       0       0
    ##              Liver78
    ## SNX18P15           0
    ## SNX18P16           0
    ## ANKRD20A12P        0
    ## ANKRD20A15P        0
    ## LOC105379272       0
    ## 54037 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors
    ## Brain74 Brain 27223813            1
    ## Brain75 Brain 31066981            1
    ## Brain76 Brain 29044162            1
    ## Lung74   Lung 31977213            1
    ## Lung75   Lung 29950710            1
    ## Lung76   Lung 29657908            1
    ## Liver74 Liver 33175303            1
    ## Liver75 Liver 30141545            1
    ## Liver78 Liver 35166457            1

Add other additional information: sex, age, part of the tissue from
which the sample was taken and also some quality information we used to
select them

``` r
y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex, colData(rse_lung_selected)$gtex.sex, colData(rse_liver_selected)$gtex.sex))

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age, colData(rse_lung_selected)$gtex.age, colData(rse_liver_selected)$gtex.age))

y$samples$part_tissue <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd, colData(rse_lung_selected)$gtex.smtsd, colData(rse_liver_selected)$gtex.smtsd))

y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_lung_selected)$gtex.smrin, colData(rse_liver_selected)$gtex.smrin))

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_lung_selected)$gtex.smrrnart, colData(rse_liver_selected)$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_lung_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(rse_liver_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_lung_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_liver_selected)$"recount_qc.aligned_reads%.chrm"))

y
```

    ## An object of class "DGEList"
    ## $counts
    ##              Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75
    ## SNX18P15           0       0       0      0      0      0       0       0
    ## SNX18P16           0       0       0      0      0      0       0       0
    ## ANKRD20A12P        0       0       0      0      0      0       0       0
    ## ANKRD20A15P        0       0       0      0      0      0       0       0
    ## LOC105379272       0       0       0      0      0      0       0       0
    ##              Liver78
    ## SNX18P15           0
    ## SNX18P16           0
    ## ANKRD20A12P        0
    ## ANKRD20A15P        0
    ## LOC105379272       0
    ## 54037 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors sex   age                     part_tissue
    ## Brain74 Brain 27223813            1   1 60-69 Brain - Putamen (basal ganglia)
    ## Brain75 Brain 31066981            1   1 60-69        Brain - Substantia nigra
    ## Brain76 Brain 29044162            1   1 60-69 Brain - Putamen (basal ganglia)
    ## Lung74   Lung 31977213            1   1 20-29                            Lung
    ## Lung75   Lung 29950710            1   1 40-49                            Lung
    ## Lung76   Lung 29657908            1   1 50-59                            Lung
    ## Liver74 Liver 33175303            1   1 60-69                           Liver
    ## Liver75 Liver 30141545            1   1 60-69                           Liver
    ## Liver78 Liver 35166457            1   2 40-49                           Liver
    ##         rin       rRNA mapped  chrm
    ## Brain74 6.8  0.0685402   89.1 27.46
    ## Brain75 6.9  0.0646821   91.8 22.85
    ## Brain76 6.8  0.0505885   91.4 25.41
    ## Lung74  8.1 0.00452413   90.3  6.21
    ## Lung75    8 0.00296733   92.1  6.36
    ## Lung76  6.9 0.00345299   89.7  4.12
    ## Liver74 7.2 0.00958842     86  8.58
    ## Liver75 6.6   0.018535   88.3 19.06
    ## Liver78 9.8 0.00737068   88.6  8.35

Remove all the genes with low or zero counts (expression): first check
how many genes do not appear in any of the 9 replicates

``` r
table(rowSums(y$counts==0)==9)
```

    ## 
    ## FALSE  TRUE 
    ## 39381 14661

**keep.exprs function** keeps all the genes that are expressed in all
the three replicates of a given tissue by removing those with zero or
low expression since they have to be ignored during normalization and
parameter estimation

``` r
keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

dim(y)
```

    ## [1] 23782     9

**23782 genes are those genes that remain after the filtering and are
those on which the DE analysis will be based on.**

Extract and store in a vector the log2 of the counts per million before
normalization and plot their distribution

``` r
logcpm_before <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("Brain74","Brain75","Brain76") , 'mediumseagreen' , 
                   ifelse(colnames(logcpm_before) %in% c("Lung74","Lung75","Lung76"), 'lightskyblue',
                          'sienna1' ) )
boxplot(logcpm_before,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) before TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Tissues",
       c("Brain","Lung","Liver"), fill=c('mediumseagreen','lightskyblue','sienna1'), horiz=TRUE, cex=0.35)
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

The medians are quite different

``` r
for (i in 1:9){
  print(median(logcpm_before[,i]))
}
```

    ## [1] 3.176728
    ## [1] 3.121314
    ## [1] 3.223566
    ## [1] 3.167707
    ## [1] 3.652474
    ## [1] 3.31125
    ## [1] 2.226591
    ## [1] 1.57835
    ## [1] 2.096291

**edgeR normalizes the counts using TMM normalization**

``` r
y <- calcNormFactors(y, method = "TMM")
y
```

    ## An object of class "DGEList"
    ## $counts
    ##           Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75 Liver78
    ## MIR6859-1       4      11       7      8     27     21      10       3       6
    ## WASH7P        996    1076     985   1001   1834   1623     638     386     503
    ## SEPT14P18      12      14       5     23     32     30      21      14      18
    ## CICP27          9       3       3     11     19     13       4       3       3
    ## LOC729737     186     217     234   2791   3832   4040     263     248     408
    ## 23777 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors sex   age                     part_tissue
    ## Brain74 Brain 27178816    1.1252452   1 60-69 Brain - Putamen (basal ganglia)
    ## Brain75 Brain 31029877    1.1167931   1 60-69        Brain - Substantia nigra
    ## Brain76 Brain 28999673    1.1464474   1 60-69 Brain - Putamen (basal ganglia)
    ## Lung74   Lung 31946342    1.2466685   1 20-29                            Lung
    ## Lung75   Lung 29896245    1.5815826   1 40-49                            Lung
    ## Lung76   Lung 29610717    1.3509198   1 50-59                            Lung
    ## Liver74 Liver 33150560    0.6999862   1 60-69                           Liver
    ## Liver75 Liver 30123390    0.5124543   1 60-69                           Liver
    ## Liver78 Liver 35145007    0.7264553   2 40-49                           Liver
    ##         rin       rRNA mapped  chrm
    ## Brain74 6.8  0.0685402   89.1 27.46
    ## Brain75 6.9  0.0646821   91.8 22.85
    ## Brain76 6.8  0.0505885   91.4 25.41
    ## Lung74  8.1 0.00452413   90.3  6.21
    ## Lung75    8 0.00296733   92.1  6.36
    ## Lung76  6.9 0.00345299   89.7  4.12
    ## Liver74 7.2 0.00958842     86  8.58
    ## Liver75 6.6   0.018535   88.3 19.06
    ## Liver78 9.8 0.00737068   88.6  8.35

**TMM normalization is based on multiplying each count value of each
sample for a constant, that is the normalization factor defined at
sample level, trying to shift the count values gene by gene in order to
have no change of expression across each pair of samples. The assumption
of TMM normalization is that most of the genes do not change their
expression in a significant way.**

Extract and store in a vector the log2 of the counts per million after
normalization and plot their distribution

``` r
logcpm_after <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_after) %in% c("Brain74","Brain75","Brain76") , 'mediumseagreen' , 
                   ifelse(colnames(logcpm_after) %in% c("Lung74","Lung75","Lung76"), 'lightskyblue',
                          'sienna1' ) )
boxplot(logcpm_after,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) after TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Tissues",
       c("Brain","Lung","Liver"), fill=c('mediumseagreen','lightskyblue','sienna1'), horiz=TRUE, cex=0.35)
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

**The effect of TMM normalization is that now the medians are very
similar**

``` r
for (i in 1:9){
  print(median(logcpm_after[,i]))
}
```

    ## [1] 3.007276
    ## [1] 2.962683
    ## [1] 3.027361
    ## [1] 2.851631
    ## [1] 2.994919
    ## [1] 2.88006
    ## [1] 2.734556
    ## [1] 2.52676
    ## [1] 2.550637

**Design the linear model** without the intercept since it makes no
sense to choose a type of tissue as reference being the three samples
independent from one another

``` r
design <- model.matrix(~0 + group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design
```

    ##         Brain Liver Lung
    ## Brain74     1     0    0
    ## Brain75     1     0    0
    ## Brain76     1     0    0
    ## Lung74      0     0    1
    ## Lung75      0     0    1
    ## Lung76      0     0    1
    ## Liver74     0     1    0
    ## Liver75     0     1    0
    ## Liver78     0     1    0
    ## attr(,"assign")
    ## [1] 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$group
    ## [1] "contr.treatment"

### Exploratory analysis

**MDS plotting the samples labeled by group**

``` r
plotMDS(logcpm_after, labels=group, main = 'Multidimensional scaling plot of distances 
        between gene expression profiles - replicate label')
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

**After normalization the different samples are projected in
bidimensional space.** The distance between points is the overall fold
ratio gene by gene between two samples and it is computed only on the
top 500 variable genes. The closer two points are the more similar are
the expression values of these two samples: we can see that the three
replicates of the same tissue are close between each other and far from
the replicates of other tissues.

In case of brain and liver replicates one sample is little farther from
the other two: by plotting the samples labeling the points with
different quality information we can try to understand which may be the
most relevant sources of variability between replicates

``` r
plotMDS(logcpm_after, labels=y$samples$rRNA, main = 'Multidimensional scaling plot of distances 
        between gene expression profiles - rRNA% label')
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

Liver, lung and brain replicates that are farther from the others of the
same tissue don’t have a very different % of reads mapping on rRNAs, so
it isn’t probably the main source of variability.

``` r
plotMDS(logcpm_after, labels=y$samples$chrm, main = 'Multidimensional scaling plot of distances 
        between gene expression profiles - chrm% label')
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Looking at the three liver replicates the behavior is quite strange
since one replicate has a content of mitochondrial genes that is doubled
than the other two but the distance between the replicates with \~8.4
value and the one with 19 value is not so higher than the distance
between the two replicated with \~8.4 value themselves. Consequently it
is difficult to say that the percentage of mitochondrial genes explains
the variability of the three liver replicates.

One of the replicate of brain tissue has a percentage that is lower than
the one of the other two replicates, that instead have more similar
values, so probably in the case of brain the % of reads mapping on
mitochondrial RNAs can explain way a replicate is little further from
the other two.

``` r
plotMDS(logcpm_after, labels=y$samples$age, main = 'Multidimensional scaling plot of distances 
        between gene expression profiles - age label')
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Plotting the MDS using other labels of each sample - age, RIN, sex or
the part of the tissue from which the sample was extracted - doesn’t
help in finding the causes of variability between replicates of the same
tissue.

**Maybe the variability between replicates is due to the combination of
more factors that can’t be simply plotted.**

**Estimate the Negative Binomial dispersion and plot the BCV (square
root of the dispersion)**

``` r
y <- estimateDisp(y, design)
plotBCV(y, main = 'Biological Coefficient of Variation with respect to the average logCPM', ylim = c(0.15,1.75))
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

**The bigger is the BCV the higher is the dispersion, and so the
variance, of the corresponding gene.** The **Common red line** is the
global NB dispersion φ: it is estimated from all the genes. But, one
single dispersion value does not fit well all the genes and, on the
other hand, the assumption is that we have too few replicates to have a
reliable estimate of the dispersion of each gene. The **Trend blue
line** is the estimated trend that tries to model the dependence between
the mean expression and the dispersion.

The **estimated trend is used to shrinkage the dispersion of each gene
towards the trend line itself so that the final gene-wise dispersion
estimate is no more the observed dispersion of that gene, but the
original dispersion value of that gene modified by pulling it towards
the estimated trend.**

As we can see the common BCV is little below 0.5: it is quite high since
we are in presence of the maximum variability possible - sex, age,
different part of the same tissues, sample preparation.

We can extract all the parameters of NB distribution applied to the
data, as well as the estimations of the common and trended dispersion,
that have been stored in y

``` r
y$common.dispersion
```

    ## [1] 0.2208488

``` r
head(y$trended.dispersion)
```

    ## [1] 0.3305180 0.1513322 0.3401495 0.3278763 0.1468383 0.1533044

``` r
head(y$tagwise.dispersion)
```

    ## [1] 0.20963232 0.07970198 0.13963519 0.17935516 0.07504895 0.09527074

### Fit our data to the “generalized linear” model we designed

``` r
fit <- glmQLFit(y, design)
```

fit object contains all the information obtained from DE analysis; now
we only have to extract the information we want.

### Design pairwise comparison liver vs brain

``` r
qlfLB <- glmQLFTest(fit, contrast = c(-1,1,0))
```

In this object there are all the information about the specified
comparison and also the table with all the results of the DE analysis
for each gene.

``` r
head(qlfLB$table)
```

    ##                    logFC     logCPM            F     PValue
    ## MIR6859-1     0.39051536 -1.3923331  0.387188898 0.54720918
    ## WASH7P       -0.36914209  4.9098727  1.823078017 0.20565942
    ## SEPT14P18     1.39718881 -0.5946152  7.372486279 0.02107863
    ## CICP27        0.03036465 -1.8110435  0.002172412 0.96370958
    ## LOC729737     1.16132071  5.1325235 19.059353895 0.00129073
    ## LOC102723897 -0.44792698  4.8244088  2.087062563 0.17806075

This is the table in which gene by gene the log2FC, the log2CPM of the
same gene between the two samples and the p-value (not adjusted!!!) are
reported.

``` r
summary(decideTests(qlfLB, p.value=0.01, lfc=1))
```

    ##        -1*Brain 1*Liver
    ## Down               4507
    ## NotSig            15628
    ## Up                 3647

Setting as thresholds a p-value of 0.01 and a value of log2FC of 1: 3647
genes result up regulated in liver than in brain; 4507 genes are up
regulated in brain with respect to liver; 15628 genes don’t change
significantly their expression.

``` r
tableLB <- topTags(qlfLB, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
```

**topTags** extract the table contained in qlfLB object, performs the
p-value adjustment for multiple testing using the method that you
specify, in this case the Benjamini-Hochberg, and eventually sorts the
genes by adjusted p-value (FDR). On the top there are the more
significant DE genes (lowest p-value adjusted).

Now tableLB object contains all the genes which expression has been
compared between the two tissues. We now convert tableLB in dataFrame
type to perform some processing

``` r
tableLB <- as.data.frame(tableLB) 
dim(tableLB)
```

    ## [1] 23782     5

Filter this table to **keep only the significantly DE genes between the
two tissues**

``` r
tableLB <- tableLB[which((tableLB$logFC > 1 | tableLB$logFC < -1) & tableLB$FDR <0.01),]
dim(tableLB)
```

    ## [1] 8154    5

Genes will be considered “DE” if their FDR \< 0.01 and their log2FC is
\< -1 (in this case the gene is considered up regulated in brain) or \>
1 (in this case the gene is considered up regulated in liver).

Add to each significantly DE gene present in tableLB object the
indication of whether it is up regulated in brain or in liver

``` r
tableLB <- cbind(tableLB, upLIVER = "", upBRAIN = "")

for (i in 1:nrow(tableLB)) {
  if (tableLB[i,]$logFC > 1) 
  {tableLB[i,]$upLIVER <- rownames(tableLB)[i]}
  else 
  {tableLB[i,]$upBRAIN <- rownames(tableLB)[i]}
}

head(tableLB)
```

    ##               logFC   logCPM        F       PValue          FDR upLIVER
    ## HNF4A     12.469956 6.695189 660.9921 9.782045e-11 1.456350e-06   HNF4A
    ## SLC24A2  -11.507029 6.012513 603.1495 1.562126e-10 1.456350e-06        
    ## ARG1      13.373487 8.691813 556.7582 2.350411e-10 1.456350e-06    ARG1
    ## RUNDC3A  -10.355847 6.503334 552.2692 2.449499e-10 1.456350e-06        
    ## MAPK8IP2  -7.851144 5.958834 455.8388 6.507596e-10 2.318712e-06        
    ## ADH1C     11.452567 7.531051 449.4883 6.988473e-10 2.318712e-06   ADH1C
    ##           upBRAIN
    ## HNF4A            
    ## SLC24A2   SLC24A2
    ## ARG1             
    ## RUNDC3A   RUNDC3A
    ## MAPK8IP2 MAPK8IP2
    ## ADH1C

Remove the genes with low expression, log2CPM \< 0, since they are more
likely false positives

``` r
tableLB <- tableLB[-(which(tableLB$logCPM<0)),] 
```

Remove all the genes with a name starting with:

-   LOC: since they are those for which the offical gene symbol is not
    avaiable
-   LINC: Long Intergenic Non-Protein Coding
-   MIR: MicroRNA
-   SNORD: Small nucleolar RNA
-   RPL: corresponding to ribosomal proteins

``` r
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'LOC'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'LINC'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'MIR'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'SNORD'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'RPL'))),]
```

Save the final table with the interesting sorted DE genes between liver
and brain in xlsx format

``` r
tableLB <- cbind (GeneName = rownames(tableLB), tableLB)
write_xlsx(tableLB, "resultsLB_nf.xlsx")
```

### Design pairwise comparison lung vs brain

``` r
qlfLungB <- glmQLFTest(fit, contrast = c(-1,0,1))
```

``` r
head(qlfLungB$table)
```

    ##                    logFC     logCPM            F       PValue
    ## MIR6859-1     0.96131673 -1.3923331 2.622321e+00 1.353292e-01
    ## WASH7P        0.16437937  4.9098727 3.638902e-01 5.593071e-01
    ## SEPT14P18     1.07379013 -0.5946152 4.522561e+00 5.837805e-02
    ## CICP27        1.09053709 -1.8110435 3.485202e+00 9.039971e-02
    ## LOC729737     3.70044284  5.1325235 1.627778e+02 1.128231e-07
    ## LOC102723897 -0.01509469  4.8244088 2.390318e-03 9.619345e-01

``` r
summary(decideTests(qlfLungB, p.value=0.01, lfc=1))
```

    ##        -1*Brain 1*Lung
    ## Down              3739
    ## NotSig           16376
    ## Up                3667

3667 genes result up regulated in lung than in brain; 3739 genes are up
regulated in brain with respect to lung; 16376 genes don’t change
significantly their expression.

``` r
tableLungB <- topTags(qlfLungB, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
tableLungB <- as.data.frame(tableLungB) 
dim(tableLungB)
```

    ## [1] 23782     5

``` r
tableLungB <- tableLungB[which((tableLungB$logFC > 1 | tableLungB$logFC < -1) & tableLungB$FDR <0.01),]
dim(tableLungB)
```

    ## [1] 7406    5

``` r
tableLungB <- cbind(tableLungB, upLUNG = "", upBRAIN = "")

for (i in 1:nrow(tableLungB)) {
  if (tableLungB[i,]$logFC > 1) 
  {tableLungB[i,]$upLUNG <- rownames(tableLungB)[i]}
  else 
  {tableLungB[i,]$upBRAIN <- rownames(tableLungB)[i]}
}

head(tableLungB)
```

    ##               logFC   logCPM        F       PValue          FDR upLUNG  upBRAIN
    ## SLC24A2  -11.076122 6.012513 631.5034 1.235307e-10 2.937806e-06         SLC24A2
    ## RUNDC3A   -8.730518 6.503334 478.4584 5.087116e-10 5.508050e-06         RUNDC3A
    ## TTYH1    -10.515483 7.064360 415.6798 1.039402e-09 5.508050e-06           TTYH1
    ## CCL21      9.712783 4.847504 395.6221 1.335575e-09 5.508050e-06  CCL21         
    ## PIGR       9.186875 5.889645 394.9431 1.347249e-09 5.508050e-06   PIGR         
    ## MAPK8IP2  -6.839764 5.958834 391.2169 1.413526e-09 5.508050e-06        MAPK8IP2

``` r
tableLungB <- tableLungB[-(which(tableLungB$logCPM<0)),] 
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'LOC'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'LINC'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'MIR'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'SNORD'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'RPL'))),]
```

``` r
tableLungB <- cbind (GeneName = rownames(tableLungB), tableLungB)
write_xlsx(tableLungB, "resultsLungB_nf.xlsx")
```

### Design pairwise comparison lung vs liver

``` r
qlfLL <- glmQLFTest(fit, contrast = c(0,-1,1))
```

``` r
head(qlfLL$table)
```

    ##                   logFC     logCPM          F       PValue
    ## MIR6859-1     0.5708014 -1.3923331  0.8990465 3.645957e-01
    ## WASH7P        0.5335215  4.9098727  3.8010498 7.874549e-02
    ## SEPT14P18    -0.3233987 -0.5946152  0.4313801 5.256101e-01
    ## CICP27        1.0601724 -1.8110435  2.9856304 1.135985e-01
    ## LOC729737     2.5391221  5.1325235 84.8071492 2.572737e-06
    ## LOC102723897  0.4328323  4.8244088  1.9505617 1.916784e-01

``` r
summary(decideTests(qlfLL, p.value=0.01, lfc=1))
```

    ##        -1*Liver 1*Lung
    ## Down              2905
    ## NotSig           17399
    ## Up                3478

3478 genes result up regulated in lung than in liver; 2905 genes are up
regulated in liver with respect to lung: 17399 genes don’t change
significantly their expression.

``` r
tableLL <- topTags(qlfLL, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
tableLL <- as.data.frame(tableLL) 
dim(tableLL)
```

    ## [1] 23782     5

``` r
tableLL <- tableLL[which((tableLL$logFC > 1 | tableLL$logFC < -1) & tableLL$FDR <0.01),]
dim(tableLL)
```

    ## [1] 6383    5

``` r
tableLL <- cbind(tableLL, upLUNG = "", upLIVER = "")

for (i in 1:nrow(tableLL)) {
  if (tableLL[i,]$logFC > 1) 
  {tableLL[i,]$upLUNG <- rownames(tableLL)[i]}
  else 
  {tableLL[i,]$upLIVER <- rownames(tableLL)[i]}
}

head(tableLL)
```

    ##            logFC   logCPM        F       PValue          FDR upLUNG upLIVER
    ## TTPA   -9.111105 5.338826 622.1597 1.333097e-10 3.170371e-06           TTPA
    ## HNF4A  -9.457624 6.695189 524.1504 3.197291e-10 3.801899e-06          HNF4A
    ## ACMSD  -9.033755 5.557636 445.6783 7.297337e-10 3.870113e-06          ACMSD
    ## PDZK1  -7.260148 5.272410 406.4062 1.165395e-09 3.870113e-06          PDZK1
    ## ABCG5  -8.233808 5.592499 400.1744 1.260347e-09 3.870113e-06          ABCG5
    ## APOF  -15.137481 7.133531 396.6887 1.317479e-09 3.870113e-06           APOF

``` r
tableLL <- tableLL[-(which(tableLL$logCPM<0)),] 
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'LOC'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'LINC'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'MIR'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'SNORD'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'RPL'))),]
```

``` r
tableLL <- cbind (GeneName = rownames(tableLL), tableLL)
write_xlsx(tableLL, "resultsLL_nf.xlsx")
```

We want to find the genes that are up regulated in brain both in
comparison with liver and lung: they will be considered the more
interesting for the subsequent enrichment analysis. Also, find the genes
that are up-regulated in liver with respect to both brain and lung and
those up-regulated in lung with respect to both brain and liver.

**Up regulated genes in brain both in comparison with liver and lung**

``` r
genes <- rownames(tableLB[which(tableLB$upBRAIN != ""),])
both_brain <- vector()
for (i in genes) {
  if (i %in% tableLungB$upBRAIN)
  {both_brain <- c(both_brain, i)}
}
```

Save the just obtained list of genes in xlsx format

``` r
write.table(both_brain, "both_brain_nf.txt")
```

**Up regulated genes in liver both in comparison with brain and lung**

``` r
genes <- rownames(tableLB[which(tableLB$upLIVER != ""),])
both_liver <- vector()
for (i in genes) {
  if (i %in% tableLL$upLIVER)
  {both_liver <- c(both_liver, i)}
}
```

``` r
write.table(both_liver, "both_liver_nf.txt")
```

**Up regulated in lung both in comparison with brain and liver**

``` r
genes <- rownames(tableLL[which(tableLL$upLUNG != ""),])
both_lung <- vector()
for (i in genes) {
  if (i %in% tableLungB$upLUNG)
  {both_lung <- c(both_lung, i)}
}
```

``` r
write.table(both_lung, "both_lung_nf.txt")
```

Looking at the three tables that contain the most significantly DE genes
for each comparison, **SLC24A2 gene is significantly up regulated in
brain both with respect to lung and to liver.** We plot the distribution
of TPM values of this gene in the three tissues, without considering the
division in replicates, to see if there is actually an evident
difference

``` r
which(rowData(rse_brain)$gene_name == "SLC24A2")
```

    ## [1] 49995

``` r
myColors <- ifelse(assays(rse_brain)$TPM[49995,] , 'mediumseagreen' , 
                   ifelse(assays(rse_liver)$TPM[49995,], 'lightskyblue',
                          'sienna1' ) )
boxplot(assays(rse_brain)$TPM[49995,], assays(rse_liver)$TPM[49995,], assays(rse_lung)$TPM[49995,], 
        xlab='Tissue Samples', ylab='TPM values', main='TPM value of SLC24A2 gene 
        considering all the replicates of each tissue',
        names=c('Brain','Liver','Lung'),outline=F, notch=F, varwidth=T , col=myColors)
legend("topright", inset=.02, title="Tissues",
       c("Brain","Lung","Liver"), fill=c('mediumseagreen','lightskyblue','sienna1'), horiz=F, cex=0.35)
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

We can see that indeed this gene is highly expressed in brain, while it
is not expressed in liver and lung. From this plot we can say that
SLC24A2 can be considered “statistically” over expressed in brain than
in the other two tissues.

Computation of the median expression of SLC24A2 in the three tissues

``` r
median(assays(rse_brain)$TPM[49995,])
```

    ## [1] 16.69491

``` r
median(assays(rse_liver)$TPM[49995,])
```

    ## [1] 0.005448531

``` r
median(assays(rse_lung)$TPM[49995,])
```

    ## [1] 0.01738319

**Assess if the difference of expression of SLC24A2 gene is
statistically significant** between brain and the other two tissues by
considering all the samples of each tissue

``` r
tpm_brain <- assays(rse_brain)$TPM[49995,]
tpm_liver <- assays(rse_liver)$TPM[49995,]
tpm_lung <- assays(rse_lung)$TPM[49995,]
```

Since the TPM values of SLC24A2 gene in brain samples are independent
from the values in the samples of the other two tissues we use the
Mann–Whitney U test performing a one-sided test to check if:

-   H0: means of TPM values of the different tissues are the same
-   H1: mean of TPM values in brain is bigger than in the other tissue

**Mann–Whitney U test for TPM values of SLC24A2 gene in brain vs liver**

``` r
wilcox.test(tpm_brain, tpm_liver, paired=F, alternative='greater')
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  tpm_brain and tpm_liver
    ## W = 735676, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is greater than 0

Since the p-value is very low, \<2.2e-16, we can conclude that the mean
TPM value of SLC24A2 gene in brain is significantly higher than in
liver.

**Mann–Whitney U test for TPM values of SLC24A2 gene in brain vs lung**

``` r
wilcox.test(tpm_brain, tpm_lung, paired=F, alternative='greater')
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  tpm_brain and tpm_lung
    ## W = 1919795, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is greater than 0

Since the p-value is very low, \<2.2e-16, we can conclude that the mean
TPM value of SLC24A2 gene in brain is significantly higher than in lung.

## Functional enrichment analysis

``` r
library('enrichR')
setEnrichrSite("Enrichr")
websiteLive <- TRUE
```

### Up regulated genes in brain with respect to both lung and liver

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(both_brain, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-73-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-73-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(both_brain, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-75-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-75-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-75-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(both_brain, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

### Up regulated genes in liver with respect to both lung and brain

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(both_liver, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-79-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-79-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(both_liver, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-81-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-81-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-81-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(both_liver, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-83-1.png)<!-- -->

### Up regulated genes in lung with respect to both liver and brain

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(both_lung, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-85-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-85-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-85-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(both_lung, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-87-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-87-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-87-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(both_lung, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-89-1.png)<!-- -->

# Differential expression analysis after filtering pseudogenes, rRNA genes, mitochondrial genes and genes of unknown type

# Data loading and pre-processing

**Load the data from GTEx portal**

``` r
rm(list=ls())

rse_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)

rse_liver <-recount3::create_rse_manual(
  project = "LIVER",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)

rse_lung <- recount3::create_rse_manual(
  project = "LUNG",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)
```

In the assays of each Large Ranged Summarized Experiment object there is
the raw_counts table that contains the overall read coverage over the
genes’ exons. **To transform it into raw counts with the read counts for
each gene:**

``` r
assays(rse_brain)$counts <- recount3::transform_counts(rse_brain)
assays(rse_liver)$counts <- recount3::transform_counts(rse_liver)
assays(rse_lung)$counts <- recount3::transform_counts(rse_lung)
```

**Add another table to the assays containing the counts converted into
TPM values:**

``` r
assays(rse_brain)$TPM <- recount::getTPM(rse_brain,length_var = 'bp_length')
assays(rse_liver)$TPM <- recount::getTPM(rse_liver,length_var = 'bp_length')
assays(rse_lung)$TPM <- recount::getTPM(rse_lung,length_var = 'bp_length')
```

**Save the number of counts before filtering:**

``` r
dim_before_brain<- dim(assays(rse_brain)$counts)[1]
dim_before_liver<- dim(assays(rse_liver)$counts)[1]
dim_before_lung<- dim(assays(rse_lung)$counts)[1]
```

**Clean the RSE objects by removing:**

-   rRNA genes
-   Pseudogenes
-   mitochondrial genes
-   genes on non-canonical chromosomes

Firstly, define a canonical-chromosomes list:

``` r
canonical <- paste("chr", seq(1,22), sep="")
canonical <- c(canonical, "chrX", "chrY")
```

Now perform the filtering:

**Brain**

``` r
rse_brain <- rse_brain[
    # Ribosomal RNA
    rowData(rse_brain)$gbkey != 'rRNA' &
    # Pseudogenes
    rowData(rse_brain)$gbkey != 'Gene' &
    # Exclude Non-canonical Chromosomes and Mitochondrial DNA
    rowRanges(rse_brain)@seqnames %in% canonical &
    # NAs
    !is.na(rowData(rse_brain)$gbkey),
]
```

**Liver**

``` r
rse_liver <- rse_liver[
    # Ribosomal RNA
    rowData(rse_liver)$gbkey != 'rRNA' &
    # Pseudogenes
    rowData(rse_liver)$gbkey != 'Gene' &
    # Exclude Non-canonical Chromosomes and Mitochondrial DNA
    rowRanges(rse_liver)@seqnames %in% canonical &
    # NAs
    !is.na(rowData(rse_liver)$gbkey),
]
```

**Lung**

``` r
rse_lung <- rse_lung[
    # Ribosomal RNA
    rowData(rse_lung)$gbkey != 'rRNA' &
    # Pseudogenes
    rowData(rse_lung)$gbkey != 'Gene' &
    # Exclude Non-canonical Chromosomes and Mitochondrial DNA
    rowRanges(rse_lung)@seqnames %in% canonical &
    # NAs
    !is.na(rowData(rse_lung)$gbkey),
]
```

**Retrieve the number of filtered genes for each tissue**

``` r
dim_before_brain-dim(assays(rse_brain)$counts)[1]
```

    ## [1] 18293

``` r
dim_before_liver-dim(assays(rse_liver)$counts)[1]
```

    ## [1] 18293

``` r
dim_before_lung-dim(assays(rse_lung)$counts)[1]
```

    ## [1] 18293

## Quality control for replicates selection

As rRNA genes have been filtered out of the RSE object, the RNA
integrity number and the percentage of mapped reads are the only
criteria considered for choosing the replicates.

-   RIN -\> 6 or higher is considered acceptable

``` r
samples <- c(74,75,76)
```

Brain samples:

``` r
for( i in samples){
  print(colData(rse_brain)$gtex.smrin[i])
}
```

    ## [1] 6.8
    ## [1] 6.9
    ## [1] 6.8

Liver samples:

``` r
for( i in samples){
  print(colData(rse_liver)$gtex.smrin[i])
}
```

    ## [1] 7.2
    ## [1] 6.6
    ## [1] 5.9

``` r
colData(rse_liver)$gtex.smrin[77]
```

    ## [1] 5.8

As for samples 76 and 77 the RIN is too low, sample 78 will be
considered.

``` r
colData(rse_liver)$gtex.smrin[78]
```

    ## [1] 9.8

``` r
samples_liv <- c(74,75,78)
```

Lung samples:

``` r
for( i in samples){
  print(colData(rse_lung)$gtex.smrin[i])
}
```

    ## [1] 8.1
    ## [1] 8
    ## [1] 6.9

-   The percentage of mapped reads -\> at least 85% of reads uniquely
    mapped

Brain samples:

``` r
for (i in samples){
  print(colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[i])
} 
```

    ## [1] 89.1
    ## [1] 91.8
    ## [1] 91.4

Liver samples:

``` r
for (i in samples_liv){
  print(colData(rse_liver)$"recount_qc.star.uniquely_mapped_reads_%_both"[i])
} 
```

    ## [1] 86
    ## [1] 88.3
    ## [1] 88.6

Lung samples:

``` r
for (i in samples){
  print(colData(rse_lung)$"recount_qc.star.uniquely_mapped_reads_%_both"[i])
} 
```

    ## [1] 90.3
    ## [1] 92.1
    ## [1] 89.7

**Build a rse object containing only the selected three replicates for
each tissue**

``` r
rse_brain_selected <- rse_brain[,samples]
rse_liver_selected <- rse_liver[,samples_liv]
rse_lung_selected <- rse_lung[,samples]
```

**Fish out the count tables since edgeR needs only them to perform DE
analysis**

``` r
counts_brain_selected <- assays(rse_brain_selected)$counts
counts_liver_selected <- assays(rse_liver_selected)$counts
counts_lung_selected <- assays(rse_lung_selected)$counts
```

**Create the single count table containing the count tables of the three
replicates of each tissue,** assign to each gene its official gene
symbol that will be useful for the final enrichment analysis and rename
the columns

``` r
final_count_table <- cbind(counts_brain_selected, counts_lung_selected, counts_liver_selected)

rownames(final_count_table) <- rowData(rse_brain_selected)$gene_name

colnames(final_count_table) <- c("Brain74", "Brain75", "Brain76", "Lung74", "Lung75", "Lung76", "Liver74", "Liver75", "Liver78")
```

Check the library sizes

``` r
size <- colSums(final_count_table)
size
```

    ##  Brain74  Brain75  Brain76   Lung74   Lung75   Lung76  Liver74  Liver75 
    ## 17023625 19380213 19200063 23055297 20695392 20587134 27305273 25023092 
    ##  Liver78 
    ## 28494356

These are the number of reads of each replicate that will be used for
quantification. We can see that the sizes are not so much variable
relative to the same tissue considered.

## Differential expression analysis

**Creation of “DGEList” object** that by the end will contain all the
information about the replicates as well as the parameters estimated
during the different steps of the DE analysis. This is the object on
which edgeR works.

``` r
y <- DGEList(counts=final_count_table)
y
```

    ## An object of class "DGEList"
    ## $counts
    ##           Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75 Liver78
    ## DDX11L1         6      10       6      8     25     17       4       4       9
    ## MIR6859-1       2       6       3      4     14     11       5       2       3
    ## MIR6859-1       4      11       7      8     27     21      10       3       6
    ## WASH7P        996    1076     985   1001   1834   1623     638     386     503
    ## MIR1302-2       0       0       0      0      0      0       0       0       0
    ## 35744 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors
    ## Brain74     1 17023625            1
    ## Brain75     1 19380213            1
    ## Brain76     1 19200063            1
    ## Lung74      1 23055297            1
    ## Lung75      1 20695392            1
    ## Lung76      1 20587134            1
    ## Liver74     1 27305273            1
    ## Liver75     1 25023092            1
    ## Liver78     1 28494356            1

Assign the group label to each replicate, that is the tissue from which
it derives

``` r
group <- as.factor(c("Brain", "Brain", "Brain", "Lung", "Lung", "Lung", "Liver", "Liver", "Liver"))
y$samples$group <- group
y
```

    ## An object of class "DGEList"
    ## $counts
    ##           Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75 Liver78
    ## DDX11L1         6      10       6      8     25     17       4       4       9
    ## MIR6859-1       2       6       3      4     14     11       5       2       3
    ## MIR6859-1       4      11       7      8     27     21      10       3       6
    ## WASH7P        996    1076     985   1001   1834   1623     638     386     503
    ## MIR1302-2       0       0       0      0      0      0       0       0       0
    ## 35744 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors
    ## Brain74 Brain 17023625            1
    ## Brain75 Brain 19380213            1
    ## Brain76 Brain 19200063            1
    ## Lung74   Lung 23055297            1
    ## Lung75   Lung 20695392            1
    ## Lung76   Lung 20587134            1
    ## Liver74 Liver 27305273            1
    ## Liver75 Liver 25023092            1
    ## Liver78 Liver 28494356            1

Add other additional information: sex, age, part of the tissue from
which the sample was taken and also some quality information we used to
select them

``` r
y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex, colData(rse_lung_selected)$gtex.sex, colData(rse_liver_selected)$gtex.sex))

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age, colData(rse_lung_selected)$gtex.age, colData(rse_liver_selected)$gtex.age))

y$samples$part_tissue <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd, colData(rse_lung_selected)$gtex.smtsd, colData(rse_liver_selected)$gtex.smtsd))

y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_lung_selected)$gtex.smrin, colData(rse_liver_selected)$gtex.smrin))

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_lung_selected)$gtex.smrrnart, colData(rse_liver_selected)$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_lung_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(rse_liver_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_lung_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_liver_selected)$"recount_qc.aligned_reads%.chrm"))

y
```

    ## An object of class "DGEList"
    ## $counts
    ##           Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75 Liver78
    ## DDX11L1         6      10       6      8     25     17       4       4       9
    ## MIR6859-1       2       6       3      4     14     11       5       2       3
    ## MIR6859-1       4      11       7      8     27     21      10       3       6
    ## WASH7P        996    1076     985   1001   1834   1623     638     386     503
    ## MIR1302-2       0       0       0      0      0      0       0       0       0
    ## 35744 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors sex   age                     part_tissue
    ## Brain74 Brain 17023625            1   1 60-69 Brain - Putamen (basal ganglia)
    ## Brain75 Brain 19380213            1   1 60-69        Brain - Substantia nigra
    ## Brain76 Brain 19200063            1   1 60-69 Brain - Putamen (basal ganglia)
    ## Lung74   Lung 23055297            1   1 20-29                            Lung
    ## Lung75   Lung 20695392            1   1 40-49                            Lung
    ## Lung76   Lung 20587134            1   1 50-59                            Lung
    ## Liver74 Liver 27305273            1   1 60-69                           Liver
    ## Liver75 Liver 25023092            1   1 60-69                           Liver
    ## Liver78 Liver 28494356            1   2 40-49                           Liver
    ##         rin       rRNA mapped  chrm
    ## Brain74 6.8  0.0685402   89.1 27.46
    ## Brain75 6.9  0.0646821   91.8 22.85
    ## Brain76 6.8  0.0505885   91.4 25.41
    ## Lung74  8.1 0.00452413   90.3  6.21
    ## Lung75    8 0.00296733   92.1  6.36
    ## Lung76  6.9 0.00345299   89.7  4.12
    ## Liver74 7.2 0.00958842     86  8.58
    ## Liver75 6.6   0.018535   88.3 19.06
    ## Liver78 9.8 0.00737068   88.6  8.35

From the data we can see that:

-   there is no anatomical reference in Lung and Liver tissues
-   only sample Liver78 is a woman
-   only sample Lung74 is 20-29 y.o.

Remove all the genes with low or zero counts (expression): first check
how many genes do not appear in any of the 9 replicates

``` r
table(rowSums(y$counts==0)==9)
```

    ## 
    ## FALSE  TRUE 
    ## 28845  6904

**keep.exprs function** keeps all the genes that are expressed in all
the three replicates of a given tissue by removing those with zero or
low expression since they have to be ignored during normalization and
parameter estimation

``` r
keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

dim(y)
```

    ## [1] 16908     9

**16908 genes are those genes that remain after the filtering and are
those on which the DE analysis will be based on.**

Extract and store in a vector the log2 of the counts per million before
normalization and plot their distribution

``` r
logcpm_before <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("Brain74","Brain75","Brain76") , 'mediumseagreen' , 
                   ifelse(colnames(logcpm_before) %in% c("Lung74","Lung75","Lung76"), 'lightskyblue',
                          'sienna1' ) )
boxplot(logcpm_before,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) before TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Tissues",
       c("Brain","Lung","Liver"), fill=c('mediumseagreen','lightskyblue','sienna1'), horiz=TRUE, cex=0.35)
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-116-1.png)<!-- -->

The medians are visibly different.

``` r
for (i in 1:9){
  print(median(logcpm_before[,i]))
}
```

    ## [1] 3.510045
    ## [1] 3.42436
    ## [1] 3.472546
    ## [1] 3.249696
    ## [1] 3.893337
    ## [1] 3.471027
    ## [1] 2.130248
    ## [1] 1.405698
    ## [1] 1.917275

**edgeR normalizes the counts using TMM normalization**

``` r
y <- calcNormFactors(y, method = "TMM")
y
```

    ## An object of class "DGEList"
    ## $counts
    ##              Brain74 Brain75 Brain76 Lung74 Lung75 Lung76 Liver74 Liver75
    ## DDX11L1            6      10       6      8     25     17       4       4
    ## MIR6859-1          4      11       7      8     27     21      10       3
    ## WASH7P           996    1076     985   1001   1834   1623     638     386
    ## LOC729737        186     217     234   2791   3832   4040     263     248
    ## LOC102723897     769    1181    1137    799   1573   1565     534     411
    ##              Liver78
    ## DDX11L1            9
    ## MIR6859-1          6
    ## WASH7P           503
    ## LOC729737        408
    ## LOC102723897     482
    ## 16903 more rows ...
    ## 
    ## $samples
    ##         group lib.size norm.factors sex   age                     part_tissue
    ## Brain74 Brain 16986665    1.3725692   1 60-69 Brain - Putamen (basal ganglia)
    ## Brain75 Brain 19349914    1.3099488   1 60-69        Brain - Substantia nigra
    ## Brain76 Brain 19163226    1.3230761   1 60-69 Brain - Putamen (basal ganglia)
    ## Lung74   Lung 23031515    1.2431830   1 20-29                            Lung
    ## Lung75   Lung 20650017    1.7055322   1 40-49                            Lung
    ## Lung76   Lung 20547338    1.4164243   1 50-59                            Lung
    ## Liver74 Liver 27284492    0.5859955   1 60-69                           Liver
    ## Liver75 Liver 25003125    0.4044998   1 60-69                           Liver
    ## Liver78 Liver 28476038    0.5905077   2 40-49                           Liver
    ##         rin       rRNA mapped  chrm
    ## Brain74 6.8  0.0685402   89.1 27.46
    ## Brain75 6.9  0.0646821   91.8 22.85
    ## Brain76 6.8  0.0505885   91.4 25.41
    ## Lung74  8.1 0.00452413   90.3  6.21
    ## Lung75    8 0.00296733   92.1  6.36
    ## Lung76  6.9 0.00345299   89.7  4.12
    ## Liver74 7.2 0.00958842     86  8.58
    ## Liver75 6.6   0.018535   88.3 19.06
    ## Liver78 9.8 0.00737068   88.6  8.35

**TMM normalization is based on multiplying each count value of each
sample for a constant, that is the normalization factor defined at
sample level, trying to shift the count values gene by gene in order to
have no change of expression across each pair of samples. The assumption
of TMM normalization is that most of the genes do not change their
expression in a significant way.**

Extract and store in a vector the log2 of the counts per million after
normalization and plot their distribution

``` r
logcpm_after <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_after) %in% c("Brain74","Brain75","Brain76") , 'mediumseagreen' , 
                   ifelse(colnames(logcpm_after) %in% c("Lung74","Lung75","Lung76"), 'lightskyblue',
                          'sienna1' ) )
boxplot(logcpm_after,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) after TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Tissues",
       c("Brain","Lung","Liver"), fill=c('mediumseagreen','lightskyblue','sienna1'), horiz=TRUE, cex=0.35)
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-119-1.png)<!-- -->

**The effect of TMM normalization is that now the medians are very
similar**

``` r
for (i in 1:9){
  print(median(logcpm_after[,i]))
}
```

    ## [1] 3.056699
    ## [1] 3.037879
    ## [1] 3.071725
    ## [1] 2.938208
    ## [1] 3.128589
    ## [1] 2.972891
    ## [1] 2.888203
    ## [1] 2.681173
    ## [1] 2.66222

**Design the linear model** without the intercept since it makes no
sense to choose a type of tissue as reference being the three samples
independent from one another

``` r
design <- model.matrix(~0 + group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design
```

    ##         Brain Liver Lung
    ## Brain74     1     0    0
    ## Brain75     1     0    0
    ## Brain76     1     0    0
    ## Lung74      0     0    1
    ## Lung75      0     0    1
    ## Lung76      0     0    1
    ## Liver74     0     1    0
    ## Liver75     0     1    0
    ## Liver78     0     1    0
    ## attr(,"assign")
    ## [1] 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$group
    ## [1] "contr.treatment"

### Exploratory analysis

**MDS plotting the samples labeled by group**

``` r
plotMDS(logcpm_after, labels=group, main = 'Multidimensional scaling plot of distances between gene expression profiles - replicate label')
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-122-1.png)<!-- -->

**After normalization the samples are projected in bidimensional
space.** The distance between points is the overall fold ratio gene by
gene between two samples and it is computed only on a subset of genes
that are more variable. The closer two points are the more similar are
the expression values of these two samples: we can see that the three
replicates of the same tissue are close between each other and far from
the replicates of other tissues.

Although age-groups are not particularly different across samples, we
can still test whether age represents a source of variability.

``` r
plotMDS(logcpm_after, labels=y$samples$age, main = 'Multidimensional scaling plot of distances between gene expression profiles - age label')
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-123-1.png)<!-- -->

Age does not seem to be an important source of variability.

Sex and part of the tissue don’t help explaining the variability between
replicates of the same tissue as there are not enough categories for
each of them.

**Estimate the Negative Binomial dispersion and plot the BCV (square
root of the dispersion)**

``` r
y <- estimateDisp(y, design)
plotBCV(y, main = 'Biological Coefficient of Variation with respect to the average logCPM',ylim=c(0.15,1.75))
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-124-1.png)<!-- -->

**The bigger is the BCV the higher is the dispersion, and so the
variance, of the corresponding gene.** The **Common red line** is the
global NB dispersion φ: it is estimated from all the genes. But, one
single dispersion value does not fit well all the genes and, on the
other hand, the assumption is that we have too few replicates to have a
reliable estimate of the dispersion of each gene. The **Trend blue
line** is the estimated trend that tries to model the dependence between
the mean expression and the dispersion.

The **estimated trend is used to shrinkage the dispersion of each gene
towards the trend line itself so that the final gene-wise dispersion
estimate is no more the observed dispersion of that gene, but the
original dispersion value of that gene modified by pulling it towards
the estimated trend.**

As we can see the common BCV is around 0.5: it is quite high since we
are in presence of the maximum variability possible - sex, age,
different part of the same tissues, sample preparation.

We can extract all the parameters of NB distribution applied to the
data, as well as the estimations of the common and trended dispersion,
that have been stored in y

``` r
y$common.dispersion
```

    ## [1] 0.2428987

``` r
head(y$trended.dispersion)
```

    ## [1] 0.3467876 0.3483031 0.1675418 0.1620653 0.1698110 0.3611342

``` r
head(y$tagwise.dispersion)
```

    ## [1] 0.1802937 0.2121563 0.0849840 0.0837677 0.1028146 0.1681516

### Fit our data to the “generalized linear” model we designed

``` r
fit <- glmQLFit(y, design)
```

fit object contains all the information obtained from DE analysis; now
we only have to extract the information we want.

### Design pairwise comparison liver vs brain

``` r
qlfLB <- glmQLFTest(fit, contrast = c(-1,1,0))
```

In this object there are all the information about the comparison and
also the table with all the results of the DE analysis for each gene.

``` r
head(qlfLB$table)
```

    ##                   logFC     logCPM          F       PValue
    ## DDX11L1       0.4050471 -1.0155373  0.4659026 0.5096557248
    ## MIR6859-1     0.5445128 -0.9196589  0.7372563 0.4097134302
    ## WASH7P       -0.1999297  5.3832616  0.5127079 0.4895864945
    ## LOC729737     1.3388569  5.6022805 22.4168163 0.0006982783
    ## LOC102723897 -0.2773662  5.2987897  0.7548063 0.4043789765
    ## LOC100996442 -0.3540817  1.2724814  0.5367770 0.4798065875

This is the table in which gene by gene the log2FC, the log2CPM of the
same gene between the two samples and the p-value (not adjusted!!!) are
reported.

``` r
summary(decideTests(qlfLB, p.value=0.01, lfc=1))
```

    ##        -1*Brain 1*Liver
    ## Down               3163
    ## NotSig            10769
    ## Up                 2976

Setting as thresholds a p-value of 0.01 and a value of log2FC of 1: 2976
genes result up regulated in liver than in brain; 3163 genes are up
regulated in brain with respect to liver; 10769 genes don’t change
significantly their expression.

``` r
tableLB <- topTags(qlfLB, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
```

**topTags** extract the table contained in qlfLB object, performs the
p-value adjustment for multiple testing using the method that you
specify, in this case the Benjamini-Hochberg, and eventually sorts the
genes by adjusted p-value (FDR). On the top there are the more
significant DE genes (lowest p-value adjusted).

Now tableLB object contain all the genes which expression has been
compared between the two tissues. We now convert tableLB in dataFrame
type to perform some processing

``` r
tableLB <- as.data.frame(tableLB) 
dim(tableLB)
```

    ## [1] 16908     5

Filter this table to **keep only the significantly DE genes between the
two tissues**

``` r
tableLB <- tableLB[which((tableLB$logFC > 1 | tableLB$logFC < -1) & tableLB$FDR <0.01),]
dim(tableLB)
```

    ## [1] 6139    5

Genes will be considered “DE” if their FDR \< 0.01 and their log2FC is
\< -1 (in this case the gene is considered up regulated in brain) or \>
1 (in this case the gene is considered up regulated in liver).

Add to each significantly DE gene present in tableLB object the
indication of whether it is up regulated in brain or in liver

``` r
tableLB <- cbind(tableLB, upLIVER = "", upBRAIN = "")

for (i in 1:nrow(tableLB)) {
  if (tableLB[i,]$logFC > 1) 
  {tableLB[i,]$upLIVER <- rownames(tableLB)[i]}
  else 
  {tableLB[i,]$upBRAIN <- rownames(tableLB)[i]}
}

head(tableLB)
```

    ##               logFC   logCPM        F       PValue          FDR upLIVER
    ## HNF4A     12.635039 7.276518 621.4574 1.099778e-10 1.560782e-06   HNF4A
    ## SLC24A2  -11.345156 6.423933 562.1621 1.846205e-10 1.560782e-06        
    ## ARG1      13.544453 9.274236 508.1327 3.108843e-10 1.752144e-06    ARG1
    ## ABCG5      9.297106 6.170104 442.5481 6.330197e-10 2.253743e-06   ABCG5
    ## ACMSD      9.205816 6.137283 431.8174 7.180900e-10 2.253743e-06   ACMSD
    ## MAPK8IP2  -7.676356 6.370645 409.6013 9.417038e-10 2.253743e-06        
    ##           upBRAIN
    ## HNF4A            
    ## SLC24A2   SLC24A2
    ## ARG1             
    ## ABCG5            
    ## ACMSD            
    ## MAPK8IP2 MAPK8IP2

Remove the genes with low expression, log2CPM \< 0, since they are more
likely false positives

``` r
tableLB <- tableLB[-(which(tableLB$logCPM<0)),] 
```

Remove all the genes with a name starting with:

-   LOC: since they are those for which the offical gene symbol is not
    avaiable
-   LINC: Long Intergenic Non-Protein Coding
-   MIR: MicroRNA
-   SNORD: Small nucleolar RNA
-   RPL: corresponding to ribosomal proteins

``` r
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'LOC'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'LINC'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'MIR'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'SNORD'))),]
tableLB <- tableLB[-(which(startsWith(rownames(tableLB), 'RPL'))),]
```

Save the final table with the interesting sorted DE genes between liver
and brain in xlsx format

``` r
tableLB <- cbind (GeneName = rownames(tableLB), tableLB)
write_xlsx(tableLB, "resultsLB_f.xlsx")
```

### Design pairwise comparison lung vs brain

``` r
qlfLungB <- glmQLFTest(fit, contrast = c(-1,0,1))
```

``` r
head(qlfLungB$table)
```

    ##                   logFC     logCPM            F       PValue
    ## DDX11L1      0.81698334 -1.0155373 2.144865e+00 1.723347e-01
    ## MIR6859-1    0.99010936 -0.9196589 2.717034e+00 1.288434e-01
    ## WASH7P       0.20347467  5.3832616 5.332033e-01 4.812369e-01
    ## LOC729737    3.74243206  5.6022805 1.477062e+02 1.616512e-07
    ## LOC102723897 0.02210563  5.2987897 4.821792e-03 9.459447e-01
    ## LOC100996442 0.18559239  1.2724814 1.528470e-01 7.036448e-01

``` r
summary(decideTests(qlfLungB, p.value=0.01, lfc=1))
```

    ##        -1*Brain 1*Lung
    ## Down              2740
    ## NotSig           11355
    ## Up                2813

2813 genes result up regulated in lung than in brain; 2740 genes are up
regulated in brain with respect to lung; 11355 genes don’t change
significantly their expression.

``` r
tableLungB <- topTags(qlfLungB, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
tableLungB <- as.data.frame(tableLungB) 
dim(tableLungB)
```

    ## [1] 16908     5

``` r
tableLungB <- tableLungB[which((tableLungB$logFC > 1 | tableLungB$logFC < -1) & tableLungB$FDR <0.01),]
dim(tableLungB)
```

    ## [1] 5553    5

``` r
tableLungB <- cbind(tableLungB, upLUNG = "", upBRAIN = "")

for (i in 1:nrow(tableLungB)) {
  if (tableLungB[i,]$logFC > 1) 
  {tableLungB[i,]$upLUNG <- rownames(tableLungB)[i]}
  else 
  {tableLungB[i,]$upBRAIN <- rownames(tableLungB)[i]}
}

head(tableLungB)
```

    ##               logFC   logCPM        F       PValue          FDR upLUNG  upBRAIN
    ## SLC24A2  -11.037620 6.423933 596.4563 1.359788e-10 2.299130e-06         SLC24A2
    ## CCL21      9.755054 5.308559 373.1396 1.518376e-09 7.844177e-06  CCL21         
    ## PIGR       9.223353 6.366964 368.4003 1.621014e-09 7.844177e-06   PIGR         
    ## MAPK8IP2  -6.797967 6.370645 358.7878 1.855731e-09 7.844177e-06        MAPK8IP2
    ## SCNN1A     9.665294 6.596118 343.1345 2.330827e-09 7.881925e-06 SCNN1A         
    ## PTPRN     -8.137321 5.897438 329.3072 2.875294e-09 8.102579e-06           PTPRN

``` r
tableLungB <- tableLungB[-(which(tableLungB$logCPM<0)),] 
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'LOC'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'LINC'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'MIR'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'SNORD'))),]
tableLungB <- tableLungB[-(which(startsWith(rownames(tableLungB), 'RPL'))),]
```

``` r
tableLungB <- cbind (GeneName = rownames(tableLungB), tableLungB)
write_xlsx(tableLungB, "resultsLungB_f.xlsx")
```

### Design pairwise comparison lung vs liver

``` r
qlfLL <- glmQLFTest(fit, contrast = c(0,-1,1))
```

``` r
head(qlfLL$table)
```

    ##                  logFC     logCPM          F       PValue
    ## DDX11L1      0.4119362 -1.0155373  0.5189933 4.869989e-01
    ## MIR6859-1    0.4455965 -0.9196589  0.5351782 4.804456e-01
    ## WASH7P       0.4034043  5.3832616  2.0842107 1.779934e-01
    ## LOC729737    2.4035751  5.6022805 68.3740657 6.442031e-06
    ## LOC102723897 0.2994719  5.2987897  0.8802172 3.692206e-01
    ## LOC100996442 0.5396741  1.2724814  1.2514648 2.882113e-01

``` r
summary(decideTests(qlfLL, p.value=0.01, lfc=1))
```

    ##        -1*Liver 1*Lung
    ## Down              2263
    ## NotSig           12331
    ## Up                2314

2314 genes result up regulated in lung than in liver; 2263 genes are up
regulated in liver with respect to lung: 12331 genes don’t change
significantly their expression.

``` r
tableLL <- topTags(qlfLL, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
tableLL <- as.data.frame(tableLL) 
dim(tableLL)
```

    ## [1] 16908     5

``` r
tableLL <- tableLL[which((tableLL$logFC > 1 | tableLL$logFC < -1) & tableLL$FDR <0.01),]
dim(tableLL)
```

    ## [1] 4577    5

``` r
tableLL <- cbind(tableLL, upLUNG = "", upLIVER = "")

for (i in 1:nrow(tableLL)) {
  if (tableLL[i,]$logFC > 1) 
  {tableLL[i,]$upLUNG <- rownames(tableLL)[i]}
  else 
  {tableLL[i,]$upLIVER <- rownames(tableLL)[i]}
}

head(tableLL)
```

    ##              logFC   logCPM        F       PValue          FDR upLUNG upLIVER
    ## TTPA     -9.244489 5.920096 590.8181 1.428175e-10 2.414758e-06           TTPA
    ## HNF4A    -9.587160 7.276518 489.4336 3.771106e-10 3.188093e-06          HNF4A
    ## ACMSD    -9.164285 6.137283 434.6639 6.942631e-10 3.793416e-06          ACMSD
    ## ABCG5    -8.363290 6.170104 397.7186 1.095154e-09 3.793416e-06          ABCG5
    ## APOF    -15.257367 7.707248 393.4908 1.156814e-09 3.793416e-06           APOF
    ## SLCO1B1 -15.824156 6.695730 738.5750 1.459749e-09 3.793416e-06        SLCO1B1

``` r
tableLL <- tableLL[-(which(tableLL$logCPM<0)),] 
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'LOC'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'LINC'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'MIR'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'SNORD'))),]
tableLL <- tableLL[-(which(startsWith(rownames(tableLL), 'RPL'))),]
```

``` r
tableLL <- cbind (GeneName = rownames(tableLL), tableLL)
write_xlsx(tableLL, "resultsLL_f.xlsx")
```

We want to find the genes that are up regulated in brain both in
comparison with liver and lung: they will be considered the more
interesting for the subsequent enrichment analysis. Also, find the genes
that are up-regulated in liver with respect to both brain and lung and
those up-regulated in lung with respect to both brain and liver.

**Up regulated genes in brain both in comparison with liver and lung**

``` r
genes <- rownames(tableLB[which(tableLB$upBRAIN != ""),])
both_brain <- vector()
for (i in genes) {
  if (i %in% tableLungB$upBRAIN)
  {both_brain <- c(both_brain, i)}
}
```

Save the just obtained list of genes in txt format

``` r
write.table(both_brain, "both_brain_f.txt")
```

**Up regulated genes in liver both in comparison with brain and lung**

``` r
genes <- rownames(tableLB[which(tableLB$upLIVER != ""),])
both_liver <- vector()
for (i in genes) {
  if (i %in% tableLL$upLIVER)
  {both_liver <- c(both_liver, i)}
}
```

``` r
write.table(both_liver, "both_liver_f.txt")
```

**Up regulated in lung both in comparison with brain and liver**

``` r
genes <- rownames(tableLL[which(tableLL$upLUNG != ""),])
both_lung <- vector()
for (i in genes) {
  if (i %in% tableLungB$upLUNG)
  {both_lung <- c(both_lung, i)}
}
```

``` r
write.table(both_lung, "both_lung_f.txt")
```

Looking at the three tables that contain the most significantly DE genes
for each comparison, **SLC24A2 gene is significantly up regulated in
brain both with respect to lung and to liver.** Plot the distribution of
TPM values of this gene in the three tissues, without considering the
division in replicates, to see if there is actually an evident
difference

``` r
assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_liver)$TPM <- recount::getTPM(rse_liver)
assays(rse_lung)$TPM <- recount::getTPM(rse_lung)
#which(rowData(rse_brain)$gene_name == "SLC24A2")
myColors <- ifelse(assays(rse_brain)$TPM[32893,] , 'mediumseagreen' , 
                   ifelse(assays(rse_liver)$TPM[32893,], 'lightskyblue',
                          'sienna1' ) )
boxplot(assays(rse_brain)$TPM[32893,], assays(rse_liver)$TPM[32893,], assays(rse_lung)$TPM[32893,], 
        xlab='Tissue Samples', ylab='TPM values', main='TPM value of SLC24A2 gene 
        considering all the replicates of each tissue',
        names=c('Brain','Liver','Lung'),outline=F, notch=F, varwidth=T , col=myColors)
legend("topright", inset=.02, title="Tissues",
       c("Brain","Lung","Liver"), fill=c('mediumseagreen','lightskyblue','sienna1'), horiz=F, cex=0.35)
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-157-1.png)<!-- -->

TPM distribution in the three tissues, without considering the division
in replicates, of one of the most relevant gene: we can see that indeed
this gene is highly expressed in brain, while it is not expressed in
liver and lung. From this plot we can say that SLC24A2 can be considered
“statistically” over expressed in brain than in the other two tissues.

**Assess if the difference of expression of SLC24A2 gene is
statistically significant** between brain and the other two tissues by
considering all the samples of each tissue

``` r
tpm_brain <- assays(rse_brain)$TPM[32893,]
tpm_liver <- assays(rse_liver)$TPM[32893,]
tpm_lung <- assays(rse_lung)$TPM[32893,]
```

Since the TPM values of SLC24A2 gene in brain samples are independent
from the values in the samples of the other two tissues we use the
Mann–Whitney U test performing a one-sided test to check if:

-   H0: means of TMV values of the different tissues are the same
-   H1: mean of TMV values in brain is bigger than in the other tissues

**Mann–Whitney U test for TPM values of SLC24A2 gene in brain vs liver**

``` r
wilcox.test(tpm_brain, tpm_liver, paired=F, alternative='greater')
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  tpm_brain and tpm_liver
    ## W = 735681, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is greater than 0

Since the p-value is very low, \<2.2e-16, we can conclude that the mean
TPM value of SLC24A2 gene in brain is significantly higher than in
liver.

**Mann–Whitney U test for TPM values of SLC24A2 gene in brain vs lung**

``` r
wilcox.test(tpm_brain, tpm_lung, paired=F, alternative='greater')
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  tpm_brain and tpm_lung
    ## W = 1919805, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is greater than 0

Since the p-value is very low, \<2.2e-16, we can conclude that the mean
TPM value of SLC24A2 gene in brain is significantly higher than in lung.

## Functional enrichment analysis

``` r
library('enrichR')
setEnrichrSite("Enrichr")
websiteLive <- TRUE
```

### Up regulated genes in brain with respect to both lung and liver

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(both_brain, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-163-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-163-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-163-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(both_brain, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-165-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-165-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-165-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(both_brain, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-167-1.png)<!-- -->

### Up regulated genes in liver with respect to both lung and brain

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(both_liver, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-169-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-169-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-169-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(both_liver, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-171-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-171-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-171-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(both_liver, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-173-1.png)<!-- -->

### Up regulated genes in lung with respect to both liver and brain

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(both_lung, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-175-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-175-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-175-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(both_lung, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-177-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-177-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-177-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(both_lung, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-179-1.png)<!-- -->

# Comparison between the lists of up-regulated genes obtained starting from not-filtered vs filtered data

First we **load the list of significantly up regulated genes in each
tissue with respect to both the other two** obtained as result of DE
analysis starting from not-filtered data

``` r
rm(list=ls())

both_brain_nf <- read.table('both_brain_nf.txt', sep = '')
both_liver_nf <- read.table('both_liver_nf.txt', sep = '')
both_lung_nf  <- read.table('both_lung_nf.txt', sep = '')
```

Then we **load the list of significantly up regulated genes in each
tissue with respect to both the other two** obtained as result of DE
analysis starting from filtered data

``` r
both_brain_f <- read.table('both_brain_f.txt', sep = '')
both_liver_f <- read.table('both_liver_f.txt', sep = '')
both_lung_f  <- read.table('both_lung_f.txt', sep = '')
```

Now we have to **extract the names of those genes that are filtered out
since they are meaningless genes.** To do it we reload the data from
GTEx portal

``` r
rse_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)

rse_liver <-recount3::create_rse_manual(
  project = "LIVER",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)

rse_lung <- recount3::create_rse_manual(
  project = "LUNG",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "refseq",
  type = "gene"
)
```

We **extract the three previously selected replicates for each tissue**

``` r
rse_brain_selected <- rse_brain[,c(74,75,76)]
rse_liver_selected <- rse_liver[,c(74,75,78)]
rse_lung_selected <- rse_lung[,c(74,75,76)]
```

**This time we keep for each tissue only the rRNA genes, mitochondrial
genes, pseudogenes or genes of unknown type**

``` r
library(edgeR)
library(recount3)
canonical <- paste("chr", seq(1,22), sep="")
canonical <- c(canonical, "chrX", "chrY")

filtered_genes_brain <- rse_brain_selected[
  # Ribosomal RNA
  rowData(rse_brain_selected)$gbkey == 'rRNA' |
    # Pseudogenes
    rowData(rse_brain_selected)$gbkey == 'Gene' |
    # Non-canonical Chromosomes and Mitochondrial DNA
    !rowRanges(rse_brain_selected)@seqnames %in% canonical |
    is.na(rowData(rse_brain_selected)$gbkey),
  ]

filtered_genes_liver <- rse_liver_selected[
  # Ribosomal RNA
  rowData(rse_liver_selected)$gbkey == 'rRNA' |
    # Pseudogenes
    rowData(rse_liver_selected)$gbkey == 'Gene' |
    # Non-canonical Chromosomes and Mitochondrial DNA
    !rowRanges(rse_liver_selected)@seqnames %in% canonical |
    # NAs
    is.na(rowData(rse_liver_selected)$gbkey),
  ]

filtered_genes_lung <- rse_lung_selected[
  # Ribosomal RNA
  rowData(rse_lung_selected)$gbkey == 'rRNA' |
    # Pseudogenes
    rowData(rse_lung_selected)$gbkey == 'Gene' |
    # Non-canonical Chromosomes and Mitochondrial DNA
    !rowRanges(rse_lung_selected)@seqnames %in% canonical |
    # NAs
    is.na(rowData(rse_lung_selected)$gbkey),
  ]
```

**We transform the coverage counts into raw counts, fish out the count
tables and extract the gene names from it**

``` r
assays(filtered_genes_brain)$counts <- transform_counts(filtered_genes_brain)
assays(filtered_genes_liver)$counts <- transform_counts(filtered_genes_liver)
assays(filtered_genes_lung)$counts <- transform_counts(filtered_genes_lung)

counts_filt_lung <- assays(filtered_genes_lung)$counts
counts_filt_brain <- assays(filtered_genes_brain)$counts
counts_filt_liver<- assays(filtered_genes_liver)$counts

rownames(counts_filt_lung) <- rowData(filtered_genes_lung)$gene_name
rownames(counts_filt_liver) <- rowData(filtered_genes_liver)$gene_name
rownames(counts_filt_brain) <- rowData(filtered_genes_brain)$gene_name

filt_name_brain <- as.vector(rownames(counts_filt_brain))
filt_name_liver <- as.vector(rownames(counts_filt_liver))
filt_name_lung <- as.vector(rownames(counts_filt_lung))
```

## Analysis of up regulated genes obtained only starting from not-filtered data

Extract the gene names of those genes that result up regulated in each
tissue (vs both the other two) only when we don’t filter the initial
data

``` r
only_nf_brain <- both_brain_nf[which(!both_brain_nf[,1] %in% both_brain_f[,1]),]
only_nf_liver <- both_liver_nf[which(!both_liver_nf[,1] %in% both_liver_f[,1]),]
only_nf_lung <- both_lung_nf[which(!both_lung_nf[,1] %in% both_lung_f[,1]),]
```

``` r
length(only_nf_brain)
```

    ## [1] 589

``` r
length(only_nf_liver)
```

    ## [1] 434

``` r
length(only_nf_lung)
```

    ## [1] 402

There are **589 brain genes, 434 liver genes and 402 lung genes that are
up regulated in brain/liver/lung only in case of not filtering.**

For each tissue we compare these genes that have been found up regulated
only when starting from not filtered data with the list of genes that we
filtered out. This is done to understand how many of those genes that
result up regulated in a tissue - when we don’t filter - are genes that
can be filtered out since they are not interesting genes for a DE
analysis

``` r
length(only_nf_brain[which(!only_nf_brain %in% filt_name_brain)])
```

    ## [1] 76

``` r
length(only_nf_liver[which(!only_nf_liver %in% filt_name_liver)])
```

    ## [1] 1

``` r
length(only_nf_lung[which(!only_nf_lung %in% filt_name_lung)])
```

    ## [1] 64

We can see that, among all the genes that are up regulated in each
tissue only starting from not-filtered data, **only 76 brain genes, 1
liver gene and 64 lung genes are interesting DE genes (so no
mitochondrial/rRNA/pseudogenes/unknown gene types).**

**This means that lot of genes that are returned as significantly up
regulated in a tissue are not tissue-specific and so it is strongly
advised to filter the data at the beginning of the analysis otherwise
you have to do it on the final list of DE genes.**

## Functional enrichment analysis

We can perform enrichment analysis of those genes - 76 brain genes, 1
liver gene and 64 lung genes - that are up regulated in each tissue with
respect to both the other two only when we don’t filter the data and
that are genes of interest of that tissue.

``` r
library('enrichR')
setEnrichrSite("Enrichr")
websiteLive <- TRUE
```

### Brain

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(only_nf_brain[which(!only_nf_brain %in% filt_name_brain)], dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-191-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-191-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-191-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(only_nf_brain[which(!only_nf_brain %in% filt_name_brain)], dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-193-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-193-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-193-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(only_nf_brain[which(!only_nf_brain %in% filt_name_brain)], dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-195-1.png)<!-- -->

### Liver

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(only_nf_liver[which(!only_nf_liver %in% filt_name_liver)], dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-197-1.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(only_nf_liver[which(!only_nf_liver %in% filt_name_liver)], dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-199-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-199-2.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(only_nf_liver[which(!only_nf_liver %in% filt_name_liver)], dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-201-1.png)<!-- -->

### Lung

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(only_nf_lung[which(!only_nf_lung %in% filt_name_lung)], dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-203-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-203-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-203-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(only_nf_lung[which(!only_nf_lung %in% filt_name_lung)], dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-205-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-205-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-205-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(only_nf_lung[which(!only_nf_lung %in% filt_name_lung)], dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-207-1.png)<!-- -->

## Analysis of up regulated genes obtained only starting from filtered data

Extract the gene names of those genes that result up regulated in each
tissue (vs both the other two) only when we filter the initial data

``` r
only_f_brain <- both_brain_f[which(!both_brain_f[,1] %in% both_brain_nf[,1]),]
only_f_liver <- both_liver_f[which(!both_liver_f[,1] %in% both_liver_nf[,1]),]
only_f_lung <- both_lung_f[which(!both_lung_f[,1] %in% both_lung_nf[,1]),]
```

``` r
length(only_f_brain)
```

    ## [1] 51

``` r
length(only_f_liver)
```

    ## [1] 115

``` r
length(only_f_lung)
```

    ## [1] 28

There are **51 brain genes, 115 liver genes and 28 lung genes that are
up regulated in brain/liver/lung only when starting from filtered
data.** These may be genes that are significantly differentially
expressed in each tissue with respect to the other two that are hidden
by the presence of lot of ‘noise’ genes and so don’t result as DE when
we don’t filter the data (this is due to the fact that starting from a
different number of genes the normalized counts can be different and as
a consequence also the result of DE analysis may change).

## Functional enrichment analysis

We can perform enrichment analysis of those genes - 51 brain genes, 115
liver gene and 28 lung genes - that are up regulated in each tissue with
respect to both the other two only when we filter the data to understand
if those genes that are lost when we don’t filter the data are actually
specific of that tissue and so important to detect as DE.

### Brain

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(only_f_brain, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-211-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-211-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-211-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(only_f_brain, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-213-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-213-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-213-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(only_f_brain, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-215-1.png)<!-- -->

### Liver

### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(only_f_liver, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-217-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-217-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-217-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(only_f_liver, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-219-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-219-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-219-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(only_f_liver, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-221-1.png)<!-- -->

### Lung

#### Pathway enrichment analysis

``` r
dbs_pathway <- c("BioPlanet_2019", "WikiPathway_2021_Human", "KEGG_2021_Human")
if (websiteLive) {
    enriched_pathway <- enrichr(only_f_lung, dbs_pathway)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of BioPlanet 2019 database", enriched_pathway[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-223-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of WikiPathway 2021 Human database", enriched_pathway[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-223-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of KEGG 2021 Human database", enriched_pathway[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-223-3.png)<!-- -->

#### Ontologies enrichment analysis

``` r
dbs_ontologies <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
if (websiteLive) {
    enriched_ontologies <- enrichr(only_f_lung, dbs_ontologies)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological Process 2021 database", enriched_ontologies[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-225-1.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Molecular Function 2021 database", enriched_ontologies[[2]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-225-2.png)<!-- -->

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of GO Cellular Component 2021 database", enriched_ontologies[[3]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-225-3.png)<!-- -->

#### Cell types enrichment analysis

``` r
dbs_celltypes <- c("Human_Gene_Atlas")
if (websiteLive) {
    enriched_celltypes <- enrichr(only_f_lung, dbs_celltypes)
}
```

``` r
if (websiteLive) plotEnrich(title = "Enriched terms of Human Gene Atlas database", enriched_celltypes[[1]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
```

![](bulk_RNA-seq_analysis_without_filtering_files/figure-gfm/unnamed-chunk-227-1.png)<!-- -->

In conclusion, from this analysis, we can say that **filtering the data
is fundamental since lot of genes that result up regulated only starting
from not-filtered data are indeed mitochondrial genes, ribosomal protein
genes, psuedogenes or genes of unknown type of which we are not
interested when performing DE analysis.**

In fact only 3.9% of the genes found to be up regulated in brain
starting from not filtered data are not detected starting from filtered
one and are interesting genes, 0.07% in liver and 5.5% in lung. From
these numbers it is clear that **not-filtering the data doesn’t give any
advantage in terms of finding meaningful DE genes - with respect to
starting from filtered data - and furthermore it requires to perform a
filtering step on the final result of all the analysis.**
