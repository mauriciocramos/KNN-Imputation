KNN-Imputation
================
Maurício Collaça Ramos
2017-03-15

Imputation of missing data using k-nearest neighbor averaging (K-NN)
--------------------------------------------------------------------

**Data source**

Khan microarray data with random missing values

**Description**

A text file containing the Khan micorarray data with random missing values introduced for illustrative purposes

**Usage**

``` r
library(impute)
data(khanmiss, verbose = TRUE)
```

**Format**

The data set khanmiss consists of 2310 rows and 65 columns. Row 1 has the sample labels, Row 2 has the class labels. The remaining rows are gene expression. Column 1 is a dummy gene number. Column 2 is the gene name. Remaining columns are gene expression.

Please note that this dataset was derived from the original by introducing some random missing values purely for the purpose of illustration.

**Source**

Khan, J. and Wei, J.S. and Ringner, M. and Saal, L. and Ladanyi, M. and Westermann, F. and Berthold, F. and Schwab, M. and Antonescu, C. and Peterson, C. and and Meltzer, P. (2001) Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural network. Nature Medicine 7, 673-679.

**References**

Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu (2002). Diagnosis of multiple cancer types by shrunken centroids of gene expression PNAS 99: 6567-6572. Available at www.pnas.org

Preparing data
--------------

``` r
gene.expression.values <- khanmiss[-1, -(1:2)]
rownames(gene.expression.values) <- 1:nrow(gene.expression.values)
colnames(gene.expression.values) <- 1:ncol(gene.expression.values)
gene.expression.values[] <- lapply(gene.expression.values[], function (x) as.numeric(levels(x)[x]))
gene.expression.values <- as.matrix(gene.expression.values)
str(gene.expression.values)
```

    ##  num [1:2308, 1:63] 0.773 -2.438 -0.483 -2.721 -1.217 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:2308] "1" "2" "3" "4" ...
    ##   ..$ : chr [1:63] "1" "2" "3" "4" ...

Checking missing values
-----------------------

``` r
missing <- is.na(gene.expression.values)
sum(missing)
```

    ## [1] 1282

``` r
sprintf("%1.2f%%", 100*mean(missing))
```

    ## [1] "0.88%"

Data patterns
-------------

### Histogram

``` r
par(mfrow=c(1, 1), mar=c(4, 4, 2, 0), xaxt="s", yaxt="s", xaxs="r", yaxs="r")
hist(gene.expression.values, pch=20, col=rgb(0,0,0,0.05), cex=0.1,
     main = "Histogram of Gene Expression Values", xlab = "Gene Expression Values")
```

![](README_files/figure-markdown_github/Histogram%20of%20Gene%20Expression%20Values-1.png)

### Matrix Plot

``` r
matplot(gene.expression.values, pch=20, cex=0.1, col=rgb(0,0,0,0.1),
        xlab = "Genes", ylab = "Gene Expression Values", xlim=c(1,nrow(gene.expression.values)),
        main = "Gene Expression Values by Gene")
```

![](README_files/figure-markdown_github/Gene%20Expression%20Values%20by%20Gene-1.png)

### Transposed Matrix Plot

``` r
matplot(t(gene.expression.values), pch=20, cex=1, col=rgb(0,0,0,0.1),
        xlab = "Samples", ylab = "Gene Expression Values", xlim=c(1,ncol(gene.expression.values)),
        main = "Gene Expression Values by Samples")
```

![](README_files/figure-markdown_github/Gene%20Expression%20Values%20by%20Samples-1.png)

### Matrix image

``` r
image(x = 1:ncol(gene.expression.values), y = 1:nrow(gene.expression.values), z = t(gene.expression.values),
      ylim=c(nrow(gene.expression.values),1), col=1, xlab = "Samples", ylab = "Gene", main="Missing Gene Expression Values")
```

![](README_files/figure-markdown_github/Gene%20Expression%20Matrix%20image-1.png)

### Hierarchical clustering

``` r
(hclust<- hclust(dist(gene.expression.values)))
```

    ## 
    ## Call:
    ## hclust(d = dist(gene.expression.values))
    ## 
    ## Cluster method   : complete 
    ## Distance         : euclidean 
    ## Number of objects: 2308

### Cluster Dendogram

``` r
plot(hclust, hang=0, labels= FALSE)
```

![](README_files/figure-markdown_github/Gene%20Expression%20Cluster%20Dendogram-1.png)

### Heatmap

``` r
heatmap(gene.expression.values)
```

![](README_files/figure-markdown_github/Gene%20Expression%20Heatmap-1.png)

Imputation
----------

``` r
if(exists(".Random.seed")) rm(.Random.seed)
imputed <- impute.knn(gene.expression.values)
```

    ## Cluster size 2308 broken into 1450 858 
    ## Done cluster 1450 
    ## Done cluster 858

``` r
matplot(imputed$data, pch=20, cex=0.1, col=rgb(missing,0,0,(missing*0.9)+0.1),
        xlab = "Genes", ylab = "Gene Expression Values", xlim=c(1,nrow(imputed$data)),
        main = "Gene Expression Values by Gene")
legend("bottomright", pch = c(20,20), col = c("black", "red"), text.col = c("black", "red"), legend = c("former", "imputed"))
```

![](README_files/figure-markdown_github/Imputed%20Gene%20Expression%20Values%20by%20Gene-1.png)

``` r
matplot(t(imputed$data), pch=20, cex=1, col=rgb(missing,0,0,(missing*0.9)+0.1),
        xlab = "Samples", ylab = "Gene Expression Values", xlim=c(1,ncol(imputed$data)),
        main = "Gene Expression Values by Sample")
legend("bottomright", pch = c(20,20), col = c("black", "red"), text.col = c("black", "red"), legend = c("former", "imputed"))
```

![](README_files/figure-markdown_github/Imputed%20Gene%20Expression%20Values%20by%20Sample-1.png)

`rng.seed` is the seed used for the random number generator (default 362436069) for reproducibility.

``` r
imputed$rng.seed # should be 362436069
```

    ## [1] 362436069

`rng.state` is the state of the random number generator, if available, prior to the call to set.seed. Otherwise, it is NULL. If necessary, this can be used in the calling code to undo the side-effect of changing the random number generator sequence.

``` r
imputed$rng.state # should be NULL
```

    ## NULL

**Restoring the random number generation sequence changed before impute.knn()**

``` r
set.seed(12345)
saved.state <- .Random.seed
# Assuming all goes well with no guarantees in case of error...
.Random.seed <- impute.knn(gene.expression.values)$rng.state
```

    ## Cluster size 2308 broken into 1450 858 
    ## Done cluster 1450 
    ## Done cluster 858

**.Random.seed validation**

``` r
sum(saved.state - .Random.seed) # should be zero!
```

    ## [1] 0
