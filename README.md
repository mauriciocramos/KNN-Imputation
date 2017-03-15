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
rownames(gene.expression.values) <- paste0("G",1:nrow(gene.expression.values))
colnames(gene.expression.values) <- paste0("S",1:ncol(gene.expression.values))
gene.expression.values[] <- lapply(gene.expression.values[], function (x) as.numeric(levels(x)[x]))
gene.expression.values <- as.matrix(gene.expression.values)#[1:20, 1:4]
str(gene.expression.values)
```

    ##  num [1:2308, 1:63] 0.773 -2.438 -0.483 -2.721 -1.217 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:2308] "G1" "G2" "G3" "G4" ...
    ##   ..$ : chr [1:63] "S1" "S2" "S3" "S4" ...

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
hist(gene.expression.values, freq = FALSE, pch=20, col=rgb(0,0,0,0.03), cex=0.1,
     main = "Histogram of Gene Expression Values", xlab = "Gene Expression Values")
## Normal Curve
x <- gene.expression.values
mean <- mean(gene.expression.values, na.rm = TRUE)
sd <- sd(gene.expression.values, na.rm = TRUE)
curve(dnorm(x, mean, sd), add = TRUE, col = "blue", lwd = 2)
## Kernel density estimation
lines(density(gene.expression.values, na.rm = TRUE), col = "red", lwd = 2)
## Legend
legend("topright", lwd = 2, col = c("blue", "red"), text.col = c("blue", "red"), legend = c("Normal Curve", "Kernel Density"))
```

![](README_files/figure-markdown_github/Histogram%20of%20Gene%20Expression%20Values-1.png)

### Matrix Plot

``` r
matplot(gene.expression.values, pch=20, cex=0.1, col=rgb(0,0,0,0.1),
        xlim=c(1,nrow(gene.expression.values)),
        main = "Gene Expression Values by Gene", 
        xlab = "Genes", ylab = "Gene Expression Values")
```

![](README_files/figure-markdown_github/Gene%20Expression%20Values%20by%20Gene-1.png)

### Transposed Matrix Plot

``` r
matplot(t(gene.expression.values), pch=20, cex=1, col=rgb(0,0,0,0.1),
        xlim=c(1,ncol(gene.expression.values)),
        main = "Gene Expression Values by Samples", 
        xlab = "Samples", ylab = "Gene Expression Values")
```

![](README_files/figure-markdown_github/Gene%20Expression%20Values%20by%20Samples-1.png)

### Matrix image

``` r
image(x = 1:ncol(gene.expression.values), y = 1:nrow(gene.expression.values), z = t(gene.expression.values),
      ylim=c(nrow(gene.expression.values),1), col=1, xlab = "Samples", ylab = "Genes", main="Missing Gene Expression Values")
```

![](README_files/figure-markdown_github/Gene%20Expression%20Matrix%20image-1.png)

### Hierarchical cluster of the Samples

``` r
(hclust<- hclust(dist(t(gene.expression.values))))
```

    ## 
    ## Call:
    ## hclust(d = dist(t(gene.expression.values)))
    ## 
    ## Cluster method   : complete 
    ## Distance         : euclidean 
    ## Number of objects: 63

``` r
plot(hclust, xlab = "Sample distance", main = "Hierarchical cluster of the Samples")
```

![](README_files/figure-markdown_github/Hierarchical%20cluster%20of%20the%20Samples-1.png)

### Hierarchical cluster of the Genes

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

``` r
plot(hclust, hang=0.1, labels = FALSE, xlab = "Gene distance", main = "Hierarchical cluster of the Genes")
```

![](README_files/figure-markdown_github/Hierarchical%20cluster%20of%20the%20Genes-1.png)

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

Matrix of imputed values only

``` r
missing.points <- missing
missing.points[!missing] <- NA
missing.points[missing] <- imputed$data[missing]
```

``` r
matplot(imputed$data, pch=20, cex=0.1, col=rgb(0,0,0,0.1),
        xlim=c(1,nrow(imputed$data)),
        main = "Gene Expression Values by Gene after imputation", 
        xlab = "Genes", ylab = "Gene Expression Values by Gene")
matpoints(missing.points, pch=20, cex=0.1, col=rgb(1,0,0,0.5))
legend("bottomright", pch = c(20,20), col = c("black", "red"), text.col = c("black", "red"), legend = c("former", "imputed"))
```

![](README_files/figure-markdown_github/Gene%20Expression%20Values%20by%20Gene%20after%20imputation-1.png)

``` r
matplot(t(imputed$data), pch=20, cex=1, col=rgb(0,0,0,0.1),
        xlim=c(1,ncol(imputed$data)),
        main = "Gene Expression Values by Sample after imputation",
        xlab = "Samples", ylab = "Gene Expression Values")
matpoints(t(missing.points), pch=20, cex=1, col=rgb(1,0,0,0.5))
legend("bottomright", pch = c(20,20), col = c("black", "red"), text.col = c("black", "red"), legend = c("former", "imputed"))
```

![](README_files/figure-markdown_github/Gene%20Expression%20Values%20by%20Sample%20after%20imputation-1.png)
