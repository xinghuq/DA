
# DA: Discriminant analysis

For a more detailed tutorial on population structure inference please visit: https://rpubs.com/xinghuq/592149

#  DAPC, DAKPC, LFDA, LFDAKPC, KLFDA

This package implements five supervised learning approaches that are suitable for ecological and evolutionary inference both visually and statistically. Five approaches in our vignette: Linear discriminant analysis of principal components (DAPC), linear discriminant analysis of kernel principal components (DAKPC), local (Fisher) linear discriminant analysis (LFDA), local (Fisher) discriminant analysis of kernel principal components (LFDAKPC) and kernel local (Fisher) linear discriminant analysis (KLFDA). 

Welcome any [feedback](https://github.com/xinghuq/DA/issues) and [pull request](https://github.com/xinghuq/DA/pulls).  


## Install the package from github:
```{R}
library(devtools)

install_github("xinghuq/DA")

library("DA")
```

## Basic examples using the irish data

```{R}
x <- iris[,-5] # this matrix contains all the predictors to be transformed
y <- iris[,5] # this should be a vector that represents different classes

```

###  Discriminant analysis of principal components (DAPC)

```{r fig1, fig.height = 5, fig.width = 10, fig.align = "center"}
library(adegenet)
iris_dapc=dapc(iris[,-5],grp=iris[,5],n.pca=3, n.da=3)

#plot the data projection on the components
library(plotly)
   cols=rainbow(length(unique(iris[,5])))
   p1 <- plot_ly(as.data.frame(iris_dapc$ind.coord), x =iris_dapc$ind.coord[,1], y =iris_dapc$ind.coord[,2], color = iris[,5],colors=cols[iris[,5]],symbol = iris[,5],symbols = 1:3L) %>% 
     add_markers() %>%
     layout(scene = list(xaxis = list(title = 'LDA1'),
                         yaxis = list(title = 'LDA2')))
p1
```

### Discriminant analysis of kernel principal components (DAKPC)

```{r fig2, fig.height = 5, fig.width = 10, fig.align = "center"}
iris_ldakpc=LDAKPC(iris[,-5],y=iris[,5],n.pc=3)

 p2 <- plot_ly(as.data.frame(iris_ldakpc$LDs), x =iris_ldakpc$LDs[,1], y =iris_ldakpc$LDs[,2], color = iris[,5],colors=cols[iris[,5]],symbol = iris[,5],symbols = 1:3L) %>% 
     add_markers() %>%
     layout(scene = list(xaxis = list(title = 'LDA1'),
                         yaxis = list(title = 'LDA2')))
p2
```


### Local Fisher Discriminant Analysis(LFDA)

```{r fig3, fig.height = 5, fig.width = 10, fig.align = "center"}
iris_lfda=LFDA(iris[,-5],y=iris[,5],r=3,tol=1)

 p3 <- plot_ly(as.data.frame(iris_lfda$Z), x =iris_lfda$Z[,1], y =iris_lfda$Z[,2], color = iris[,5],colors=cols[iris[,5]],symbol = iris[,5],symbols = 1:3L) %>% 
     add_markers() %>%
     layout(scene = list(xaxis = list(title = 'LDA1'),
                         yaxis = list(title = 'LDA2')))
p3

```
##  Local (fisher) discriminant analysis of kernel principal components (LFDAKPC)
```{r fig4, fig.height = 5, fig.width = 10, fig.align = "center"}
iris_lfdakpc=LFDAKPC(iris[,-5],y=iris[,5],n.pc=3,tol=1)

 p4 <- plot_ly(as.data.frame(iris_lfdakpc$$LFDAKPC$Z), x =iris_lfdakpc$$LFDAKPC$Z[,1], y =iris_lfdakpc$$LFDAKPC$Z[,2], color = iris[,5],colors=cols[iris[,5]],symbol = iris[,5],symbols = 1:3L) %>% 
     add_markers() %>%
     layout(scene = list(xaxis = list(title = 'LDA1'),
                         yaxis = list(title = 'LDA2')))
p4

```

### Kernel Local Fisher Discriminant Analysis(KLFDA)
The default kernel is polydot(degree = 1, scale = 1, offset = 1). Users can set the kernel based on their own purpose.
 
```{r fig5, fig.height = 5, fig.width = 10, fig.align = "center"}
iris_klfda=klfda_1(as.matrix(iris[,-5]),as.matrix(iris[,5]),r=3,tol=1E-10,prior = NULL)


p5 <- plot_ly(as.data.frame(iris_klfda$Z), x =iris_klfda$Z[,1], y =iris_klfda$Z[,2], color = iris[,5],colors=cols[iris[,5]],symbol = iris[,5],symbols = 1:3L) %>% 
     add_markers() %>%
     layout(scene = list(xaxis = list(title = 'LDA1'),
                         yaxis = list(title = 'LDA2')))
p5
```
Note we did not show the discrimination results,users can look into the discriminated classes and their posterior possibility from the results.



## Citation

Xinghu Qin, Mary Wu, T. Ryan Lock, Robert L. Kallenbach. 2020. DA: Ecological and evolutionary inference using machine learning approaches both visually and statistically.R package version, 0.5.

