melanoma_ATAC-RNA_integration
================
Konnie Guo
3/7/2022

``` r
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir ='C:/Users/kqguo/Desktop/khavari-lab')
library(GenomicRanges)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomeInfoDb

``` r
library(VennDiagram)
```

    ## Loading required package: grid

    ## Loading required package: futile.logger

``` r
library(RColorBrewer)
library(foreach)
```

## ATAC Gene Lists

Read in lists of genes corresponding to ATAC peaks, based on GREAT
analysis

``` r
files<-list.files("data/ATAC-seq", pattern = "*diffPeaks*")

diffpeaks<-foreach(file=files) %do% read.table(paste("data/ATAC-seq/",file,sep=""), header = T, sep = ",")
names(diffpeaks)<-basename(files) 

diffpeaks_p0.01<-lapply(diffpeaks,subset,pvalue<0.01)

diffpeaks_p0.05<-lapply(diffpeaks,subset,pvalue<0.05)
diffpeaks_coords<-lapply(diffpeaks_p0.05, '[', , c(2:4,1))

names(diffpeaks_coords)<-c("melanomavsMelanocyte_peaks_p0.05","COLOvsMelanocyte_peaks_0.05","PLX_peaks_0.05","shMITF_peaks_0.05","SKMEL5vsMelanocyte_peaks_0.05","WM2664vsMelanocyte_peaks_0.05")

#write out bed files for differential peaks with p-value < 0.01, for input into GREAT
# sapply(names(diffpeaks_coords), function (x) write.table(diffpeaks_coords[[x]], file=paste(x, "bed", sep="."), sep="\t", col.names = F, row.names = F, quote = F )   )


########################################################################################
# Alternate code to read out .bed files for GREAT - used if coords are out of bounds
########################################################################################

# diffpeaks_combined<-do.call("rbind",diffpeaks_coords)
# diffpeaks_GRanges<-makeGRangesFromDataFrame(diffpeaks_combined,
#                                 seqinfo = Seqinfo(genome="hg38"),
#                                 keep.extra.columns = T)
# diffpeaks_GRanges<-trim(diffpeaks_GRanges)
# 
# diffpeaks_coords<-data.frame(seqnames=seqnames(diffpeaks_GRanges),
#   starts=start(diffpeaks_GRanges)-1,
#   ends=end(diffpeaks_GRanges),
#   names=c(rep(".", length(diffpeaks_GRanges))),
#   scores=c(rep(".", length(diffpeaks_GRanges))),
#   strands=strand(diffpeaks_GRanges),
#   sample=diffpeaks_GRanges@ranges@NAMES,
#   peakID=diffpeaks_GRanges$GeneID)
# 
# melanomapeaks<-diffpeaks_coords[grep("melanomavsMelanocyte",diffpeaks_coords$sample),]
# PLXpeaks<-diffpeaks_coords[grep("PLX",diffpeaks_coords$sample),]
# shMITFpeaks<-diffpeaks_coords[grep("MITF",diffpeaks_coords$sample),]


# write.table(melanomapeaks[,c(1:3,8)],"results/ATAC-seq/melanoma_diffpeaks_p0.05.bed",sep="\t",col.names = F,row.names = F,quote = F)
# write.table(PLXpeaks[,c(1:3,8)],"results/ATAC-seq/PLX_diffpeaks_p0.05.bed",sep="\t",col.names = F,row.names = F,quote = F)
# write.table(shMITFpeaks[,c(1:3,8)],"results/ATAC-seq/shMITF_diffpeaks_p0.05.bed",sep="\t",col.names = F,row.names = F,quote = F)
```

## RNA-seq vs ATAC GREAT genes

``` r
RNA_PLX<-read.csv("results/RNA-seq/DEG_tables/DGE_drug.csv")
RNA_PLX_filt<-subset(RNA_PLX,abs(RNA_PLX$log2FoldChange)>=1.5)
RNA_PLX_filt<-RNA_PLX_filt[!is.na(RNA_PLX_filt$SYMBOL),]
ATAC_PLX_GREAT<-read.table("results/ATAC-seq/GREAT/genes_GREATanalysis_PLXtreated_p0.05_prefilt.txt", sep = "\t", header = F)

genes<-as.data.frame(intersect(RNA_PLX$SYMBOL,ATAC_PLX_GREAT$V1))

draw.pairwise.venn(nrow(RNA_PLX_filt),nrow(ATAC_PLX_GREAT), nrow(genes), category = c("RNA-seq", "ATAC-seq"), 
             lwd=2,
             fill=brewer.pal(3, "Pastel2")[c(1, 2)],
             cex = 1.5,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 1.5,
             cat.fontface = "bold",
             cat.fontfamily = "sans"
)
```

![](melanoma_ATAC-RNA_integration_files/figure-gfm/overlap_PLX-1.png)<!-- -->

    ## (polygon[GRID.polygon.1], polygon[GRID.polygon.2], polygon[GRID.polygon.3], polygon[GRID.polygon.4], text[GRID.text.5], text[GRID.text.6], text[GRID.text.7], text[GRID.text.8], text[GRID.text.9])

``` r
# write.table(genes,"results/commongenes_PLX_RNA_ATAC.txt", col.names = F, row.names = F, quote = F)
```

``` r
RNA_shMITF<-read.csv("results/RNA-seq/DEG_tables/DGE_KD.csv")
RNA_shMITF_filt<-subset(RNA_shMITF,abs(RNA_shMITF$log2FoldChange)>=1.5)
RNA_shMITF_filt<-RNA_shMITF_filt[!is.na(RNA_shMITF_filt$SYMBOL),]
ATAC_shMITF_GREAT<-read.table("results/ATAC-seq/GREAT/genes_GREATanalysis_shMITF_p0.05_prefilt.txt", sep = "\t", header = F)

genes<-as.data.frame(intersect(RNA_shMITF$SYMBOL,ATAC_shMITF_GREAT$V1))

draw.pairwise.venn(nrow(RNA_shMITF_filt),nrow(ATAC_shMITF_GREAT), nrow(genes), category = c("RNA-seq", "ATAC-seq"), 
             lwd=2,
             fill=brewer.pal(3, "Pastel2")[c(1, 2)],
             cex = 1.5,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 1.5,
             cat.fontface = "bold",
             cat.fontfamily = "sans"
)
```

![](melanoma_ATAC-RNA_integration_files/figure-gfm/overlap_shMITF-1.png)<!-- -->

    ## (polygon[GRID.polygon.10], polygon[GRID.polygon.11], polygon[GRID.polygon.12], polygon[GRID.polygon.13], text[GRID.text.14], text[GRID.text.15], lines[GRID.lines.16], text[GRID.text.17], text[GRID.text.18], text[GRID.text.19])

``` r
# write.table(genes,"results/commongenes_shMITF_RNA_ATAC.txt", col.names = F, row.names = F, quote = F)
```
