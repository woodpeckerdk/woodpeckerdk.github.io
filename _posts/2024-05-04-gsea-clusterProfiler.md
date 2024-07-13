---
title: "clusterProfiler包进行GSEA分析"
layout: post
post-image: "https://s2.loli.net/2024/05/04/IxqtJrfM34o82AW.png"
tags:
- R语言
- 富集分析
- clusterProfiler
- GSEA
---

clusterProfiler包是南方医科大学余光创教授开发的R包，用于分析基因富集分析结果，包括GO富集分析、KEGG富集分析、GSEA富集分析等。本文将介绍如何使用clusterProfiler包进行GSEA分析。

# 通过clusterProfiler包进行GSEA分析
## 环境的配置

清空环境变量，释放内存，加载包

```{r}
rm(list = ls());gc()
library(clusterProfiler) 
library(readxl) #读取excel文件
```

在MSigDB数据库下载需要的gmt文件，其实就是基因集名称和所包含的基因对应的表，因为每个基因集包含的基因数目不一样，所以不能以data.frame这样的格式读取，如果将gmt文件以文本格式读取进来会是一个列表。clusterProfiler包定义了一个`read.gmt`函数，可以读取gmt文件，生成每一个基因对应自己的基因集名称的数据框。也就是gene这一列是唯一的，但是基因集term这一列可以是重复的，两列对齐。

![](images/clipboard-3950127237.png)

```{r}
res <-read.gmt("c5.go.bp.v2023.2.Hs.symbols.gmt")
head(res)
```

## 数据的读取

读取差异分析结果文件,并进行排序,需要注意的是，最后输入的genelist是一个有name的向量，而不是一列数据框。一定注意,避免踩坑。

```{r}
geneList<-read_excel("差异分析.xlsx")
genelist <- geneList$log2FoldChange
names(genelist) <- geneList$id
genelist <- sort(genelist,decreasing = T)
class(genelist)
```

开始分析，这里面的参数，genelist就是一个数值型向量，拥有geneid作为name,`TERM2GENE= res`则是上面我们由gmt文件生成的两列数据框,起到注释的作用。

```{r}
GSEA_enrichment <- GSEA(genelist,#待富集的基因列表
                    TERM2GENE = res,  #基因集
                    pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                    pAdjustMethod = 'BH',
                    minGSSize=5,
                    maxGSSize = 500,
                    eps = 1e-10)
```

分析完成，可以将结果保存为csv文件，查看富集到的通路上调下调情况

```{r}
GSEA_enrichment@result %>% head()
write.csv(GSEA_enrichment@result,"result.csv")
```

分析完成后进行可视化,这里示例取3个基因集进行可视化，这是没有进行排序的，实际操作中自己根据p值和FoldChange值进行筛选和选择展示的通路，可以是一条或多条。

```{r}
geneSET <- GSEA_enrichment@result$ID[1:3]
geneSET
```

## 可视化展示

```{r}
library(enrichplot)
gseaplot2(
  GSEA_enrichment, #gseaResult object，即GSEA结果
  geneSetID = geneSET,#富集的ID编号
  title = "", #标题
  color = "green",#GSEA线条颜色
  base_size = 11,#基础字体大小
  rel_heights = c(1.5, 0.5,0),#副图的相对高度
  subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
  pvalue_table = T, #是否添加 pvalue table
  ES_geom = "line" #running enrichment score用先还是用点ES_geom = "dot"
)
```
