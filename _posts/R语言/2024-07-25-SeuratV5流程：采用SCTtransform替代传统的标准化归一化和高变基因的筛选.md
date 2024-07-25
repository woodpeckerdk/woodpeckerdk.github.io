---
title: SeuratV5流程：采用SCTtransform替代传统的标准化归一化和高变基因的筛选
layout: post
post-image: "../assets/images/image_in_post/integration_seurat5.webp"
description: Seurat流程是单细胞分析的最基础的一步，几乎所有的分析都建立在其基础之上，目前Seurat从V4升级到了V5版本，数据结构增加了layer层的概念，标准分析流程也有一些小的改动，比如采用SCTransform包替代 NormalizeData(), ScaleData()和FindVariableFeatures()进行标准化。
tags:
- 单细胞
- Seurat
- 生信分析
- SCT
- 小鼠基因
---



# Seurat V5版小鼠单细胞流程（采用SCTransform包进行标准化）

和之前的分析流程基本类似，区别在于：
- 这次是小鼠单细胞
- 用SCTransform()替代 NormalizeData(), ScaleData()和 FindVariableFeatures()，按照Seurat官方的说法，SCTransform方法应该是更优的，具体暂不去深究
- 修改了之前笔记中存在的错误
- 按照vscode语法检查优化了代码书写格式.
## 加载包

```R
rm(list = ls()); gc() # nolint
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
```
## 批量读取标准10X数据
```r
full_name <- list.files("filtered_matrix", full.names = TRUE)
base_name <- basename(full_name)

#读取数据
sce_list <- list()
for (i in seq_along(full_name)) {
  print(paste0("正在读取第", i, "个样本数据"))
  counts <- Read10X(data.dir = full_name[i])
  sce_list[[i]] <- CreateSeuratObject(counts = counts,
                                      project = base_name[i],
                                      min.features = 200,
                                      min.cells = 3)
}
```
## 质控
- 计算线粒体比例的时候注意匹配的字符串是 "^mt-" 而不是“^MT-”，以防万一，最好先grep(）一下确认匹配表达式后再进行计算，否则有可能到后面出现线粒体比例为0的情况。

```r
# 计算线粒体比例
grep(pattern = "^mt-", rownames(sce_list[[1]]), value = TRUE)
for (i in seq_along(sce_list)) {
  print(paste0("正在计算第", i, "个样本的线粒体比例"))
  sc <- sce_list[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  sce_list[[i]] <- sc
  rm(sc)
}
```
分别绘制小提琴图可能更能观察每个样本的测序质量情况，虽然本质上与合并后再行质控是一样的，这里选择先绘制质控前小提琴图，根据测序质量调整参数到最佳（保证符合质控标准的前提下尽量保留合格的细胞）

```r
# 分别绘制小提琴图
violin_before <- list()
for (i in seq_along(sce_list)){
  violin_before[[i]] <- VlnPlot(sce_list[[i]],
                                features = c("nFeature_RNA",
                                             "nCount_RNA",
                                             "mt_percent"),
                                layer = "counts",
                                pt.size = 0.01,
                                ncol = 3)
}
violin_before_merge <- wrap_plots(plots = violin_before,
                                  nrow = length(sce_list),
                                  legend = "none")
# 保存图片，指定pdf画板的长宽。
ggsave("seurat流程/figures/violin_before_merge.pdf",
       plot = violin_before_merge,
       width = 15,
       height = 30)
# 质控
sce_list <- lapply(sce_list, function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 &
                nFeature_RNA < 5000 &
                mt_percent < 10 &
                nCount_RNA  < quantile(nCount_RNA, 0.97) &
                nCount_RNA > 500)
}
)
# 合并样本
sce <- merge(x = sce_list[[1]], y = sce_list[-1])
table(sce[[]]$orig.ident) # 统计细胞数
# 过滤后可视化,这是经过上一步QC后的小提琴图，展示目前的质控标准后保留的细胞。
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "mt_percent"),
        layer = "counts",
        pt.size = 0.01,
        ncol = 3)
```
## SCTransform转换矩阵

转换后，就要记住之后都是assay = "SCT",而不是之前常用的 assay = "RNA",这二者是同时存在的，目前"RNA"矩阵中只有“counts”一个layer,而“SCT”中是3个。其他的用法和之前是一样的，需要assay和layer参数时记得改过来就可以正常运行，以往的参数中有“slot”参数，相当于现在的“layer”。
```r
sce <- SCTransform(sce, vars.to.regress = "mt_percent", verbose = FALSE)

# seuratv5使用sct后的新标准流程
sce <- RunPCA(sce, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:30, verbose = FALSE)
sce <- RunTSNE(sce, dims = 1:30, verbose = FALSE)

sce <- FindNeighbors(sce, dims = 1:30, verbose = FALSE)
sce <- FindClusters(sce, verbose = FALSE, resolution = 0.5)
DimPlot(sce, reduction = "umap", label = TRUE)
DimPlot(sce, group.by = "orig.ident")
DimPlot(sce, reduction = "tsne")
```
## 差异基因分析
这里注意：在进行FindAllMarkers（）寻找差异基因之前，需要进行SCT数据预处理

```r
table(sce@meta.data$orig.ident)
sce <- PrepSCTFindMarkers(sce, assay = "SCT", verbose = TRUE)
all_markers <- FindAllMarkers(object = sce, test.use = "wilcox",
                              only.pos = TRUE,
                              logfc.threshold = 0.25) # FC值小于0.25的marker被删除
# 筛选出p < 0.05的marker
all_markers <- all_markers %>%
  select(gene, everything()) %>%
  subset(p_val < 0.05)
# 再筛选出每个cluster中log2FC值最大的top10
top10 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
write.csv(all_markers, "seurat流程/output/all_markers.csv", row.names = TRUE)
write.csv(top10, "seurat流程/output/top10.csv", row.names = TRUE)
# 可视化marker
DoHeatmap(sce, features = top10$gene) + NoLegend()

VlnPlot(sce, features = c("Adgre1"), assay = "SCT", layer = "data")
FeaturePlot(sce, features = c("Adgre1", "Cd14"), pt.size = 1)
```
## 细胞类型注释

细胞注释的时候首选DotPlot（）气泡图，能够快速注释，如果是比较熟悉的测序数据，可以直接经验性的将一些marker放进去验证，对于不熟悉的情况，比如这次小鼠的数据一些基因名称和人类不太一样，可以借助前一步获得的 top10_markers来进行注释，先确认哪些基因在特定cluster中特征性表达，再将基因去 [Cell Marker2.0](http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.html) 中查询属于哪些细胞类型或亚群。
### 确定各个celltype的marker基因
```r
genes_to_check <- c("Cd14", "Lyz2", "Itgam", "Adgre1", # 髓系细胞
                   "Fscn1", "Ccl22", "Cd209a", "Cd209c", "Clec9a",      "Ccser1", # DC细胞 # nolint
                   "Cd3e", "Cd4", "Cd8a",  # T淋巴细胞
                   "Klrb1c", "Itga2", "Klrk1", # NK细胞
                   "Col6a3", "C1s1", # Fibroblast
                   "Cd19", "Cd22", # B细胞
                   "Mcpt4", "Cma1", "Mrgprb1" # Mast cell
)


DotPlot(sce, features = genes_to_check, assay = "SCT",
        group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
### 注释

```r
celltype <- data.frame(clusterID = 0:18,
                       celltype = "unknown")
celltype[celltype$clusterID %in% c(0, 1, 2, 4, 5, 8, 10), 2] <- "Myeloid"
celltype[celltype$clusterID %in% c(3, 6, 7, 15, 17), 2] <- "T cell"
celltype[celltype$clusterID %in% c(9), 2] <- "NK cell"
celltype[celltype$clusterID %in% c(11), 2] <- "Fibroblast"
celltype[celltype$clusterID %in% c(12, 13, 14), 2] <- "Dendritic cell"
celltype[celltype$clusterID %in% c(16), 2] <- "Mast cell"
celltype[celltype$clusterID %in% c(18), 2] <- "B cell"
celltype
table(celltype$celltype)
sce_in <- sce
# 先加一列celltype所有值为空
sce_in@meta.data$celltype <-  "NA"
# 注释
for (i in seq_along(celltype$clusterID)) {
  sce_in@meta.data[which(sce_in@meta.data$seurat_clusters == celltype$clusterID[i]), # nolint
                   "celltype"] <- celltype$celltype[i]
}
table(sce_in@meta.data$celltype)
sce <- sce_in
```
注释完后记得保存所需的图

```r
DimPlot(sce, reduction = "tsne",
        group.by = "celltype",
        label = TRUE, repel = TRUE)
DimPlot(sce, reduction = "umap",
        group.by = "celltype",
        label = TRUE)

DefaultAssay(sce) <- "SCT"  # DefaultAssay设置为RNA意味着接下来的分析将基于原始值
```
注释完成后，可以再整理一下我们的元数据meta.data,这里添加一下sce$group实验组和对照组信息,为后续的分析做好准备，也可以根据情况加入其他的信息，比如各种评分、各种其他分组都可以使用相同的方法将相关信息添加进meta.data

```r
sce@meta.data %>% colnames()
sce$orig.ident %>% table()
sce$group <- ifelse(sce$orig.ident %in%
                      c("ctr_15_filtered", "ctr_16_filtered"), "control", "IKE")
sce$group %>% table()

saveRDS(sce, file = "seurat流程/rda/sce.rds")
```
### 注释后气泡图
```r
DotPlot(sce, features = genes_to_check, assay = "SCT",
        group.by = "celltype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

