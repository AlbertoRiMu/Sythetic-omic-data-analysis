# Leer los nombres de los genes desde el archivo .tsv
gene_names <- read.table("cancerGeneList.tsv", header=FALSE, stringsAsFactors = FALSE, sep="\t", skip=1)[,1]

# Cargar las librerías necesarias
library(tidyverse)

# Creación de las tres matrices
n_genes <- length(gene_names)
n_patients <- 500
expr_matrix <- matrix(rnorm(n_genes * n_patients), nrow = n_genes)
met_matrix <- matrix(rnorm(n_genes * n_patients), nrow = n_genes)
prot_matrix <- matrix(rnorm(n_genes * n_patients), nrow = n_genes)

# Nombres de las filas y columnas de las matrices
rownames(expr_matrix) <- gene_names
rownames(met_matrix) <- gene_names
rownames(prot_matrix) <- gene_names
colnames(expr_matrix) <- paste0("patient", 1:n_patients)
colnames(met_matrix) <- paste0("patient", 1:n_patients)
colnames(prot_matrix) <- paste0("patient", 1:n_patients)

# PCA para la matriz de expresión génica
pca_expr <- prcomp(t(expr_matrix), scale = TRUE)

# PCA para la matriz de metilación
pca_met <- prcomp(t(met_matrix), scale = TRUE)

# PCA para la matriz de proteómica
pca_prot <- prcomp(t(prot_matrix), scale = TRUE)

# Visualización de los componentes principales en gráficos de dispersión

ggplot(as.data.frame(pca_expr$x), aes(PC1, PC2)) +
  geom_point(color = "blue") +
  xlab(paste0("PC1 (", round(pca_expr$sdev[1] / sum(pca_expr$sdev) * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(pca_expr$sdev[2] / sum(pca_expr$sdev) * 100, 2), "%)")) +
  ggtitle("PCA para la matriz de expresión génica")

ggplot(as.data.frame(pca_met$x), aes(PC1, PC2)) +
  geom_point(color = "red") +
  xlab(paste0("PC1 (", round(pca_met$sdev[1] / sum(pca_met$sdev) * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(pca_met$sdev[2] / sum(pca_met$sdev) * 100, 2), "%)")) +
  ggtitle("PCA para la matriz de metilación")

ggplot(as.data.frame(pca_prot$x), aes(PC1, PC2)) +
  geom_point(color = "green") +
  xlab(paste0("PC1 (", round(pca_prot$sdev[1] / sum(pca_prot$sdev) * 100, 2), "%)")) +
  ylab(paste0("PC2 (",round(pca_prot$sdev[2] / sum(pca_prot$sdev) * 100, 2), "%)")) +
  ggtitle("PCA para la matriz de proteómica")


# Matriz de correlación omicas
cor_matrix_expr_met <- cor(cbind(t(expr_matrix), t(met_matrix)))
cor_matrix_expr_prot <- cor(cbind(t(expr_matrix), t(prot_matrix)))
cor_matrix_prot_met <- cor(cbind(t(prot_matrix), t(met_matrix)))

# Heatmap de correlación

library(gplots)

jpeg("correlation_heatmap_expr_met.jpeg")
heatmap.2(cor_matrix_expr_met, trace="none")
dev.off()

## Lasso model de las tres ómicas integradas 

# Load the necessary packages

library(glmnet)


# Merge the three matrices into a single data frame
merged_matrix <- cbind(expr_matrix, met_matrix, prot_matrix)

# Normalize the data
normalized_data <- scale(merged_matrix)

# Assuming the data has the first 350 patients with cancer and the rest 150 without cancer
y1 <- rep(1, 350)
y0 <- rep(0, 150)
y <- c(y1, y0)
y_long <- c(y, y, y)


# Perform feature fusion using the Lasso regularization method
lasso_model <- cv.glmnet(x = t(normalized_data), y = y_long, alpha = 1, nfolds = 5)

# Select the genes that have non-zero coefficients
selected_genes <- lasso_model[["nzero"]]


# These are the DEGs identified by the Lasso regularization method
DEGs <- merged_matrix[selected_genes, ]
DEGs_Unique <- unique(DEGs, fromLast = TRUE)



## Análisis de expresión diferencial con DESEQ2

library(DESeq2)

# Generate the synthetic dataset 
library(compcodeR)

data <- generateSyntheticData(dataset = "Cancer_genes_Expr", 
                              n.vars = 1088, 
                              samples.per.cond = 250, 
                              n.diffexp = 110,
                              seqdepth = 1e+07,
                              minfact = 0.7,
                              maxfact = 1.4,
                              relmeans = "auto",
                              dispersions = "auto",
                              fraction.upregulated = 1,
                              between.group.diffdisp = FALSE,
                              filter.threshold.total = 1,
                              filter.threshold.mediancpm = 0,
                              fraction.non.overdispersed = 0,
                              random.outlier.high.prob = 0,
                              random.outlier.low.prob = 0,
                              single.outlier.high.prob = 0,
                              single.outlier.low.prob = 0,
                              effect.size = 1.5,
                              output.file = NULL)

# Crear un objeto de tipo DESeqDataSet a partir de los datos generados con compcodeR 

condition <- c(rep("cancer", 250), rep("no_cancer", 250))
colData <- data.frame(condition=condition)
dds <- DESeqDataSetFromMatrix(countData = data@count.matrix,
                              colData = colData,
                              design = ~condition)

# Normalizar datos

dds <- DESeq(dds)

# Análisis de expresión diferencial 
res <- results(dds)

# Convertir resultados a data.frame
res_df <- as.data.frame(res)
rownames(res_df) <- gene_names

library(ggplot2)

# Crear una columna para identificar a cada gene
res_df$regulacion <- ifelse(res_df$log2FoldChange > 1.5 & res_df$padj < 0.05, "Alza",
                            ifelse(res_df$log2FoldChange < -1.5 & res_df$padj < 0.05, "Reprimido", "No regulado"))

#Agrega la columna "GeneName" a res_df

res_df$GeneName <- rownames(res_df)
# Crear el plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulacion)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(0, 100)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Padj") +
  ggtitle("Volcano plot") +
  scale_color_manual(values = c("red" , "black" )) +
  geom_text(data = head(res_df[res_df$regulacion=="Alza",], 10), aes(label = GeneName), size = 3, hjust = 0, vjust = 0, color = "red")
  theme_classic()

  
## Enriquecimiento funcional 
  
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)

# Nos quedamos con los genes regulados al alza solamente, vamos a tomar un fold change de 1 esta vez

up_reg <- res_df[res_df$log2FoldChange > 1,]
genes_to_test <- rownames(up_reg)

# Biological process
GO_results_BP <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

df_BP <- as.data.frame(GO_results_BP)

fit_BP <- plot(barplot(GO_results_BP, showCategory = 10, font.size = 10))
fit_BP


# Molecular functions
GO_results_MF <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF")

df_MF <- as.data.frame(GO_results_MF)

fit_MF <- plot(barplot(GO_results_MF, showCategory = 10, font.size = 10))
fit_MF


# Celular Components

GO_results_CC <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC")

df_CC <- as.data.frame(GO_results_CC)

fit_CC <- plot(barplot(GO_results_CC, showCategory = 10, font.size = 10))
fit_CC

# Save results
jpeg("BP.jpeg", res = 250, width = 1200, height = 1000)
print(fit_BP)
dev.off()

jpeg("MF.jpeg", res = 250, width = 1200, height = 1000)
print(fit_MF)
dev.off()

jpeg("CC.jpeg", res = 250, width = 1200, height = 1000)
print(fit_CC)
dev.off()