library(HMP2Data) #Bioconductor 
library(phyloseq) #Bioconductor
library(SummarizedExperiment) #Bioconductor
library(MultiAssayExperiment) #Bioconductor
library(dplyr)
library(ggplot2)
library(UpSetR)
library(DESeq2) #Bioconductor
library(igraph)
library(Hmisc)
library(data.table)


#Extraemos datos del proyecto IBD del Human Microbiome Project (HMP2) usando la 
#librería de HMP2Data
#------------------------------------------------------------------------------
# 1.- Cargamos todos los datos incluyendo Metadatos clínicos, Datos Taxonómicos, y Matriz de conteos de OTUs.

IBD <- IBD16S() #%>% tax_glom(., taxrank = "Genus")

clinical_data <- sample_data(IBD) %>%  as.data.frame()
tax_data <- tax_table(IBD) %>%  as.data.frame()
otus <- otu_table(IBD) %>%  as.data.frame()

otus[1:5,1:5]
clinical_data[1:5,1:5]
tax_data[1:5,1:6]
#-------------------------------------------------------------------------------
# 2. Preprocesamiento y normalización de datos
# Filtramos otus por número de conteos
threshold <- round(dim(otus)[2]/20) #Umbral de otus !=0 del 5%
print(paste0("Threshold: ",threshold))
otus <- otus[rowSums(otus != 0) >= threshold,  ]

# Normalización

#Primero obtenemos los ids de las muestras que pertenecen a cada condición
CD_samples_ids <- clinical_data[clinical_data$diagnosis == "CD",] %>% pull(sample_id)
UC_samples_ids <- clinical_data[clinical_data$diagnosis == "UC",] %>% pull(sample_id)
Ctrl_samples_ids <-  clinical_data[clinical_data$diagnosis == "nonIBD",] %>% pull(sample_id)

#Damos formato al experimento para DeSEQ2

coldata_CD <- cbind(CD_samples_ids, rep("CD", length(CD_samples_ids)))
coldata_UC <- cbind(UC_samples_ids, rep("UC", length(UC_samples_ids)))
coldata_ctrl <- cbind(Ctrl_samples_ids, rep("NonIBD", length(Ctrl_samples_ids)))

colnames(coldata_CD) <- c("sample_id", "condition")
colnames(coldata_UC) <-  c("sample_id", "condition")
colnames(coldata_ctrl) <-  c("sample_id", "condition")

coldata <- rbind.data.frame(coldata_CD,coldata_UC,coldata_ctrl) 

coldata$condition <- as.factor(coldata$condition)

#Obtenemos objeto de DESeq
dds=DESeqDataSetFromMatrix(countData = otus + 1, 
                           colData = coldata, design = ~ condition) 

#Normalizamos con Variance Stabilizing transformation
vst_matrix <- assay(varianceStabilizingTransformation(dds))

#------------------------------------------------------------------------------
# 3. Abundancia diferencial

DA <- DESeq(dds)

CD_vs_ctrl <- results(DA, alpha = 0.01, contrast = c("condition","CD","NonIBD")) 
UC_vs_ctrl <- results(DA, alpha = 0.01, contrast = c("condition","UC","NonIBD"))
CD_vs_UC <- results(DA, alpha = 0.01, contrast = c("condition","UC","CD"))

summary(CD_vs_ctrl)
summary(UC_vs_ctrl)
summary(CD_vs_UC)

as.data.frame(CD_vs_ctrl) %>% filter(., log2FoldChange >1)


#------------------------------------------------------------------------------
# 4. redes de co-abundancia o co-ocurrencia

#Queremos separar las matrices que correspondan a CD, UC y nonIBD

CD_otus <- vst_matrix[,colnames(vst_matrix) %in% CD_samples_ids]
UC_otus <- vst_matrix[,colnames(vst_matrix) %in% UC_samples_ids]
Ctrl_otus <- vst_matrix[,colnames(vst_matrix) %in% Ctrl_samples_ids]

#Calculamos matriz de correlación y matriz de significancia estadística para cada caso
CD_corr <- rcorr(t(CD_otus), type = "pearson")$r
CD_pval <- rcorr(t(CD_otus), type = "pearson")$P

UC_corr <- rcorr(t(UC_otus), type = "pearson")$r
UC_pval <- rcorr(t(UC_otus), type = "pearson")$P

Ctrl_corr <- rcorr(t(Ctrl_otus), type = "pearson")$r
Ctrl_pval <- rcorr(t(Ctrl_otus), type = "pearson")$P

#Creamos función para aplanar las matrices de correlación y significancia en una lista de interacciones
####################################################
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
###################################################
#Necesitamos corregir valores de p para pruebas múltiples 

CD_flat <- flattenCorrMatrix(CD_corr,CD_pval) %>% mutate(p_adjust = p.adjust(p)) %>% filter(., p_adjust <= 0.05 )
UC_flat <- flattenCorrMatrix(UC_corr,UC_pval) %>% mutate(p_adjust = p.adjust(p)) %>% filter(., p_adjust <= 0.05 )
Ctrl_flat <- flattenCorrMatrix(Ctrl_corr,Ctrl_pval) %>% mutate(p_adjust = p.adjust(p)) %>% filter(., p_adjust <= 0.05 )

CD_flat

fwrite(CD_flat, "CD_genus.tsv", sep = "\t")
fwrite(UC_flat, "UC_genus.tsv", sep = "\t")
fwrite(Ctrl_flat, "Ctrl_genus.tsv", sep = "\t")

tax_data$id  <- rownames(tax_data) 
fwrite(tax_data, "tax_data.tsv", sep = "\t")

#Creamos redes como objetos de igraph
g_CD <- graph_from_data_frame(CD_flat[,1:3], directed = F)
g_UC <- graph_from_data_frame(UC_flat[,1:3], directed = F)
g_Ctrl <- graph_from_data_frame(Ctrl_flat[,1:3], directed = F)

g_CD

V(g_CD)$label <- NA
V(g_UC)$label <- NA
V(g_Ctrl)$label <- NA

V(g_CD)$size <- 3
V(g_UC)$size <- 3
V(g_Ctrl)$size <- 3

plot(g_CD)
plot(g_UC)
plot(g_Ctrl)

#------------------------------------------------------------------------------
#5.Análisis de redes

#Obtenemos distribución del grado de cada red

#Calcular grado de cada nodo en cada red
V(g_CD)$degree <- degree(g_CD)
V(g_UC)$degree <- degree(g_UC)
V(g_Ctrl)$degree <- degree(g_Ctrl)


#Calculamos frecuencia de grado
 deg_freq_CD <- table(V(g_CD)$degree) %>% as.data.frame() %>% mutate(Var1 = as.numeric(Var1))
 deg_freq_UC <- table(V(g_UC)$degree) %>% as.data.frame() %>% mutate(Var1 = as.numeric(Var1))
 deg_freq_Ctrl <- table(V(g_Ctrl)$degree) %>% as.data.frame() %>% mutate(Var1 = as.numeric(Var1))
 
 
#Visualizamos distribuciones de grado log-log

ggplot(deg_freq_CD) +
 aes(x = Var1, y = Freq) +
 geom_point(shape = "bullet", size = 3.6, colour = "#FF07D4") +
 scale_x_continuous(trans = "log") +
 scale_y_continuous(trans = "log") +
 labs(x = "Grado", y = "Frecuencia del grado") +
 theme_minimal() +
  geom_smooth(method = "loess", se=F)+
 theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"))

ggplot(deg_freq_UC) +
  aes(x = Var1, y = Freq) +
  geom_point(shape = "bullet", size = 3.6, colour = "blue") +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(x = "Grado", y = "Frecuencia del grado") +
  theme_minimal() +
  geom_smooth(method = "loess", se=F)+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"))

ggplot(deg_freq_Ctrl) +
  aes(x = Var1, y = Freq) +
  geom_point(shape = "bullet", size = 3.6, colour = "tomato") +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(x = "Grado", y = "Frecuencia del grado") +
  theme_minimal() +
  geom_smooth(method = "loess", se=F)+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"))




#Distribuciones sin log-log

ggplot(deg_freq_CD) +
 aes(x = Var1, y = Freq) +
 geom_point(shape = "bullet", size = 4L, colour = "#D601FF") +
 geom_smooth(span = 1L) +
 labs(x = "Grado", y = "Frecuencia de Grado") +
 theme_minimal() +
 theme(axis.title.y = element_text(face = "bold"), 
 axis.title.x = element_text(face = "bold"))

ggplot(deg_freq_UC) +
  aes(x = Var1, y = Freq) +
  geom_point(shape = "bullet", size = 4L, colour = "blue") +
  geom_smooth(span = 1L) +
  labs(x = "Grado", y = "Frecuencia de Grado") +
  theme_minimal() +
  theme(axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"))


ggplot(deg_freq_Ctrl) +
  aes(x = Var1, y = Freq) +
  geom_point(shape = "bullet", size = 4L, colour = "tomato") +
  geom_smooth(span = 1L) +
  labs(x = "Grado", y = "Frecuencia de Grado") +
  theme_minimal() +
  theme(axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"))



 
 
 
















