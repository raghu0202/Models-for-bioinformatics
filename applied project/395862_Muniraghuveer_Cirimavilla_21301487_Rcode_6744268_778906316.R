getwd()
setwd("C:\\Users\\raghu\\OneDrive\\Desktop\\MB assignment\\applied project")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.17")
#BiocManager::install("genefilter")
#BiocManager::install("limma")
#BiocManager::install("qvalue")
library(genefilter)
library(limma)
library(qvalue)



#install.packages(c("cluster", "gplots", "RColorBrewer", "survival"))
library(cluster)
library(gplots)
library(RColorBrewer)
library(survival)


#BiocManager::install("edgeR")
install.packages("edgeR", force = TRUE)

library(edgeR)
getwd()
setwd("C:\\Users\\raghu\\OneDrive\\Desktop\\MB assignment\\applied project")

Bc_clinical <- readRDS("bc_clinical.rds")
Bc_expression <- readRDS("BC_data_with_gene_info.rds")

#cluster analysis
dim(Bc_expression)
head(Bc_expression)
gene_data <- data.matrix(Bc_expression[, -(1:3)])
head(gene_data)


#scaling
scale(gene_data,center = TRUE)
scaled.gene_data<-t(scale(t(gene_data),center = TRUE))
dim(scaled.gene_data)

# Normalizing between arrays
normalise_Gene <- normalizeBetweenArrays(gene_data, method="scale")
par(mfrow=c(1,2))





#complete linkage cluster analysis
d <- dist(t(scaled.gene_data))
hc <- hclust(d, method="complete")
plot(hc, cex=0.5)
abline(h=220, col=2)



#selecting the optimal no of clusters
K<-2:5
sh_complete<-NULL
for (i in K) {
  sh_complete<- c(sh_complete, median(silhouette(cutree(hc, k = i), dist = d)[, 3], na.rm = T))
}
plot(K,sh_complete,type = "l", main = "Median silhouette (complete linkage)", xlab = "Number of clusters")
cl <- cutree(hc, k = K[which.max(sh_complete)])
cl

install.packages("BiocManager")
BiocManager::install("genefilter")
library(genefilter)
install.packages("gplots")
library(gplots)

#heat map
rv <- rowVars(scaled.gene_data)
idx<-order(-rv)[1:1200]
cols <- colors()[seq(8, length(colors()), len = length(unique(cl)))]
head(cbind(colnames(E), cols), n = 3)
heatmap.2(scaled.gene_data[idx, ], labCol = cl, trace = "none", ColSideColors = cols[cl], col = redgreen(100))


#Principal component analysis
par(bg = "white")
pc <- princomp(scaled.gene_data[idx, ])
plot(pc$load[, 1:2], col = cl)
title("PCs 1 and 2 of Breast cancer data, coloured by clusters")

library(qvalue)

#gene expression analysis
design <- model.matrix(~as.factor(cl))
DE.object <- lmFit(gene_data, design)
DE.object <- eBayes(DE.object)
qval.complete<-qvalue(DE.object$p.value[, 2], fdr.level = 0.05)
sum(qval.complete$significant)

#gene ontology test
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GO.db")
library(org.Hs.eg.db)
library(GO.db)
BiocManager::install("clusterProfiler")
library(clusterProfiler)

homo_sapiens <- org.Hs.eg.db 
my_symbols <- as.character(sig_genes) 
conversion_tab <- AnnotationDbi::select(homo_sapiens, 
                                        keys = my_symbols,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")



DE_genes <- conversion_tab$ENTREZID[conversion_tab$SYMBOL %in% rownames(DE.object)[which(qval.complete$significant)]]
library(limma)
GO_terms_present  <- goana(DE_genes, species='Hs')
GO_top_terms <-  topGO(GO_terms_present, n=200)

library(knitr)
GO_BP_terms <- GO_top_terms[GO_top_terms$Ont=='BP',][1:15,]
knitr::kable(GO_BP_terms)




install.packages("survival")
library(survival)

#survival anlysis
Bc_clinical<- readRDS("bc_clinical.rds")
gene.score<- colSums(gene_data[qval.complete$sig, ])
gene.score<- scale(gene.score)
cox.model <- coxph(Surv(Bc_clinical$Surv_time, Bc_clinical$event) ~ gene.score, data = Bc_clinical)
summary(cox.model)

fit_surv <- survfit(cox.model)
plot(fit_surv, main="Survival Curves", xlab="Time", ylab="Survival Probability")





head(Bc_clinical$Surv_Time)
colnames(Bc_clinical)



#Kaplan-Meier Survival Curves
Y <- Surv(Bc_clinical$Surv_time, Bc_clinical$event == TRUE)
kmfit <- survfit(Y ~Bc_clinical$age)
plot(kmfit, lty = c("solid", "dashed"), col = c("blue", "orange"),
     xlab = "Survival Time In Days", ylab = "Survival Probabilities")

Y <- Surv(Bc_clinical$Surv_time, Bc_clinical$event == TRUE)
kmfit <- survfit(Y ~Bc_clinical$tumor_size_mm)
plot(kmfit, lty = c("solid", "dashed"), col = c("green", "red"),
     xlab = "Survival Time In Days", ylab = "Survival Probabilities")


Y <- Surv(Bc_clinical$Surv_time, Bc_clinical$event == TRUE)
kmfit <- survfit(Y ~Bc_clinical$LNstatus)
plot(kmfit, lty = c("solid", "dashed"), col = c("red", "black"),
     xlab = "Survival Time In Days", ylab = "Survival Probabilities")


install.packages(c("ggplot2", "ggfortify"))
library(ggplot2)
library(ggfortify)


autoplot(kmfit) + 
  labs(x = "\n Survival Time (Years) since diagnosis ", y = "Survival Probabilities \n", 
       title = "K-M Survival Curves for the Breast Cancer Patients \n")


install.packages(c("plotly", "survminer"))
library(plotly)
library(survminer)

model_fit <- survfit(Surv(time = Surv_time, event = event) ~ LNstatus, data =Bc_clinical)
p1 <- ggsurvplot(model_fit)
plotly::ggplotly(p1[[1]])%>%
  layout(title = "\n Kaplan-Meier Survival Curves for the Breast Cancer Patients",
         xaxis = list(title = "Survival Time (Years) since diagnosis"),
         legend=list(title=list(text='LNstatus')))


model_fit <- survfit(Surv(time = Surv_time, event = event) ~ age, data =Bc_clinical)
p1 <- ggsurvplot(model_fit)
plotly::ggplotly(p1[[1]])%>%
  layout(title = "\n Kaplan-Meier Survival Curves for the Breast Cancer Patients",
         xaxis = list(title = "Survival Time (Years) since diagnosis"),
         legend=list(title=list(text='age')))


model_fit <- survfit(Surv(time = Surv_time, event = event) ~ tumor_size_mm, data =Bc_clinical)
p1 <- ggsurvplot(model_fit)
plotly::ggplotly(p1[[1]])%>%
  layout(title = "\n Kaplan-Meier Survival Curves for the Breast Cancer Patients",
         xaxis = list(title = "Survival Time (Years) since diagnosis"),
         legend=list(title=list(text='tumor_size_mm')))










