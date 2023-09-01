library(datapasta)
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)
library(plotly)
library(DT)
library(RColorBrewer)
library(heatmaply)

data_for_DGE <-as_tibble(read_tsv("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Karen_Cancer_Metabolomics\\GSE153713_hNKrawcount.txt"))%>% 
  column_to_rownames("...1")
myDGEList <- DGEList(data_for_DGE)
myDGEList<-myDGEList #%>% column_to_rownames("...1")
sampleLabels <- colnames(myDGEList)

data_for_DGE

log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")%>% column_to_rownames("geneID")
#colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = Healthy1:TB6, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

log2.cpm.df.pivot

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 95,
               size = 10,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=5 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = Healthy1:TB6, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 95,
               size = 10,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = Healthy1:TB6, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               shape = 95,
               size = 10,
               color = "black",
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

group <- c("Healthy","Healthy","Healthy","Healthy","Healthy","Healthy",
           "Patient","Patient","Patient","Patient","Patient","Patient",
           "Ascites1","Ascites1","Ascites1","Ascites1","Ascites1","Ascites1",
           "Ascites2","Ascites2","Ascites2","Ascites2","Ascites2","Ascites2")
group <- factor(group)
group

View(log2.cpm.filtered.norm)
#PCA ----
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA Plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
  
ggplotly(pca.plot)

#Gene Search ----
mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = (Healthy1 + Healthy2 + Healthy3 + Healthy4 + Healthy5 + Healthy6)/6, 
                    patient.AVG = (Patient1 + Patient2 + Patient3 + Patient4 + Patient5 + Patient6)/6,
                    TA.AVG = (TA1 + TA2 + TA3 + TA4 + TA5 + TA6)/6,
                    TB.AVG = (TB1 + TB2 + TB3 + TB4 + TB5 + TB6)/6,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (patient.AVG - healthy.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

datatable(mydata.df[,c(1,26:28)], 
          extensions = c('KeyTable', "FixedHeader"), 
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))


# Volcano ----
#group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = Patient - Healthy,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "NK cells Healthy vs. Patients",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

# Heatmaps ----
#v.DEGList.filtered.norm$E[rownames(v.DEGList.filtered.norm$E) %in% c("GPKOW" ,"SH3BGRL3" ),]
#rownames(v.DEGList.filtered.norm$E) %in% abc$`NK FUNCTION` )
metabolic_gene_lists<-read_csv("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Karen_Cancer_Metabolomics\\gene list.csv",locale=locale(encoding="latin1"))
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))

#metabolic_gene_lists
myplots <- vector('list', ncol(metabolic_gene_lists))

current_number =1 
(heatmaply(v.DEGList.filtered.norm$E[rownames(v.DEGList.filtered.norm$E) 
                                     %in% (metabolic_gene_lists[[colnames(metabolic_gene_lists)[current_number]]]) ,],dendrogram = "row",scale = "row",colors = myheatcolors, main = colnames(metabolic_gene_lists)[current_number]))


NKs <- c("GZMB", "PRF1", "LAMP1", "FCGR3A", "IFNG", "TNF", "IL2", "TNFSF10", "FASLG")
Lipid_Uptake <- c("CD36", "FABP", "FABP2", "FABP3", "FABP4", "FABP5", "FABP6", "FABP7", "FABP8", "FABP9", "FABP12", "LDLR", "FFAR1", "FFAR2", "FFAR3", "FFAR4")
Mito_FAO <- c("CPT1A", "CPT1B", "ACSl1", "SLC25A20", "ACOT13", "ACOX1", "HADHA", "HADHB", "ECHS1", "HADH", "ACADL", "ACADM", "ACADS", "ACADVL", "ACSM6", "MECR", "ACSM3", "MMAA", "MMUT", "PCCA", "PCCB", "MCEE", "ACAA2", "ACOT2", "ACOT7", "THEM4", "ACOT12", "DBI", "ECI1", "DECR1", "ACOT9", "ACOT11", "MCAT", "THEM5", "ACBD7", "NDUFAB1", "ACOT13", "PCTP", "ACOT1", "ACSF2", "ACAD10", "ACAD11", "ACBD6")

NKs_hm <- v.DEGList.filtered.norm$E[rownames(v.DEGList.filtered.norm$E) %in% NKs, ]

heatmaply(NKs_hm, dendrogram = "row",scale = "row",colors = myheatcolors, main = 'NK Genes')

group_labels = c("Healthy","Healthy","Healthy","Healthy","Healthy","Healthy",
                  "Patient","Patient","Patient","Patient","Patient","Patient",
                  "Ascites1","Ascites1","Ascites1","Ascites1","Ascites1","Ascites1",
                  "Ascites2","Ascites2","Ascites2","Ascites2","Ascites2","Ascites2")

v.DEGList.filtered.norm$E <- rbind(v.DEGList.filtered.norm$E, group_labels)
v.DEGList.filtered.norm$E
#write.csv(v.DEGList.filtered.norm$E, 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Karen_Cancer_Metabolomics\\normal_counts.csv', row.names = T)

#write.csv(mydata.df, 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Karen_Cancer_Metabolomics\\average counts.csv', row.names = T)
