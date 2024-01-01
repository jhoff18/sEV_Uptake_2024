library(tidyverse)
library(VennDiagram)
library(mdatools)
library(pheatmap)
library(RColorBrewer)
library(enrichR)
theme_set(theme_classic())
theme_update(text = element_text(size = 16))

# read in uptake values --> similarity heatmap -----------------------------------------------------------

up = read.csv("./data/uptakevals.csv") %>% 
  dplyr::rename(Sample = X)
upnorm = up %>% filter(grepl("Normoxia",Sample)) %>% 
  pivot_longer(cols = !Sample) %>% 
  mutate(name = str_sub(name,3,-1),
         name = gsub("[.]", "", name),
         cell = str_sub(name,-3,-1),
         cell = ifelse(cell == "RCF", "CF", cell),
         method = str_sub(name,1,-4))

matrix_uptake = upnorm %>% 
  dplyr::select(Sample, name, value) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  column_to_rownames('Sample')

annotc = data.frame(name = colnames(matrix_uptake)) %>% 
  left_join(upnorm %>% dplyr:::select(name, cell, method) %>% distinct()) %>% 
  column_to_rownames('name') %>% 
  mutate(method = ifelse(method == "CaveolaeLipidRaft", "Caveolae/Lipid Raft", method))

annotation_colors = list(cell = c("CEC" = "#9A7787", CF = "#E4AFB0"),
                         method = c("Clathrin" = "#FFCE3F", "Macropinocytosis" = "#06304B", 
                                    "Caveolae/Lipid Raft" = "#0094B6"))

pheatmap(as.matrix(dist(scale(t(matrix_uptake)),method = "euclidean")),
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(16),
         annotation_col = annotc,
         annotation_row = annotc,
         annotation_colors = annotation_colors,
         #width = 6.5, height = 5, filename = "outcomedistancehm_111122.png"
         show_rownames = F,
         show_colnames = F)


# read in lipid and protein data ------------------------------------------

lipid = read.csv("./data/clean_lipid_features.csv")
protein = read.csv("./data/clean_protein_features.csv") 

# PLSR -----------------------------------------------------------------------

# CENTER AND SCALE PROTEINS AND LIPIDS
scaled_lipid = prep.autoscale(lipid %>% 
                                column_to_rownames('Name') %>% 
                                log10() %>% 
                                dplyr::select(ends_with('Normoxia')) %>% 
                                t(), 
                              center = TRUE, scale = TRUE)
rownames(scaled_lipid) = gsub("[.]", "-", rownames(scaled_lipid))
scaled_lipid[1:5,1:5]


scaled_prot = prep.autoscale(protein %>% 
                               filter(grepl("Normoxia",X)) %>% 
                               column_to_rownames('X'),
                             center = TRUE, scale = TRUE)
rownames(scaled_prot) = gsub("[.]", "-", rownames(scaled_prot))
scaled_prot[1:5,1:5]

# combine protein and lipid --> X 
bothnormoxia = rbind(t(scaled_lipid), t(scaled_prot)) %>% 
  as.data.frame() %>% 
  dplyr::select(!`MSC-F-Normoxia`)

# uptake values --> Y
y_normoxia = data.frame(Sample = colnames(bothnormoxia)) %>% 
  left_join(up %>% mutate(Sample = gsub("[.]", "-", Sample))) %>% 
  column_to_rownames('Sample')
colnames(y_normoxia) = str_sub(gsub("[.]", "", colnames(y_normoxia)), 2,-1)


# plot outcomes before and after center & scale ---------------------------

y_grid = grid.arrange(ggplot(y_normoxia %>% 
                               rownames_to_column('sample') %>% 
                               pivot_longer(!sample) %>% 
                               mutate(recep_cell = str_sub(name,-3,-1),
                                      recep_cell = ifelse(recep_cell == "RCF", "CF", recep_cell),
                                      cell = str_sub(sample,1,3),
                                      cell = ifelse(cell == "RCF", "CF", cell),
                                      method = str_sub(name, 1,-4)), 
                             aes(x = paste0(recep_cell, "-", str_sub(method,1,4)), y = value)) + 
                        geom_boxplot(aes(color = method), outlier.shape = NA) +
                        geom_jitter(alpha = 0.8, aes(color = cell), size = 2) +
                        scale_color_manual(values = c("Clathrin" = "#FFCE3F", "Macropinocytosis" = "#06304B", "CaveolaeLipidRaft" = "#0094B6",
                                                      CPC = "#F3C18A", MSC = "#A8A8CD", CEC = "#95C1AC", CF = "#E389B6")) +
                        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                              legend.position = "none") + 
                        labs(x = "", subtitle = "Original"),
                      
                      ggplot(prep.autoscale(y_normoxia, center = TRUE, scale = TRUE) %>% 
                               as.data.frame() %>% 
                               rownames_to_column('sample') %>% 
                               pivot_longer(!sample) %>% 
                               mutate(recep_cell = str_sub(name,-3,-1),
                                      recep_cell = ifelse(recep_cell == "RCF", "CF", recep_cell),
                                      cell = str_sub(sample,1,3),
                                      cell = ifelse(cell == "RCF", "CF", cell),
                                      method = str_sub(name, 1,-4)),
                             aes(x = paste0(recep_cell, "-", str_sub(method,1,4)), y = value)) + 
                        geom_boxplot(aes(color = method), outlier.shape = NA) +
                        geom_jitter(alpha = 0.8, aes(color = cell), size = 2) +
                        scale_color_manual(values = c("Clathrin" = "#FFCE3F", "Macropinocytosis" = "#06304B", "CaveolaeLipidRaft" = "#0094B6",
                                                      CPC = "#F3C18A", MSC = "#A8A8CD", CEC = "#95C1AC", CF = "#E389B6")) +
                        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                              legend.position = "none") + 
                        labs(x = "", subtitle = "Centered & Scaled"), ncol = 2)

ggsave(y_grid, file = paste0("./plots/PLSR/outcomes_scaled_boxplot_",Sys.Date(),".png"), height = 3.5, width = 7)



# create PLSR model -------------------------------------------------------

# check matching names
table(rownames(y_normoxia) == colnames(bothnormoxia))

mod = pls(x = t(bothnormoxia), y = prep.autoscale(y_normoxia, scale = TRUE, center = TRUE), 10, cv = 1)

makereg = function(x,y){
  mod = pls(x = x, y = prep.autoscale(y, scale = TRUE, center = TRUE), 10, cv = 1)
  return(mod)}

plottraining = function(mod,ncomp,row,column){
  par(mfrow = c(row,column))
  for (i in 1:nrow(mod$yloadings)) {plotRMSE(mod, ny = i)}# each is #450x450
  dev.off()
  par(mfrow = c(row,column))
  for (i in 1:nrow(mod$yloadings)) {plotPredictions(mod, ny = i)}
  dev.off()
  plotYCumVariance(mod, type = "h", show.labels = TRUE)
  plotXCumVariance(mod, type = "h", show.labels = TRUE)
  plotXResiduals(mod, show.labels=TRUE, ncomp=ncomp)
  summary(mod, ncomp = ncomp)}

getVIP = function(mod,x,y,ncomp,threshold, row,column){
  vip = vipscores(mod, ncomp = ncomp)
  m3 = pls(x = prep.autoscale(x, scale = TRUE, center = TRUE), 
           y = prep.autoscale(y, scale = TRUE, center = TRUE), 10, cv = 1, exclcols = (rowMeans(vip) < threshold))
  paste0("# VIPs: ",table(rowMeans(vipscores(mod, ncomp=3))<1.5)[1], " - threshold = ", threshold)
  summary(m3, ncomp = ncomp)
  return(m3)}

# test! 
# CEC_uptake = makereg(t(bothnormoxia), y_normoxia %>% dplyr::select(ends_with('CEC')))
# summary(CEC_uptake, ncomp = 3)

uptake = makereg(t(bothnormoxia), prep.autoscale(y_normoxia, center = TRUE, scale = TRUE))
summary(uptake, ncomp = 3) # results from full model - 3 components
plotYCumVariance(uptake, type = "h", show.labels = TRUE)
plotXCumVariance(uptake, type = "h", show.labels = TRUE)
plotXScores(uptake,show.label = TRUE)
plotXYLoadings(uptake)


# make pretty ggplot PLSR plots -------------------------------------------

scoreplot = function(mod, titl){
  scores = as.data.frame(mod$res$cal$xdecomp$scores)[,1:2] %>% rownames_to_column("Sample") %>% 
    mutate(Group = str_sub(Sample,1,3),
           Group = ifelse(Group == "RCF", "CF", Group),
           letter = str_sub(Sample,5,5)) %>% 
    dplyr::rename("Component 1" = "Comp 1",
                  "Component 2" = "Comp 2") %>% 
    left_join(rep_df) %>% 
    mutate(Sample = paste0(Group, "-", rep))
  expvar = c(mod$res$cal$xdecomp$expvar[[1]], mod$res$cal$xdecomp$expvar[[2]])
  
  pcaplot = ggplot(scores, aes(x = `Component 1`, y = `Component 2`, label = Sample)) + 
    geom_point(size = 4, aes(color = Group)) + 
    scale_color_manual(values=c(CPC = "#F3C18A", MSC = "#A8A8CD", CEC = "#95C1AC", CF = "#E389B6")) +
    theme(plot.title = element_text(hjust = 0.5), 
          text = element_text(size = 18), 
          legend.key = element_blank(),
          legend.position = 'top')+
    geom_hline(yintercept=0, color = "gray45", size=0.25) + 
    geom_vline(xintercept=0, color = "gray45", size=0.25) + 
    geom_text_repel(size = 4)+
    labs(title = titl,
         x = paste0("Component 1 (",round(expvar[1],2),"%)"),
         y = paste0("Component 2 (",round(expvar[2],2),"%)"),
         color = "")
  return(pcaplot)
}

plotloadings = function(mod, titl){
  yloadings = as.data.frame(mod$yloadings)[,1:2] %>% 
    rownames_to_column("Label") %>% 
    mutate(method = str_sub(Label, 1, -4),
           cell = str_sub(Label,-3,-1),
           cell = ifelse(cell == "RCF", "CF", cell),
           Type = "Uptake")
  
  xloadings = as.data.frame(mod$xloadings)[,1:2] %>% 
    rownames_to_column("Label") %>% 
    mutate(Type = ifelse(Label %in% colnames(scaled_lipid), "Lipid", "Protein"),
           method = NA,
           cell = NA)
  loadings = rbind(yloadings,xloadings) %>% 
    dplyr::rename("Component 1" = "Comp 1",
                  "Component 2" = "Comp 2")
  expvar = c(mod$res$cal$xdecomp$expvar[[1]], mod$res$cal$xdecomp$expvar[[2]])
  
  loadingspl = ggplot(loadings, aes(x = `Component 1`, y = `Component 2`, label = paste0(cell,"-",str_sub(Label,1,4)))) + 
    geom_point(data = loadings %>% filter(is.na(cell)), size = 1, aes(color = Type)) + 
    geom_point(data = loadings %>% filter(!is.na(cell)), size = 3, aes(color = Type)) + 
    scale_color_manual(values=c(Protein="grey", Lipid="black", `Uptake` ="red")) +
    theme(plot.title = element_text(hjust = 0.5), 
          text = element_text(size = 18), 
          legend.key = element_blank(),
          legend.position = 'top')+
    geom_hline(yintercept=0, color = "gray45", size=0.25) + 
    geom_vline(xintercept=0, color = "gray45", size=0.25) + 
    ggrepel::geom_text_repel(data = loadings %>% filter(!is.na(cell)), size = 4)+
    labs(title = titl,
         x = paste0("Component 1 (",round(expvar[1],2),"%)"),
         y = paste0("Component 2 (",round(expvar[2],2),"%)"), 
         color = "")
  return(loadingspl)
}


scoreplot(uptake, titl = "Full Model")
ggsave(scoreplot(uptake, titl = "Full Model"), file = paste0("./plots/PLSR/fullPLSRmodel_PCA_",Sys.Date(),".png"), width = 5, height = 5)

plotloadings(uptake, titl = "Full Model")
ggsave(plotloadings(uptake, titl = "Full Model"), file = paste0("./plots/PLSR/fullPLSRmodel_loadings_PCA_",Sys.Date(),".png"), width = 5, height = 5)


# get VIPs ----------------------------------------------------------------

nrow(bothnormoxia)
threshold = 1

vip = vipscores(uptake, ncomp = 3)
newuptake = pls(x = t(bothnormoxia), 
                y = prep.autoscale(y_normoxia, scale = TRUE, center = TRUE), 
                10, cv = 1, 
                exclcols = (rowMeans(vip) < threshold))
paste0("# VIPs: ",table(rowMeans(vipscores(uptake, ncomp=3))<threshold)[1], " - threshold = ", threshold)
summary(newuptake, ncomp = 2)

scoreplot(newuptake, titl = "Reduced Model")
ggsave(scoreplot(newuptake, titl = "Reduced Model"), file = paste0("./plots/PLSR/reducedPLSRmodel_PCA_",Sys.Date(),".png"), width = 5, height = 5)

plotloadings(newuptake, titl = "Reduced Model")
ggsave(plotloadings(newuptake, titl = "Reduced Model"), file = paste0("./plots/PLSR/reducedPLSRmodel_loadings_PCA_",Sys.Date(),".png"), width = 5, height = 5)


# plot top VIP heatmap ----------------------------------------------------

vip_df = data.frame(vip)
vip_df$m = rowMeans(vip_df)
sel2 = data.frame(vip_df) %>% arrange(-m) 
sel = sel2 %>% select(!m)

annotc = data.frame(sample = colnames(sel),
                    Cell = str_sub(colnames(sel),-3,-1),
                    Uptake = str_sub(colnames(sel),1,-4)) %>% 
  mutate(Cell = ifelse(Cell == "RCF", "CF", Cell),
         Uptake = ifelse(Uptake == "CaveolaeLipidRaft", "Caveolae/Lipid Raft", Uptake)) %>% 
  column_to_rownames('sample')

annotation_colors = list(Cell = c("CEC" = "#9A7787", CF = "#E4AFB0"),
                         Uptake = c("Clathrin" = "#FFCE3F", "Macropinocytosis" = "#06304B", 
                                    "Caveolae/Lipid Raft" = "#0094B6"),
                         Feature = c(Protein="grey", Lipid="black"))
annotr = sel2 %>% dplyr::select(m) %>% dplyr::rename(`mean(VIP)` = m) %>% 
  rownames_to_column('feat') %>% 
  mutate(Feature = ifelse(feat %in% colnames(scaled_lipid), "Lipid", "Protein")) %>% 
  column_to_rownames('feat')

pheatmap(t(sel[1:25,]), 
         annotation_col = annotr %>% head(25),
         annotation_row = annotc,
         annotation_colors = annotation_colors,
         show_rownames = F,
         fontsize = 9,
         color = colorRampPalette(c('white','#ffc9bb', '#FF6242'))(10))



# plsr predict hypoxic samples-----------------------------------------------------------

# CENTER AND SCALE PROTEINS AND LIPIDS
scaled_lipid = prep.autoscale(lipidmergeclean %>% 
                                column_to_rownames('Name') %>% 
                                log10() %>% 
                                dplyr::select(ends_with('Hypoxia')) %>% 
                                t(), 
                              center = TRUE, scale = TRUE)

scaled_prot = prep.autoscale(combo_wide %>% as.data.frame() %>% 
                               rownames_to_column('sample') %>% 
                               filter(grepl("Hypoxia",sample)) %>% 
                               column_to_rownames('sample'),
                             center = TRUE, scale = TRUE)
rownames(scaled_prot) = gsub("[.]", "-", rownames(scaled_prot))

bothhypoxia = rbind(t(scaled_lipid), t(scaled_prot)) %>% 
  as.data.frame() %>% 
  dplyr::select(!`MSC-F-Hypoxia`)

y_hypomoxia = data.frame(Sample = colnames(bothhypoxia)) %>% 
  left_join(up %>% mutate(Sample = gsub("[.]", "-", Sample))) %>% 
  column_to_rownames('Sample')
colnames(y_hypomoxia) = str_sub(gsub("[.]", "", colnames(y_hypomoxia)), 2,-1)

# check matching names
table(rownames(y_hypomoxia) == colnames(bothhypoxia))

predictions_hypoxia = predict(newuptake, t(bothhypoxia), prep.autoscale(y_hypomoxia, center = TRUE, scale = TRUE))

predicted = data.frame(predictions_hypoxia$y.pred) %>% rownames_to_column('sample') %>% 
  pivot_longer(!sample, values_to = 'predicted') %>% 
  mutate(Component = as.numeric(str_sub(name,6,6)),
         outcome = str_sub(name, 8,-1)) %>% 
  filter(Component <= 3)
observed = data.frame(predictions_hypoxia$y.ref) %>% rownames_to_column('sample') %>% 
  pivot_longer(!sample, values_to = 'observed', names_to = 'outcome')

pred_obs = full_join(predicted, observed) %>% filter(Component == 3) %>% 
  mutate(recipient = str_sub(outcome, -3, -1),
         recipient = ifelse(recipient == "RCF", "CF", recipient),
         method = str_sub(outcome, 1, -4), 
         method = ifelse(method == "CaveolaeLipidRaft","Caveolae/Lipid Raft",method))
# from MDA TOOLS
plotPredictions(predictions_hypoxia, ny = 6, ncomp = 3, show.stat = T, show.lines = c(0,0))

# GGPLOT OBS VS. PRED.
redpred = ggplot(pred_obs, aes(x = observed, y = predicted, color = method)) + 
  geom_hline(yintercept = 0, linewidth = 0.5, color = 'grey') + 
  geom_vline(xintercept = 0, linewidth = 0.5, color = 'grey') +
  geom_point(size = 3, alpha = 0.5) + 
  facet_wrap(~paste0(recipient, ": ", method)) + 
  stat_poly_eq(col = "black") + 
  stat_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + 
  theme(legend.position = "none") + 
  scale_color_manual(values = c("Clathrin" = "#FFCE3F", "Macropinocytosis" = "#06304B", "Caveolae/Lipid Raft" = "#0094B6")) + 
  labs(x = "Observed", y = "Predicted") + 
  scale_x_continuous(breaks = seq(-3,3,1), limits = c(-3,3)) +
  scale_y_continuous(breaks = seq(-3,3,1), limits = c(-2,2))

redpred
ggsave(redpred, file = paste0("./plots/PLSR/reducedPLSRmodel_predictions_",Sys.Date(),".png"), width = 6, height = 4)
