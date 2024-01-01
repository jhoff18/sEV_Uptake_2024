library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
library(enrichR)
theme_set(theme_classic())
theme_update(text = element_text(size = 16))

# Proteomics --------------------------------------------------------------

### human data (CPC, MSC) ###
# grab abundance (PSMs and normalized abundances)
human = read.csv("./data/mdavis_sruti_24_01to12_filebased_multiconcensus.csv")

# get # PSMs for each protein across samples
humantidyPSMs = human[,c(4:5,27:38)] %>% 
  pivot_longer(starts_with('X')) %>% 
  mutate(symbol = str_match(Description, "GN=\\s*(.*?)\\s*PE=")[,2]) %>% 
  group_by(Accession, symbol) %>% 
  summarize(min_PSM = min(value,na.rm=T), max_PSM = max(value,na.rm=T), mean_PSM = mean(value,na.rm=T))

# set cut-off at >= 2 PSM in at least one sample
meets_cutoff_human = humantidyPSMs %>% filter(max_PSM > 1) %>% drop_na(symbol)

humanab = human %>% select(Accession, Description, contains("Abundances.."))
# assign sample names to abundances
colnames(humanab)[3:14] = c("CPC.A.Normoxia", "CPC.A.Hypoxia", 
                            "CPC.B.Normoxia", "CPC.B.Hypoxia", 
                            "CPC.C.Normoxia", "CPC.C.Hypoxia",
                            "MSC.D.Normoxia","MSC.D.Hypoxia",
                            "MSC.E.Normoxia","MSC.E.Hypoxia",
                            "MSC.F.Normoxia","MSC.F.Hypoxia")
humanab$symbol = str_match(humanab$Description, "GN=\\s*(.*?)\\s*PE=")[,2]

### rat data (CEC, CF) ###
rat = read.csv("./data/mdavis_sruti_24_13to24_filebased_multiconcensus.csv")
rattidyPSMs = rat[,c(4:5,27:38)] %>% 
  pivot_longer(starts_with('X')) %>% 
  mutate(symbol = toupper(str_match(Description, "GN=\\s*(.*?)\\s*PE=")[,2])) %>% 
  group_by(Accession, symbol) %>% 
  summarize(min_PSM = min(value,na.rm=T), max_PSM = max(value,na.rm=T), mean_PSM = mean(value,na.rm=T))

ratab = rat %>% select(Accession, Description, contains("Abundances.."))
colnames(ratab)[3:14] = c("CEC.D.Normoxia", "CEC.D.Hypoxia", 
                          "CEC.E.Normoxia", "CEC.E.Hypoxia", 
                          "CEC.F.Normoxia", "CEC.F.Hypoxia",
                          "RCF.A.Normoxia","RCF.A.Hypoxia",
                          "RCF.B.Normoxia","RCF.B.Hypoxia",
                          "RCF.C.Normoxia","RCF.C.Hypoxia")
ratab$symbol = str_match(ratab$Description, "GN=\\s*(.*?)\\s*PE=")[,2]

# set cut-off at >= 2 PSM in at least one sample
meets_cutoff_rat = rattidyPSMs %>% filter(max_PSM > 1) %>% drop_na(symbol)

# Venn diagram before

forvenn_human = humanab %>% 
  dplyr::select(!c(Description, Accession)) %>%
  pivot_longer(!symbol) %>% drop_na()

forvenn_rat = ratab %>% 
  dplyr::select(!c(Description, Accession)) %>% 
  pivot_longer(!symbol) %>% drop_na()

venn.diagram(list(CEC = (forvenn_rat %>% filter(str_sub(name, 1,3) == "CEC"))$symbol,
                  RCF = (forvenn_rat %>% filter(str_sub(name, 1,3) == "RCF"))$symbol,
                  CPC = (forvenn_human %>% filter(str_sub(name, 1,3) == "CPC"))$symbol,
                  MSC = (forvenn_human %>% filter(str_sub(name, 1,3) == "MSC"))$symbol),
             category.names = c("CEC", "CF", "CPC", "MSC"),
             filename = paste0("./plots/protein/venn_proteins_beforefilter", Sys.Date(),".png"),
             disable.logging = T,
             imagetype="png" ,
             height = 480 , 
             width = 480 , 
             resolution = 300,
             cat.cex = 0.7,
             cex = 0.5,
             fontfamily = "sans",
             cat.fontfamily = "sans",
             #print.mode=c("raw","percent"),
             col=c("#95C1AC", '#E389B6', '#F3C18A', '#A8A8CD'),
             fill = c(alpha("#95C1AC",0.3), alpha('#E389B6',0.3), alpha('#F3C18A',0.3), alpha('#A8A8CD',0.3)))


# combine all proteomics data ---------------------------------------------

##### AT LEAST 2 PSM IN A HUMAN AND RAT SAMPLE ##### 
intersect(meets_cutoff_rat$symbol, meets_cutoff_human$symbol)

combo = rbind(humanab %>% pivot_longer(cols = ends_with('oxia'), names_to = 'sample'),
              ratab %>% pivot_longer(cols = ends_with('oxia'), names_to = 'sample') %>% 
                mutate(symbol = toupper(symbol))) %>% 
  dplyr::select(!c(Accession, Description)) %>% 
  pivot_wider(names_from = sample, values_from = value, values_fn = function(x)mean(x,na.rm=T)) %>% 
  filter(symbol %in% intersect(meets_cutoff_rat$symbol, meets_cutoff_human$symbol)) %>% # PSM filter
  filter(str_sub(symbol,1,3) != "KRT") %>% # remove keratins
  column_to_rownames('symbol')

dim(combo)

# VENN DIAGRAM AFTER PSM FILTER

forvenn = combo %>% rownames_to_column('symbol') %>% 
  pivot_longer(cols = !symbol) %>% drop_na()

venn.diagram(list(CEC = (forvenn %>% filter(str_sub(name, 1,3) == "CEC"))$symbol,
                  RCF = (forvenn %>% filter(str_sub(name, 1,3) == "RCF"))$symbol,
                  CPC = (forvenn %>% filter(str_sub(name, 1,3) == "CPC"))$symbol,
                  MSC = (forvenn %>% filter(str_sub(name, 1,3) == "MSC"))$symbol),
             category.names = c("CEC", "CF", "CPC", "MSC"),
             filename = paste0("./plots/protein/venn_proteins_afterfilter", Sys.Date(),".png"),
             disable.logging = T,
             imagetype="png" ,
             height = 480 , 
             width = 480 , 
             resolution = 300,
             cat.cex = 0.7,
             cex = 0.5,
             fontfamily = "sans",
             cat.fontfamily = "sans",
             #print.mode=c("raw","percent"),
             col=c("#95C1AC", '#E389B6', '#F3C18A', '#A8A8CD'),
             fill = c(alpha("#95C1AC",0.3), alpha('#E389B6',0.3), alpha('#F3C18A',0.3), alpha('#A8A8CD',0.3)))

# keep proteins with less than 12 NAs
combo = combo[rowSums(is.na(combo))<12,]
dim(combo)

combo_wide = t(combo)

#replace NA with ⅕ minimum value for each sample
for(i in 1:ncol(combo_wide)){
  colum = combo_wide[,i]
  colum[is.na(colum)] = min(colum, na.rm = TRUE)/5
  combo_wide[,i] = colum}


# PCA ---------------------------------------------------------------------

PCA = prcomp(log10(combo_wide), scale = TRUE)
percent_var = round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
data_pca = data.frame(PC1 = PCA$x[,1], 
                    PC2 = PCA$x[,2],
                    sample = rownames(combo_wide),
                    Oxygen = str_sub(rownames(combo_wide),7,-1),
                    Cell = str_sub(rownames(combo_wide),1,3)) %>% 
  mutate(Cell = ifelse(Cell == "RCF", "CF", Cell))

proteinpca = ggplot(data_pca, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell, shape = Oxygen), size = 3) +
  ggtitle("Proteomics") +
  xlab(paste0("PC1, VarExp: ", percent_var[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percent_var[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16)) + 
  scale_color_manual(values = c(CPC = "#F3C18A", MSC = "#A8A8CD", CEC = "#95C1AC", CF = "#E389B6"))+
  geom_hline(yintercept=0, color = "grey", size=0.5) + 
  geom_vline(xintercept=0, color = "grey", size=0.5) +
  scale_x_continuous(breaks = seq(-90,90,5)) +
  scale_y_continuous(breaks = seq(-90,90,5)) 
proteinpca
ggsave(proteinpca, file = paste0("./plots/protein/proteinPCA_",Sys.Date(),".png"), width = 5, height = 4)


# EV markers -----------------------------------------------------------------

# Get boxplots for EV MARKERS (MISEV 2018)
markers = combo %>% as.data.frame() %>% #dplyr::select(CD63, CD81, HSP90AB1, HSPA8, APOB, ALB) %>% 
  rownames_to_column('Protein') %>% 
  pivot_longer(cols = !Protein) %>% 
  group_by(name) %>% 
  mutate(sum = sum(value, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(med = median(sum, na.rm = T),
         normfactor = med/sum, 
         value_norm = value*normfactor,
         value_norm_log = log10(value_norm)) %>% 
  mutate(cell = str_sub(name,1,3), 
         cell = ifelse(cell == "RCF", "CF", cell),
         oxygen = str_sub(name,7,-1))

markerpl = markers %>% filter(Protein %in% c('CD63', 'CD81', 'HSP90AB1', 'HSPA8', 'RPL14', 'RPS13')) %>% 
  ggplot(aes(x = str_sub(oxygen,1,3), y = log10(value), color = cell)) + 
  geom_hline(yintercept = c(7), linetype = "dotted") +
  facet_wrap(~Protein, ncol = 6) +
  geom_boxplot() + 
  theme(legend.position = "top") +
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "", y = "log10(intensity)", color = "") + 
  scale_color_manual(values = c(CPC = "#F3C18A", MSC = "#A8A8CD", CEC = "#95C1AC", CF = "#E389B6"))
markerpl
ggsave(markerpl, file =  paste0("./plots/protein/marker_boxplot_",Sys.Date(),".png"), width = 9, height = 3)

# protein DEA for normoxic samples -------------------------------------------------------------

# Get long data frame - normalized - for DEA
combo_long = combo_wide %>% t() %>% as.data.frame() %>% rownames_to_column('Protein') %>% 
  pivot_longer(cols = !Protein) %>% 
  group_by(name) %>% 
  mutate(sum = sum(value)) %>% 
  ungroup() %>% 
  mutate(med = median(sum),
         normfactor = med/sum, 
         value_norm = value*normfactor,
         value_norm_log = log10(value_norm)) %>%
  mutate(cell = str_sub(name,1,3),
         oxygen = str_sub(name,7,-1))

getDEA = function(condition1, condition2, condition1color, condition2color, fccutoff, pcutoff){
  res = combo_long %>% dplyr::filter(oxygen == "Normoxia" & (cell == condition1 | cell == condition2)) %>% 
    group_by(Protein, cell) %>% 
    mutate(mean = mean(value_norm,na.rm = T),
           sd = sd(value_norm,na.rm = T), 
           iqr = IQR(value_norm, na.rm = T),
           upperbound = mean + 1.5*iqr,
           lowerbound = mean - 1.5*iqr) %>% 
    dplyr::filter(value_norm < upperbound & value_norm > lowerbound) %>% 
    summarise(value = list(value_norm)) %>% 
    spread(cell, value) %>% 
    group_by(Protein)
  
  colnames(res)[2:3] = c("Condtion1", "Condition2")
  res = res %>% 
    mutate(p_value = t.test(unlist(Condtion1), unlist(Condition2), paired = F, alternative = 'two.sided', var.equal = T)$p.value,
           t_value = t.test(unlist(Condtion1), unlist(Condition2), paired = F, alternative = 'two.sided', var.equal = T)$statistic,
           mean_C1 = mean(unlist(Condtion1),na.rm = T),
           mean_C2 = mean(unlist(Condition2),na.rm = T),
           log2FC = log2(mean_C2/mean_C1),
           log2FC = ifelse(log2FC > 8, 8, log2FC), # set bounds on fold-change
           log2FC = ifelse(log2FC < -8, -8, log2FC),
           siggroup = ifelse(log2FC > fccutoff & p_value < pcutoff, 'condition2', "NS"),
           siggroup = ifelse(log2FC < -1*fccutoff & p_value < pcutoff, 'condition1', siggroup),
           siggroup = ifelse(abs(log2FC) > fccutoff & p_value > pcutoff, "*p only", siggroup)) %>% 
    arrange(p_value)
  colnames(res)[2:3] = c(condition1, condition2)
  colnames(res)[6:7] = c(paste0("mean_",condition1), paste0("mean_",condition2))
  
  labeled_points = rbind(res %>% filter(!(siggroup %in% c('NS','*p only'))) %>% ungroup() %>% slice_min(n=5, order_by = log2FC),
                         res %>% filter(!(siggroup %in% c('NS','*p only'))) %>% ungroup() %>% slice_max(n=5, order_by = log2FC))
  
  plotvp = ggplot(res, aes(x = log2FC, y = -log10(p_value), color = siggroup)) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = 'grey', size = 1) +
    geom_vline(xintercept = c(-1,1), linetype = "dotted", color = 'grey', size = 1) +
    geom_point() + 
    scale_color_manual(values = c('condition1' = condition1color, 'condition2' = condition2color,'NS' = 'grey', '*p only' = 'black'))+
    theme(text = element_text(size = 18),
          legend.position = "none") + 
    labs(x = "log2(fold-change)",
         y = "-log10(p-value)",
         title = paste0(condition1, " vs. ", condition2)) +
    ggrepel::geom_text_repel(data = labeled_points, aes(label = Protein), color = "grey30") +
    scale_x_continuous(breaks = seq(-10,10,2), limits = c(-8.1, 8.1))
  
  return(list(df = res, pl = plotvp))
}

getDEA('CPC', 'MSC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$pl
getDEA('CEC', 'RCF', "#95C1AC", "#E389B6", log2(2), 0.05)$pl

ggsave(getDEA('CPC', 'MSC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$pl,
       file =  paste0("./plots/protein/volcano_CPCvsMSC_",Sys.Date(),".png"), width = 4, height = 4)
ggsave(getDEA('CEC', 'RCF', "#95C1AC", "#E389B6", log2(2), 0.05)$pl,
       file =  paste0("./plots/protein/volcano_CECvsCF_",Sys.Date(),".png"), width = 4, height = 4)


# pathway -----------------------------------------------------------------

# grab proteins enriched in each cell type by pairwise comparisons (fold-change > 2)
# enrichR bug :( -- get identical results if you submit multiple sets at once --> have to run each individually

# color is just dummy variable here
up_in_CPC = rbind(getDEA('CPC', 'MSC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "MSC") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                  getDEA('CPC', 'RCF', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "RCF") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                  getDEA('CPC', 'CEC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "CEC") %>% dplyr::select(Protein, log2FC, p_value, comparison))

up_in_MSC = rbind(getDEA('MSC', 'CPC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "CPC") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                  getDEA('MSC', 'RCF', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "RCF") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                  getDEA('MSC', 'CEC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "CEC") %>% dplyr::select(Protein, log2FC, p_value, comparison))

up_in_CF = rbind(getDEA('RCF', 'MSC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "MSC") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                 getDEA('RCF', 'CPC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "CPC") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                 getDEA('RCF', 'CEC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "CEC") %>% dplyr::select(Protein, log2FC, p_value, comparison))

up_in_CEC = rbind(getDEA('CEC', 'MSC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "MSC") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                  getDEA('CEC', 'RCF', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "RCF") %>% dplyr::select(Protein, log2FC, p_value, comparison),
                  getDEA('CEC', 'CPC', "#F3C18A", "#A8A8CD", log2(2), 0.05)$df %>% filter(log2FC <= -1) %>% mutate(comparison = "CPC") %>% dplyr::select(Protein, log2FC, p_value, comparison))


pw_CPC = enrichR::enrichr(genes = up_in_CPC$Protein, databases = "GO_Biological_Process_2021")$GO_Biological_Process_2021 %>% 
  dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score) %>% 
  mutate(group = "CPC")

pw_CEC = enrichR::enrichr(genes = up_in_CEC$Protein, databases = "GO_Biological_Process_2021")$GO_Biological_Process_2021 %>% 
  dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score) %>% 
  mutate(group = "CEC")

pw_CF = enrichR::enrichr(genes = up_in_CF$Protein, databases = "GO_Biological_Process_2021")$GO_Biological_Process_2021 %>% 
  dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score) %>% 
  mutate(group = "CF")

pw_MSC = enrichR::enrichr(genes = up_in_MSC$Protein, databases = "GO_Biological_Process_2021")$GO_Biological_Process_2021 %>% 
  dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score) %>% 
  mutate(group = "MSC")

# Get top 15 GO:BP per cell group

pwpl = full_join(pw_MSC, full_join(pw_CPC, full_join(pw_CEC, pw_CF))) %>% 
  mutate(count = as.numeric(gsub("/.*$","",Overlap)),
         set_size = as.numeric(gsub(".*/","",Overlap)),
         pct = 100*(count/set_size)) %>% 
  group_by(group) %>% 
  slice_min(n=15, order_by = Adjusted.P.value) %>% 
  dplyr::select(Term, Adjusted.P.value, pct, group) %>% 
  ggplot(aes(x = group, y = Term, size = pct, color = -log10(Adjusted.P.value))) +
  geom_point() + 
  scale_size_binned(range = c(1, 10), limits = c(10,50)) +
  scale_color_viridis() + 
  # theme(legend.position = "top", legend.box="vertical") +
  labs(x = "", y = "", color = "-log10(padj)", size = "% of Proteins \n in GO pathway") 

pwpl
ggsave(pwpl, file =  paste0("./plots/protein/proteinDEA_pathway_",Sys.Date(),".png"), width = 13, height = 6)

# protein box plot --------------------------------------------------------

rep_df = data.frame(letter = toupper(letters[1:6]), rep = rep(1:3,2))

proteinbox = ggplot(combo %>% 
                      pivot_longer(cols = ends_with('oxia')) %>% 
                      mutate(oxygen = str_sub(name,7,-1),
                             cell = str_sub(name, 1, 3),
                             letter = str_sub(name, 5, 5),
                             cell = ifelse(cell == "RCF", "CF", cell)) %>% 
                      left_join(rep_df),
                    aes(x = rep, y = log10(value), color = oxygen, group  = interaction(oxygen, rep))) +
  geom_boxplot(aes(x = rep), outlier.shape = NA)+
  geom_point(position = position_jitterdodge(), size = 0.25, alpha = 0.3) +
  theme(legend.position = "top") +
  facet_wrap(~cell, ncol = 1) +
  labs(x = "Replicate", y = "log10(intensity)", color = "") +
  scale_color_manual(values = c("Normoxia" = "#EE857D", "Hypoxia" = "#6DAAF8"))
proteinbox
ggsave(proteinbox, file = paste0("./plots/protein/proteinbox_12", Sys.Date(),".png"), height = 8, width = 3)
