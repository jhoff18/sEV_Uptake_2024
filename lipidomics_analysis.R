library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
theme_set(theme_classic())
theme_update(text = element_text(size = 16))

# read in data ------------------------------------------------------------

read_lipidomics = function(sheetnumber){
  columns = c(3,4,6,7,10:21)
  df = readxl::read_excel("./data/lipidomics_annotated.xlsx",sheetnumber)[,columns]
  colnames(df) = gsub("\\..*","",colnames(df))
  df = df %>% 
    dplyr::rename(RT = `RT [min]`, MW = `Molecular Weight`) %>% 
    mutate(LCMode = ifelse(LCMode == "RP_Pos", "RPPos", LCMode),
           name_mode = paste0(Name,"_",LCMode)) %>% 
    pivot_longer(cols = starts_with("QC"), names_to = "QC", values_to = "QC_value") %>% 
    group_by(across(c(-QC, -QC_value))) %>% 
    summarize(sd_QC = sd(QC_value, na.rm = T),
              mean_QC = mean(QC_value, na.rm = T),
              cv_QC = 100*(sd_QC/mean_QC)) %>% 
    dplyr::select(!c(sd_QC, mean_QC)) %>% 
    ungroup() %>% 
    group_by(name_mode) %>% mutate(n = n()) %>% 
    ungroup()
  
  # Take care of repeats and isomers
  # if RT is the same --> take the lowest CV value
  # if RT is different --> isomers --> take the peak sum
  
  df_split = df %>% filter(n != 1) %>% group_by(name_mode) %>% mutate(range = max(RT) - min(RT))
  
  lowestcv = df_split %>% filter(range <= 0.3) %>% arrange(name_mode, cv_QC) %>% slice_head(n=1)
  
  summed = df_split %>% 
    filter(range > 0.3) %>%  
    mutate(MW = round(MW,3)) %>% 
    pivot_longer(cols = ends_with("oxia"), names_to = 'sample') %>% 
    group_by(name_mode, sample) %>% 
    dplyr::select(!c(RT, cv_QC)) %>% 
    pivot_wider(names_from = sample, values_from = value, values_fn = function(x)sum(x,na.rm=T)) %>% 
    distinct()
  
  df = rbind(lowestcv, summed, df %>% filter(n == 1)) %>% 
    dplyr::select(!c(cv_QC, n, range)) %>% 
    pivot_longer(cols = ends_with("oxia"), names_to = 'sample') %>% 
    mutate(cell = str_split_i(sample, "-", 1),
           rep = str_split_i(sample, "-", 2),
           oxygen = str_split_i(sample, "-", 3))
  
  return(df)}

# example
read_lipidomics(4)

df_all = rbind(read_lipidomics(1),
               read_lipidomics(2),
               read_lipidomics(3),
               read_lipidomics(4))

# venn diagram of detected lipids ------------------------------------------------------------
venn.diagram(list(CEC = unique((df_all %>% filter(cell == "CEC"))$name_mode),
                  RCF = unique((df_all %>% filter(cell == "RCF"))$name_mode),
                  CPC = unique((df_all %>% filter(cell == "CPC"))$name_mode),
                  MSC = unique((df_all %>% filter(cell == "MSC"))$name_mode)),
             category.names = c("CEC", "CF", "CPC", "MSC"),
             filename = paste0("./plots/lipid/venn_lipids", Sys.Date(),".png"),
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

# filter and clean -----------------------------------------------------------------------

# Merge lipidomics
# 1. merge cell type sheets by lipid name+mode
# 2. Filter: keep lipids < 12 NAs
# 3. clean = no "similar to" unknowns
# 4. replace NA with ⅕ minimum value for each sample

lipidmergelong = df_all %>% group_by(LCMode, Name, name_mode) %>% 
  mutate(n = n_distinct(sample),
         n_NA = 24 - n) %>% 
  arrange(-n_NA) %>% 
  filter(n_NA < 12)

print(paste0("# features, w/ 'similar to' unknowns: ", length(unique(df_all$name_mode))))

lipidmerge = lipidmergelong %>% ungroup() %>% 
  dplyr::select(Name, name_mode, sample, value) %>% 
  pivot_wider(names_from = sample, values_from = value, values_fn = function(x)mean(x,na.rm = T))

#replace NA with ⅕ minimum value for each sample
for(i in 3:ncol(lipidmerge)){
  coll = lipidmerge[,i]
  coll[is.na(coll)] = min(coll, na.rm = TRUE)/5
  lipidmerge[,i] = coll}

lipidmergeclean = lipidmerge %>% filter(!grepl("Similar to",Name))
print(paste0("# features, w/o 'similar to' unknowns: ", nrow(lipidmergeclean)))

# Get common names
lipidlynx = readxl::read_excel("./data/LipidLynx_concertednames.xlsx") %>% 
  dplyr::rename(Name = TextInput, Standard_Name = TextInput_converted)

lipidmergeclean_long = lipidmergeclean %>% 
  pivot_longer(cols = ends_with('oxia'), names_to = 'sample') %>% 
  mutate(cell = str_split_i(sample, "-", 1),
         rep = str_split_i(sample, "-", 2),
         oxygen = str_split_i(sample, "-", 3),
         cell = ifelse(cell == "RCF", "CF", cell)) %>% 
  left_join(lipidlynx) %>% 
  mutate(Name = ifelse(!is.na(Standard_Name), Standard_Name, Name)) %>% 
  filter(!(Name %in% c('Penicillin G', 'Pyridoxine', 'Pantothenate', 'JWH 018 7-hydroxyindole metabolite')))

lipidmergeclean = lipidmergeclean_long %>% dplyr::select(Name, sample, value) %>% 
  pivot_wider(names_from = sample, values_from = value)

rep_df = data.frame(letter = toupper(letters[1:6]), rep = rep(1:3,2))

lipidbox = ggplot(lipidmergeclean_long %>% dplyr::rename(letter = rep) %>% left_join(rep_df),
                  aes(x = rep, y = log10(value), color = oxygen, group = interaction(oxygen, rep))) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(), size = 0.5, alpha = 0.3) +
  facet_wrap(~cell, ncol = 1) +
  theme(legend.position = "top") +
  labs(x = "Replicate", y = "log10(intensity)", color = "") +
  scale_color_manual(values = c("Normoxia" = "#EE857D", "Hypoxia" = "#6DAAF8"))
lipidbox
ggsave(lipidbox, file = paste0("./plots/lipid/lipidbox_12", Sys.Date(),".png"), height = 8, width = 3)

# heatmap -----------------------------------------------------------------

lipidlynx = readxl::read_excel("./data/LipidLynx_concertednames.xlsx") %>% 
  dplyr::rename(Name = TextInput, Standard_Name = TextInput_converted)

lipf = left_join(lipidmergeclean, 
                 lipidlynx, 
                 by = "Name") %>% 
  filter(!is.na(Standard_Name)) %>% 
  column_to_rownames('Standard_Name') %>% 
  dplyr::select(!Name)

lipf = log10(lipf)
head(lipf)

annotation_samples = data.frame(sample = colnames(lipf)) %>% 
  mutate(Cell = str_split_i(sample,"-",1),
         Cell = ifelse(Cell == "RCF", "CF", Cell),
         Oxygen = str_split_i(sample,"-",3)) %>% 
  column_to_rownames('sample')

annotation_lipids = data.frame(lipid = rownames(lipf)) %>% 
  mutate(Class = sub("\\s*\\(.*", "", lipid)) %>% 
  column_to_rownames('lipid')

annotation_colors = list(Oxygen = c("Normoxia" = "#EE857D", "Hypoxia" = "#6DAAF8"),
                         Cell = c(CPC = "#F3C18A", MSC = "#A8A8CD", CEC = "#95C1AC", CF = "#E389B6"))

pheatmap(lipf,
         annotation_col = annotation_samples,
         annotation_colors = annotation_colors,
         annotation_row = annotation_lipids,
         show_colnames = F,
         color = rev(brewer.pal(n = 11, name = "RdYlBu")))


# pca ---------------------------------------------------------------------

tlipid = log10(t(lipidmergeclean %>% select(!Name)))
PCA = prcomp(tlipid, scale = TRUE)
percent_var = round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
data_PCA = data.frame(PC1 = PCA$x[,1], 
                    PC2 = PCA$x[,2],
                    sample = colnames(lipidmergeclean)[2:ncol(lipidmergeclean)],
                    Oxygen = str_split_i(colnames(lipidmergeclean)[2:ncol(lipidmergeclean)],"-",3),
                    Cell = str_split_i(colnames(lipidmergeclean)[2:ncol(lipidmergeclean)],"-",1)) %>% 
  mutate(Cell = ifelse(Cell == "RCF", "CF", Cell))

lipidpca = ggplot(data_PCA, aes(PC1, PC2)) +
  geom_point(aes(colour = Cell, shape = Oxygen), size = 3) +
  ggtitle("Lipidomics") +
  xlab(paste0("PC1, VarExp: ", percent_var[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percent_var[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16)) + 
  scale_color_manual(values = c(CPC = "#F3C18A", MSC = "#A8A8CD", CEC = "#95C1AC", CF = "#E389B6"))+
  geom_hline(yintercept=0, color = "grey", size=0.5) + 
  geom_vline(xintercept=0, color = "grey", size=0.5) +
  scale_x_continuous(breaks = seq(-90,90,3)) +
  scale_y_continuous(breaks = seq(-90,90,3)) 
lipidpca
ggsave(lipidpca, file = paste0("./plots/lipid/lipidPCA_",Sys.Date(),".png"), width = 5, height = 4)

