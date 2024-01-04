library(ggplot2)
library(gggenes)
library(tidyverse)
library(patchwork)

orthomyxo <- read_tsv("data/orthomyxo_genome_maps.tsv")

split_df_list <- split(orthomyxo, orthomyxo$Name)

WMV4 <- split_df_list$`Wuhan Mosquito virus 4`
WMV6 <- split_df_list$`Wuhan Mosquito virus 6`
AJOV <- split_df_list$`Aedes japonicus orthomyxovirus`

get_annotations <- function(df){
  df <- df %>% 
    mutate(Protein=gsub(" \\[.*\\]$", "", Protein),
           Protein=gsub("^.*: ", "", Protein),
           color=case_when(Protein %in% c("PA", "PB1", "PB2") ~ "RdRP",
                         str_detect(Protein, "hypothetical") ~ "hypothetical proteins",
                         T ~ Protein),
           Contig = reorder(Contig, -as.numeric(length))
    )
  return(df)
}

plot_genomes <- function(df){
  p <- ggplot(df, aes(y = Contig, fill=color, label=Protein, xmin=start, xmax=end)) +
    geom_segment(aes(x=0, xend=length, yend=Contig))+
    geom_gene_arrow() +
    geom_gene_label()+
    facet_wrap(~ Contig, scales="free_y", ncol = 1)+
    labs(x="", y="")+
    theme_genes()+
    theme(axis.text.y = element_blank(), 
          panel.grid.major.y = element_blank())+
    scale_fill_brewer(palette = "Accent", name="")
  return(p)
}

WMV4 <- get_annotations(WMV4)
a <- plot_genomes(WMV4)+
  ggtitle("Wuhan Mosquito Virus 4")

WMV6 <- get_annotations(WMV6)
b <- plot_genomes(WMV6)+
  ggtitle("Wuhan Mosquito Virus 6")

AJOV <- get_annotations(AJOV)
c <- plot_genomes(AJOV)+
  ggtitle("Aedes japonicus orthomyxovirus")

combined <- a | b | c
combined + 
  plot_layout(guides = "collect")
ggsave("figures/orthomyxo_genome_layout.pdf", dpi=300, width=11, height = 4)

#library(ORFik)
#
#get_annotations <- function(file){
#  annotations <- read_tsv(file)
#  annotations <- annotations %>% 
#    select(Query, Description) %>% 
#    mutate(Description=gsub(" \\[.*\\]$", "", Description),
#           Description=gsub("^.*: ", "", Description),
#           color=case_when(Description %in% c("PA", "PB1", "PB2") ~ "RdRP",
#                         str_detect(Description, "hypothetical") ~ "hypothetical proteins",
#                         T ~ Description))
#  return(annotations)
#}
#
#get_orfs <- function(fasta, minimumLength=100){
#  orfs <- findORFsFasta(filePath = fasta, minimumLength = minimumLength)
#  
#  df <- as.data.frame(orfs)
#  return(df)
#}
#
#merge_orf_anno <- function(df, anno){
#  df %<>% 
#    group_by(seqnames) %>% 
#    dplyr::slice(which.max(width)) %>% 
#    mutate(length=as.integer(str_extract(seqnames, "(?<=_length_)[0-9]+")))
#  
#  df %<>%
#    left_join(anno, by=join_by("seqnames"=="Query"))
#  
#  df$seqnames <- reorder(df$seqnames, rev(df$length))
#  
#  return(df)
#}
#
#plot_genomes <- function(df){
#  p <- ggplot(df, aes(y = seqnames, fill=color, label=Description, xmin=start, xmax=end)) +
#    geom_segment(aes(x=0, xend=length, yend=seqnames))+
#    geom_gene_arrow() +
#    geom_gene_label()+
#    facet_wrap(~ seqnames, scales="free_y", ncol = 1)+
#    labs(x="", y="")+
#    theme_genes()+
#    theme(axis.text.y = element_blank(), 
#          panel.grid.major.y = element_blank())
#  return(p)
#}
#
##WMV6
#WMV6_orfs <- get_orfs("data/WMV6_MEMO067.fasta")
#WMV6_anno <- get_annotations("data/")
#WMV6 <- merge_orf_anno()
#plot_genomes(WMV6)
#
##WMV4
##MEMO061
#WMV4_orfs <- get_orfs("data/WMV4_MEMO061.fasta")
#WMV4_anno <- get_annotations("data/WMV4_MEMO061.tsv")
#WMV4 <- merge_orf_anno(WMV4_orfs, WMV4_anno)
#plot_genomes(WMV4)
#
##MEMO126
#plot_genomes("data/WMV4_MEMO126.fasta")
#
##AJOV
#AJOV_anno <- get_annotations("data/AAOLV.tsv")
#plot_genomes("data/AAOLV_MEMO117.fasta", anno = AJOV_anno)
#
#AJOV_orfs <- findORFsFasta(filePath = "data/AAOLV.fasta", stopCodon = c("TAG|TAA|TGA"))
#
#AJOV_df <- as.data.frame(AJOV_orfs)
#
#AJOV_df %<>% 
#  group_by(seqnames) %>% 
#  slice(which.max(width)) %>% 
#  mutate(length=as.integer(str_extract(seqnames, "(?<=_length_)[0-9]+")), 
#         start=case_when(length==2300 ~ 126, 
#                         T ~ start),
#         end=case_when(length==2300 ~ 2300, 
#                       T ~ end))
#AJOV_df <- left_join(AJOV_df, AJOV_anno, by=join_by("seqnames"=="Query"))
#
#AJOV_df$seqnames <- reorder(AJOV_df$seqnames, rev(AJOV_df$length))
#
#plot_genomes(AJOV_df)+
#  scale_fill_brewer(palette = "Set3")

