library(ggplot2)
library(gggenes)
library(ORFik)

plot_genomes <- function(fasta, minimumLength=100){
  orfs <- findORFsFasta(filePath = fasta, minimumLength = minimumLength)
  
  df <- as.data.frame(orfs)
  
  df %<>% 
    group_by(seqnames) %>% 
    slice(which.max(width)) %>% 
    mutate(length=as.integer(str_extract(seqnames, "(?<=_length_)[0-9]+")))
  
  p <- ggplot(df, aes(y = seqnames)) +
    geom_segment(aes(x=0, xend=length, yend=seqnames))+
    geom_gene_arrow(aes(xmin=start, xmax=end)) +
    facet_wrap(~ seqnames, scales="free_y", ncol = 1)+
    labs(x="", y="")+
    theme_genes()+
    theme(axis.text.y = element_blank(), 
          panel.grid.major.y = element_blank())
  return(p)
}

#WMV6
plot_genomes("data/WMV6_MEMO067.fasta")

#WMV4
#MEMO061
plot_genomes("data/WMV4_MEMO061.fasta")

#MEMO126
plot_genomes("data/WMV4_MEMO126.fasta")

#AJOV
plot_genomes("data/AAOLV_MEMO117.fasta")

AJOV_orfs <- findORFsFasta(filePath = "data/AAOLV.fasta", stopCodon = c("TAG|TAA|TGA"))

AJOV_df <- as.data.frame(AJOV_orfs)

AJOV_df %<>% 
  group_by(seqnames) %>% 
  slice(which.max(width)) %>% 
  mutate(length=as.integer(str_extract(seqnames, "(?<=_length_)[0-9]+")), 
         start=case_when(length==2300 ~ 126, 
                         T ~ start),
         end=case_when(length==2300 ~ 2300, 
                       T ~ end))

ggplot(AJOV_df, aes(y = seqnames)) +
  geom_segment(aes(x=0, xend=length, yend=seqnames))+
  geom_gene_arrow(aes(xmin=start, xmax=end)) +
  facet_wrap(~ seqnames, scales="free", ncol = 1)+
  theme_genes()
