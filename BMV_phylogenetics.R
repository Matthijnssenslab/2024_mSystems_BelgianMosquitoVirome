phylo_path <- here("phylogenetics")
phylo_packages <- c('phytools', 'ggtree', 'treeio', 'rphylopic',
                    'ggimage', 'aplot', 'rcartocolor')
lapply(phylo_packages, library, character.only = TRUE)

#' Determine shape based on which classification was used
shape <- c(ICTV=16, NCBI=15, NODE=17)
shape
#+ include=FALSE
# full_genomes <- c("Culex Iflavi-like virus 4", #Iflaviridae
#                   "Xanthi chryso-like virus", #Chrysoviridae
#                   "Daeseongdong virus 2", "Culex mosquito virus 1", #Nodaviridae
#                   "Wuhan Mosquito Virus 6", "Wuhan Mosquito Virus 4", "Guadeloupe mosquito quaranja-like virus 1", #Orthomyxo
#                   "Culex orthophasmavirus sp.", "Aedes orthophasmavirus sp.", #Phasmaviridae
#                   "Orbivirus sp.", "Valmbacken virus", #Reoviridae
#                   "Hubei mosquito virus 4", #Tymoviridae
#                   "Yongsan picorna-like virus 2", #Picorna
#                   "Culex-associated Tombus-like virus",
#                   "Culex negev-like virus 1", "Negev-like virus 174", "Rinkaby virus") #Negeviridae
# fg <- list()
# for (i in full_genomes) {
#   print(i)
#   print(taxa_names(phyloseq::subset_taxa(BMV.V2, Species == i)))
#   fg[[i]] <- taxa_names(phyloseq::subset_taxa(BMV.V2, Species == i))
# }
# 
# taxa_names(phyloseq::subset_taxa(BMV.V2, Order=="Picornavirales"))
# taxa_names(phyloseq::subset_taxa(BMV.V2, Species=="Hubei odonate virus 15"))
# 
# BMV_taxtable <- psmelt(BMV_final) %>% 
#   select(OTU, Phylum, Class, Order, Family, Genus, Species) %>% 
#   distinct()
# virfam <- BMV_taxtable %>% 
#   select(Family, Species)
# for (i in BMV_taxtable$Species) {
#   print(i)
#   print(taxa_names(phyloseq::subset_taxa(BMV.V2, Species == i)))
#   fg[[i]] <- taxa_names(phyloseq::subset_taxa(BMV.V2, Species == i))
# }
# fg$`Culex Iflavi-like virus 4`
# segments <- stack(fg) %>% 
#   mutate(Length=gsub(x=gsub(x=values, pattern = ".*_length_", replacement = "", perl = T), 
#                      pattern = "_cov.*", replacement = "", perl = T)) %>% 
#   rename(Species=ind, OTU=values) 
# 
# seg_count <- segments %>% group_by(Species) %>% count()
# 
# seg_max <- segments%>% 
#   select(OTU, Species, Length) %>% 
#   group_by(Species) %>% 
#   slice(which.max(Length)) %>% 
#   arrange(Species)
#BMV_taxtable <- as.data.frame(tax_table(BMV_final))

#setwd("~/OneDrive - KU Leuven/Documents/Manuscripts/Belgian mosquitoes (2021)/analysis/phylogenetics/")


# pal.safe(stepped(n=24))
# pal.safe(viridis::plasma(n = 5, begin=0.2, end=0.9))
# pal.safe(viridis::viridis(n = 5, begin=0.2, end=0.9))
# pal.safe(viridis::mako(n = 5, begin=0.2, end=0.9))
# pal.safe(viridis::cividis(n = 5, begin=0, end=1))
# pal.safe(viridis::rocket(n = 5, begin=0.3, end=0.9))
# pal.safe(viridis::turbo(n = 5, begin=0.6, end=0.9))

#' ## Bunyavirales

meta.bunya <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Bunyavirales/"), pattern = "*.csv", full.names = T), read_csv))
bunya.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Bunyavirales/"), pattern = "*_ictv.csv", full.names = T), read_csv))
bunya <- read.iqtree(here(phylo_path, "Bunyavirales/bunyavirales.treefile"))
bunya@phylo <- phytools::midpoint.root(bunya@phylo)

meta.bunya <- meta.bunya %>% 
  mutate(Classification = case_when(Accession %in% bunya.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

bunya@phylo$tip.label <- gsub(bunya@phylo$tip.label, pattern = "_", replacement=" ")
bunya@phylo$tip.label <- gsub(bunya@phylo$tip.label, pattern = "P ", replacement="P_")
bunya@phylo$tip.label <- gsub(bunya@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

bunya_clean <- clean_phylometa(bunya, metadata = meta.bunya)
bunya_ictv_clean <- clean_phylometa(bunya, metadata = bunya.ictv)

bunya_tree <- group_phylotaxa(bunya, bunya_clean, group = c("Family","Genus", "Host", "Geo_Location", "Classification"))
bunya_tree

#' Give clean metadata ictv!
bunyacladedf <- mrca_ictv(tree = bunya_tree, meta = bunya_ictv_clean, 
                          group="Genus")

famvec <- sort(unique(bunya_clean$Family))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.4, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

bunya_tree@phylo$tip.label <- gsub(bunya_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(bunya_tree, col=col, shape=shape, cladedf = bunyacladedf, 
                    labels=c(famvec, "NODE"),
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F, )+
  #geom_nodepoint(aes(subset= node %in% bunyacladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,6)

bunyaplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4))
bunyaplot

ggsave(filename = here(phylo_path, "Bunyavirales/bunyaviridae.pdf"), bunyaplot, dpi=1200, width=6.875)

bunyaplot_minimal <- plot_minimal_phylotree(bunya_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,6.5)
bunyaplot_minimal

ggsave(filename = here(phylo_path, "Bunyavirales/bunyaviridae_minimal.pdf"), bunyaplot_minimal, dpi=1200, width=6.875)

#' ## Endornaviridae

meta.endorna <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Endornaviridae/"), pattern = "*.csv", full.names = T), read_csv))
endorna.ictv <- read_csv(here(phylo_path, "Endornaviridae/endornaviridae_ictv.csv"))
endorna <- read.iqtree(here(phylo_path, "Endornaviridae/endornaviridae.treefile"))
endorna@phylo <- phytools::midpoint.root(endorna@phylo)

meta.endorna <- meta.endorna %>% 
  mutate(Classification = case_when(Accession %in% endorna.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

endorna@phylo$tip.label <- gsub(endorna@phylo$tip.label, pattern = "_", replacement=" ")
endorna@phylo$tip.label <- gsub(endorna@phylo$tip.label, pattern = "P ", replacement="P_")
endorna@phylo$tip.label <- gsub(endorna@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

endorna_clean <- clean_phylometa(endorna, metadata = meta.endorna)
endorna_ictv_clean <- clean_phylometa(endorna, metadata = endorna.ictv)

endorna_tree <- group_phylotaxa(endorna, endorna_clean, group = c("Family","Genus", "Host", "Geo_Location", "Classification"))
endorna_tree

#' Give clean metadata ictv!
endornacladedf <- mrca_ictv(tree= endorna_tree, meta = endorna_ictv_clean,
                            group = "Genus")

famvec <- sort(unique(endorna_clean$Genus))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.7, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

endorna_tree@phylo$tip.label <- gsub(endorna_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(endorna_tree, col=col, shape=shape, cladedf = endornacladedf, labels=c(famvec, "NODE"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F, group="Genus")+
  #geom_nodepoint(aes(subset= node %in% endornacladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,4)

endornaplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(0.1,0.15))
endornaplot

ggsave(filename = here(phylo_path, "Endornaviridae/endornaviridae.pdf"), endornaplot, dpi=1200, width=6.875)

endornaplot_minimal <- plot_minimal_phylotree(endorna_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,4.5)
endornaplot_minimal

ggsave(filename = here(phylo_path, "Endornaviridae/endornaviridae_minimal.pdf"), endornaplot_minimal, dpi=1200, width=6.875)

#' ## Ghabrivirales
meta.ghabri <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Toti_Chryso/"), pattern = "*.csv", full.names = T), read_csv))
ghabri.ictv <- read_csv(here(phylo_path, "Toti_Chryso/ghabrivirales_ictv.csv"))
ghabri <- read.iqtree(here(phylo_path, "Toti_Chryso/ghabrivirales.treefile"))
ghabri@phylo <- phytools::midpoint.root(ghabri@phylo)

meta.ghabri <- meta.ghabri %>% 
  mutate(Classification = case_when(Accession %in% ghabri.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

ghabri@phylo$tip.label <- gsub(ghabri@phylo$tip.label, pattern = "_", replacement=" ")
ghabri@phylo$tip.label <- gsub(ghabri@phylo$tip.label, pattern = "P ", replacement="P_")
ghabri@phylo$tip.label <- gsub(ghabri@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

ghabri_clean <- clean_phylometa(ghabri, metadata = meta.ghabri)
ghabri_ictv_clean <- clean_phylometa(ghabri, metadata = ghabri.ictv)

ghabri_tree <- group_phylotaxa(ghabri, ghabri_clean, group = c("Family","Genus", "Host", "Geo_Location", "Classification"))
ghabri_tree

#' Give clean metadata ictv!
ghabricladedf <- mrca_ictv(tree=ghabri_tree, meta =ghabri_ictv_clean, 
                           group="Genus", subset = c("unclassified", "Chrysovirus"))

famvec <- sort(unique(ghabri_clean$Family))
famvec
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.2, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

ghabri_tree@phylo$tip.label <- gsub(ghabri_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(ghabri_tree, col=col, shape=shape, cladedf = ghabricladedf, labels=c(famvec, "NODE"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F)+
  #geom_nodepoint(aes(subset= node %in% ghabricladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,4)
ghabriplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(0.15,0.15))
ghabriplot
ggsave(filename = here(phylo_path, "Toti_Chryso/ghabrivirales.pdf"), ghabriplot, dpi=1200, width=6.875)

ghabriplot_minimal <- plot_minimal_phylotree(ghabri_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,4.5)
ghabriplot_minimal

ggsave(filename = here(phylo_path, "Toti_Chryso/ghabrivirales_minimal.pdf"), ghabriplot_minimal, dpi=1200, width=6.875)

#' ## Nodamuvirales
meta.noda <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Nodaviridae/"), pattern = "*.csv", full.names = T), read_csv))
noda.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Nodaviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
noda <- read.iqtree(here(phylo_path, "Nodaviridae/nodamuvirales.treefile"))
noda@phylo <- phytools::midpoint.root(noda@phylo)

meta.noda <- meta.noda %>% 
  mutate(Classification = case_when(Accession %in% noda.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

noda@phylo$tip.label <- gsub(noda@phylo$tip.label, pattern = "_", replacement=" ")
noda@phylo$tip.label <- gsub(noda@phylo$tip.label, pattern = "P ", replacement="P_")
noda@phylo$tip.label <- gsub(noda@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

noda_clean <- clean_phylometa(noda, metadata = meta.noda)
noda_ictv_clean <- clean_phylometa(noda, metadata = noda.ictv)

noda_tree <- group_phylotaxa(noda, noda_clean, group = c("Family","Genus", "Host", "Geo_Location", "Classification"))
noda_tree

#' Give clean metadata ictv!
nodacladedf <- mrca_ictv(tree = noda_tree, meta = noda_ictv_clean, group = "Genus")

famvec <- sort(unique(noda_clean$Family))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.4, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

noda_tree@phylo$tip.label <- gsub(noda_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(noda_tree, col=col, shape=shape, cladedf = nodacladedf, labels=c(famvec, "NODE"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F)+
  #geom_nodepoint(aes(subset= node %in% nodacladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,4.5)

nodaplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(0.1,0.15))
nodaplot
ggsave(filename = here(phylo_path, "Nodaviridae/nodamuvirales.pdf"), nodaplot, dpi=1200, width=6.875)

nodaplot_minimal <- plot_minimal_phylotree(noda_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,5)
nodaplot_minimal

ggsave(filename = here(phylo_path, "Nodaviridae/nodamuvirales_minimal.pdf"), nodaplot_minimal, dpi=1200, width=6.875)

#' ## Orthomyxoviridae

meta.orthomyxo <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Orthomyxoviridae/"), pattern = "*.csv", full.names = T), read_csv))
orthomyxo.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Orthomyxoviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
orthomyxo <- read.iqtree(here(phylo_path, "Orthomyxoviridae/orthomyxoviridae.treefile"))
orthomyxo@phylo <- phytools::midpoint.root(orthomyxo@phylo)

meta.orthomyxo <- meta.orthomyxo %>% 
  mutate(Classification = case_when(Accession %in% orthomyxo.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

orthomyxo@phylo$tip.label <- gsub(orthomyxo@phylo$tip.label, pattern = "_", replacement=" ")
orthomyxo@phylo$tip.label <- gsub(orthomyxo@phylo$tip.label, pattern = "P ", replacement="P_")
orthomyxo@phylo$tip.label <- gsub(orthomyxo@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

orthomyxo_clean <- clean_phylometa(orthomyxo, metadata = meta.orthomyxo)
orthomyxo_ictv_clean <- clean_phylometa(orthomyxo, metadata = orthomyxo.ictv)

orthomyxo_tree <- group_phylotaxa(orthomyxo, orthomyxo_clean, group = c("Family","Genus", "Host", "Geo_Location", "Classification"))
orthomyxo_tree

#' Give clean metadata ictv!
orthocladedf <- mrca_ictv(tree=orthomyxo_tree, meta = orthomyxo_ictv_clean, 
                          group="Genus")

famvec <- sort(unique(orthomyxo_clean$Genus))
#col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
col <- c(rcartocolor::carto_pal(n = length(famvec), "Vivid"))
col <- rev(col)[-1]
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

orthomyxo_tree@phylo$tip.label <- gsub(orthomyxo_tree@phylo$tip.label, pattern = "lcl [ORF]*[0-9]* *", replacement="", perl = T)

p <- plot_phylotree(orthomyxo_tree, col=col, shape=shape, cladedf = orthocladedf, labels=c(famvec, "NODE"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F, group="Genus")+
  #geom_nodepoint(aes(subset= node %in% orthocladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,4.5)

orthoplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4))
orthoplot
ggsave(filename = here(phylo_path, "Orthomyxoviridae/orthomyxoviridae.pdf"), orthoplot, dpi=1200, width=6.875)

orthoplot_minimal <- plot_minimal_phylotree(orthomyxo_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,5)
orthoplot_minimal

ggsave(filename = here(phylo_path, "Orthomyxoviridae/orthomyxoviridae_minimal.pdf"), orthoplot_minimal, dpi=1200, width=6.875)

#' ## Picornavirales

meta.picorna <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Picorna/"), pattern = "*.csv", full.names = T), read_csv))
picorna.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Picorna/"), pattern = "*_ictv.csv", full.names = T), read_csv))
picorna <- read.iqtree(here(phylo_path, "Picorna/picornavirales.treefile"))
picorna@phylo <- phytools::midpoint.root(picorna@phylo)

meta.picorna <- meta.picorna %>% 
  mutate(Classification = case_when(Accession %in% picorna.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

picorna@phylo$tip.label <- gsub(picorna@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

picorna_clean <- clean_phylometa(picorna, metadata = meta.picorna)
picorna_ictv_clean <- clean_phylometa(picorna, metadata = picorna.ictv)

picorna_tree <- group_phylotaxa(picorna, picorna_clean, 
                                group = c("Family","Genus", "Host", "Geo_Location", "Classification"))

#' Give clean metadata ictv!
picornacladedf <- mrca_ictv(tree=picorna_tree, meta = picorna_ictv_clean, 
                            group="Genus", subset="unclassified")#, subset = "Cripavirus")

famvec <- sort(unique(picorna_clean$Family))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.2, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

picorna_tree@phylo$tip.label <- gsub(picorna_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)


p <- plot_phylotree(picorna_tree, col=col, shape=shape, cladedf = picornacladedf, labels=c(famvec, "NODE"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F)+
  #geom_nodepoint(aes(subset= node %in% picornacladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,4)

picornaplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(.8,.15))
picornaplot
#+ echo=FALSE
# picplot <- picornaplot+
#   geom_tiplab(data=picorna_tree, aes(subset=Genus == "Cripavirus", label="Cripavirus"), 
#                size=2.5, color="red", offset = .05)+
#   geom_nodelab(aes(label=as.numeric(label), subset = !is.na(as.numeric(label)) & as.numeric(label) > 0), 
#                size=2, hjust=-.01)+
#   theme(legend.position = c(.8,.15))
# picplot

ggsave(filename = here(phylo_path, "Picorna/picornavirales.pdf"), picornaplot, dpi=1200, width=6.875)

#+ echo=TRUE
picornaplot_minimal <- plot_minimal_phylotree(picorna_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,4.5)
picornaplot_minimal

ggsave(filename = here(phylo_path, "Picorna/picornavirales_minimal.pdf"), picornaplot_minimal, dpi=1200, width=6.875)

#' ## Reovirales

meta.reo <- do.call(rbind, lapply(list.files(path=here::here(phylo_path, "Reoviridae/"), pattern = "*.csv", full.names = T), read_csv))
reo.ictv <- do.call(rbind, lapply(list.files(path=here::here(phylo_path, "Reoviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
reo <- read.iqtree(here(phylo_path, "Reoviridae/reoviridae.treefile"))
reo@phylo <- phytools::midpoint.root(reo@phylo)

meta.reo <- meta.reo %>% 
  mutate(Classification = case_when(Accession %in% reo.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

reo@phylo$tip.label <- gsub(reo@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

reo_clean <- clean_phylometa(reo, metadata = meta.reo)
reo_ictv_clean <- clean_phylometa(reo, metadata = reo.ictv)

reo_tree <- group_phylotaxa(reo, reo_clean, group = c("Family","Genus", "Host", "Geo_Location", "Classification"))
reo_tree@phylo$tip.label

#' Give clean metadata ictv!
reocladedf <- mrca_ictv(tree=reo_tree, meta =reo_ictv_clean, 
                        group="Genus")

famvec <- sort(unique(reo_clean$Family))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
#col <- c(rcartocolor::carto_pal(n = length(famvec), "Vivid"))
#col <- rev(col)[-1]
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

reo_tree@phylo$tip.label <- gsub(reo_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)
#+echo=FALSE
reo_tree@phylo$tip.label <- gsub(reo_tree@phylo$tip.label, pattern = " pasted", replacement="", perl = T)
#+echo=TRUE
p <- plot_phylotree(reo_tree, col=col, shape=shape, cladedf = reocladedf, labels=c(famvec, "NODE"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F, group="Family")+
  #geom_nodepoint(aes(subset = node %in% reocladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,7)
reoplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(0.7,.3))
reoplot

ggsave(filename = here(phylo_path, "Reoviridae/reoviridae.pdf"), reoplot, dpi=1200, width=6.875)

reoplot_minimal <- plot_minimal_phylotree(reo_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,9)
reoplot_minimal

ggsave(filename = here(phylo_path, "Reoviridae/reoviridae_minimal.pdf"), reoplot_minimal, dpi=1200, width=6.875)

#' ## Tombusviridae
#' Rerun tree with hmv4!
meta.tombus <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Tombusviridae/"), pattern = "*.csv", full.names = T), read_csv))
tombus.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Tombusviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
tombus <- read.iqtree(here(phylo_path, "Tombusviridae/tombusviridae2.treefile"))
tombus@phylo <- phytools::midpoint.root(tombus@phylo)

tombus@phylo$tip.label <- gsub(tombus@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")


meta.tombus <- meta.tombus %>% 
  mutate(Classification = case_when(Accession %in% tombus.ictv$Accession ~ "ICTV",
                                    TRUE ~ "NCBI"))

tombus_clean <- clean_phylometa(tombus, metadata = meta.tombus)
tombus_ictv_clean <- clean_phylometa(tombus, metadata = tombus.ictv)

tombus_tree <- group_phylotaxa(tombus, tombus_clean, 
                               group = c("Family","Genus", "Host", "Geo_Location", "Classification"))

#' Give clean metadata ictv!
tombuscladedf <- mrca_ictv(tree = tombus_tree, meta = tombus_ictv_clean, 
                           group="Genus")

famvec <- sort(unique(tombus_clean$Genus))

#col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
col <- c(rcartocolor::carto_pal(n = length(famvec), "Prism"))
col <- rev(col)[-1]
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

tombus_tree@phylo$tip.label <- gsub(tombus_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)
tombus_tree

p <- plot_phylotree(tombus_tree, col=col, shape=shape, cladedf = tombuscladedf,
                    labels=c(famvec, "NODE"),
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F,
                    group="Genus")+
  #geom_nodepoint(aes(subset= node %in% tombuscladedf$node), size=1, shape=18)+
  #geom_nodepoint(aes(fill=as.numeric(label), subset = !is.na(as.numeric(label))), shape=23, color="transparent", size=1)+
  #scale_fill_gradientn(colours = c("red2","orange","gold1","forestgreen")) +
  scale_y_reverse()+
  xlim(-.05,4)
tombusplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(0.1,.2))+
  guides(shape="none")
tombusplot

ggsave(filename = here(phylo_path, "Tombusviridae/tombusviridae.pdf"), tombusplot, dpi=1200, width=6.875)
#ggsave(filename = here(phylo_path, "Tombusviridae/tombusviridae_meeting.pdf"), tombusplot+theme(legend.position="left"), dpi=300, width=5.6, height=3.5)

tombusplot_minimal <- plot_minimal_phylotree(tombus_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,5)
tombusplot_minimal

ggsave(filename = here(phylo_path, "Tombusviridae/tombusviridae_minimal.pdf"), tombusplot_minimal, dpi=1200, width=6.875)

#' ## Negevirus

meta.negev <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Negevirus/"), pattern = "*.csv", full.names = T), read_csv))
#tombus.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Tombusviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
negev <- read.iqtree(here(phylo_path, "Negevirus/negeviruses.treefile"))
negev@phylo <- phytools::midpoint.root(negev@phylo)

meta.negev <- meta.negev %>% 
  mutate(Classification = "NCBI")

negev@phylo$tip.label <- gsub(negev@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

negev_clean <- clean_phylometa_negev(negev, metadata = meta.negev)

#negev_ictv_clean <- clean_phylometa(negev, metadata = negev.ictv)

negev_tree <- group_phylotaxa(negev, negev_clean, 
                              group = c("Family","Genus", "Host", "Geo_Location", "Classification"))

#' Give clean metadata ictv!
# negevcladedf <- mrca_ictv(tree = negev_tree, meta = negev_ictv_clean, 
#                            group="Genus")

famvec <- sort(unique(negev_clean$Genus))

col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")

negev_tree@phylo$tip.label <- gsub(negev_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)
negev_tree

p <- plot_phylotree(negev_tree, col=col, shape=shape,
                    labels=c(famvec, "NODE"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F,
                    group="Genus")+
  #geom_nodepoint(aes(subset= node %in% negevcladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,3.5)

negevplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(0.2,.2))
negevplot

ggsave(filename = here(phylo_path, "Negevirus/negeviruses.pdf"), negevplot, dpi=1200, width=6.875)

negevplot_minimal <- plot_minimal_phylotree(negev_tree, align=T)+
  scale_y_reverse()+
  xlim(-.05,5)
negevplot_minimal

ggsave(filename = here(phylo_path, "Negevirus/negeviruses_minimal.pdf"), negevplot_minimal, dpi=1200, width=6.875)


#' ## Combine phylogenetic trees
# p1 <- ggarrange(endornaplot, orthoplot,
#                 bunyaplot, picornaplot,
#                 ghabriplot, reoplot,
#                 ncol=2, nrow=3, labels="AUTO",
#                 heights = c(.2,.4,.4))

# p1 <- ggarrange(endornaplot, orthoplot,
#                 bunyaplot, picornaplot,
#                 ghabriplot, reoplot,
#                 tombusplot, negevplot,
#                 ncol=2, nrow=4, labels="AUTO",
#                 heights = c(.15,.25,.25,.25))
# p1

p1 <- ggarrange(negevplot, endornaplot, orthoplot,
                picornaplot, bunyaplot, ghabriplot, 
                reoplot, tombusplot, nodaplot,
                ncol=3, nrow=3, labels="AUTO",
                heights = c(.25, .37, .37),
                font.label = list(size=10))
#p1

#+ echo=TRUE, fig.width=20, fig.height=20
# ggarrange(p1, 
#           ggarrange(NULL, tombusplot+theme(legend.position = "right"), NULL, 
#                     widths=c(0.1,0.8,0.1), labels=c("","G",""), ncol=3), 
#           ncol=1, heights = c(.75, .25))
# p1 <- ggarrange(p1, 
#           ggarrange(reoplot, tombusplot, nodaplot,#+theme(legend.position = "right"), 
#                     labels=c("G","H", "I"), ncol=3), 
#           ncol=1, heights = c(.6, .4)) #0.6, 0.4, widths=c(0.4,0.4, 0.2)

#ggsave("figures/treeplots.pdf", width=17, height = 28, units = "in", dpi=300)
ggsave("figures/treeplots.pdf", width=20, height = 20, units = "in", dpi=600)
ggsave("figures/treeplots_small.pdf", width=6.875, height = 9.0625, units = "in", dpi=1200)

#' Minimal phylogenetic trees:
# p2 <- ggarrange(endornaplot_minimal, orthoplot_minimal,
#                 bunyaplot_minimal, picornaplot_minimal,
#                 ghabriplot_minimal, reoplot_minimal,
#                 ncol=2, nrow=3, labels="AUTO",
#                 heights = c(.2,.4,.4))

p2 <- ggarrange(endornaplot_minimal, orthoplot_minimal, negevplot_minimal,
                bunyaplot_minimal, picornaplot_minimal, ghabriplot_minimal,
                ncol=3, nrow=2, labels="AUTO",
                heights = c(.4,.6))

#+ echo=TRUE, fig.width=17, fig.height=28
# ggarrange(p2, 
#           ggarrange(NULL, tombusplot_minimal, NULL, 
#                     widths=c(0.1,0.8,0.1), labels=c("","G",""), ncol=3), 
#           ncol=1, heights = c(.75, .25))
#ggarrange(p2, 
#          ggarrange(reoplot_minimal, tombusplot_minimal, 
#                    widths=c(0.5,0.5), labels=c("G","H"), ncol=2), 
#          ncol=1, heights = c(.7, .3))

p2 <- ggarrange(negevplot_minimal, endornaplot_minimal, orthoplot_minimal,
                picornaplot_minimal, bunyaplot_minimal, ghabriplot_minimal, 
                reoplot_minimal, tombusplot_minimal, nodaplot_minimal,
                ncol=3, nrow=3, labels="AUTO",
                heights = c(.25, .37, .37),
                font.label = list(size=10))

#ggsave("figures/treeplots_minimal.pdf", width=17, height = 28, units = "in", dpi=300)
ggsave("figures/treeplots_minimal.pdf", width=20, height = 20, units = "in", dpi=600)
ggsave(plot=p2,filename="figures/treeplots_minimal_small.pdf", width=6.875, height = 9.0625, units = "in", dpi=600)