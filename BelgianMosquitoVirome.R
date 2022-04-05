#' ---
#' title: "Belgian Mosquito Virome"
#' author: "Lander De Coninck"
#' date: ""
#' output:
#'   html_document:
#'       df_print: paged
#'       number_sections: TRUE
#'       keep_md: no
#'       theme: default
#'       highlight: kate
#'       toc: true
#'       toc_depth: 3 
#'       toc_float:
#'          smooth_scroll: true
#' ---
#' R script to analyse the data on virome sequencing of single Belgian mosquitoes.
#+ include=FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, 
                      fig.width = 11.2, fig.height = 6.92, cache = F)#, fig.path='figures/')
#+ include=TRUE, results="hide"
necessary_packages <- c('tidyverse', 'data.table', 'scales', 'knitr', 'grid', 'magrittr', 'DT', 'here', 'gridtext',
                        'ggpubr', 'ggthemes', 'ggpmisc', 'reshape2', 'viridis', 'pals', 'circlize',
                        'vegan', 'phyloseq', 'metagenomeSeq', 'ComplexHeatmap', 'decontam')
lapply(necessary_packages, library, character.only = TRUE)
i_am("BelgianMosquitoVirome.R")
#setwd(here())
#setwd("/Users/lander/OneDrive - KU Leuven/Documents/Manuscripts/Belgian mosquitoes (2022)/BMV-analysis/")
source("BMV_functions.R")
#+ echo=FALSE
si <- sessioninfo::session_info()
pckgs <- map2(si$packages$package[si$packages$attached == T], 
              si$packages$loadedversion[si$packages$attached == T],
              ~ paste0(.x, " ", .y)) %>% 
  simplify()
print(si$platform$version)
print(paste0("Running under: ", si$platform$os))
print(pckgs)

#' ***
#' # Eukaryotic virome analysis
#' ## Prepare the data
#' ### Load the OTU table, taxonomy file and metadata into R
OTU <- read.table("data/abundance-table.txt", header=TRUE, row.names=1, sep="\t", dec=".")
names(OTU) <- gsub(x = names(OTU), pattern = "\\.", replacement = "-")
summary(rowSums(OTU))
summary(colSums(OTU))

tax <- read.table("data/BEmosq_classification-1000nt.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
meta <- read.table("data/BEmosq_metadata.csv", header=TRUE, row.names = 1, sep=";", dec=".")
meta <- cbind(Sample=rownames(meta), meta)
#+ echo=FALSE
datatable(meta)

#' ### Make a phyloseq object
OTU.UF <- otu_table(as.matrix(OTU), taxa_are_rows=T)
tax.UF <- tax_table(as.matrix(tax))
meta.UF <- sample_data(meta)

BMV <- phyloseq(OTU.UF, tax.UF, meta.UF)
BMV

#' ### Remove contamination
#' #### Visualize library sizes of samples and negative controls
decontam <- as.data.frame(sample_data(BMV))
decontam$LibrarySize <- sample_sums(BMV)
decontam <- decontam[order(decontam$LibrarySize),]
decontam$Index <- seq(nrow(decontam))
ggplot(data=decontam, aes(x=Index, y=LibrarySize, color=Control)) + 
  geom_point()+
  ggtitle("Library sizes")

#' #### Detect contaminants
sample_data(BMV)$is.neg <- sample_data(BMV)$Control == "Yes"
contamdf.prev <- isContaminant(BMV, method="prevalence", neg="is.neg")

#' **Number of contaminating contigs:**
table(contamdf.prev$contaminant)

#' #### Visualize prevalence of contaminants in samples and negative controls
ps.pa <- transform_sample_counts(BMV, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "Yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "No", ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  ggtitle("Prevalence of contaminants in samples vs NCs")

#' #### Remove negative controls and contaminants from phyloseq object
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, BMV)
ps.noncontam
ps.noncontam <- prune_samples(sample_data(BMV)$Control!='Yes', ps.noncontam)
BMV <- ps.noncontam
BMV

#' Remove negative controls from metadata table
meta <- meta[meta$SKA_Subspecies!="control",]

#' ### Subset virome data
#' #### Subset only viruses and remove prokaryotic viruses and EVEs
BMV.V <- subset_taxa(BMV, Kingdom=="Viruses")
BMV.V

EVE_phage <- c("Atrato Retro-like virus", "Gurupi chuvirus-like 1", "Aedes aegypti To virus 1",
  "Aedes aegypti To virus 2", "Guato virus", "Kaiowa virus", "Atrato Chu-like virus 1",
  "Chuvirus Mos8Chu0", "Chibugado virus", "Prokaryotic dsDNA virus sp.")

BMV.V2 <- subset_taxa(BMV.V, !is.element(Species, EVE_phage))
BMV.V2 <- subset_taxa(BMV.V2, Order!="Caudovirales")
BMV.V2 <- subset_taxa(BMV.V2, Phylum!="Phage")
BMV.V2

#' **Info on sample variables and level of taxonomic ranks:**
sample_variables(BMV.V2)
rank_names(BMV.V2)

#' #### Agglomerate taxa on viral species level
BMV_species <- BMV.V2 %>%
  tax_glom(taxrank = "Species")
BMV_species

vcount <- as.data.frame(colSums(otu_table(BMV_species) != 0))

names(vcount) <- "viral_count"

vcount %>% 
  group_by(viral_count) %>% 
  count() %>% 
  ggplot(aes(x=viral_count, y=n))+
  geom_col()+
  theme_bw()

df <- merge(meta, vcount, by=0) %>% 
  group_by(viral_count, SKA_Subspecies) %>% 
  count()

p <- ggplot(df, aes(x=viral_count, y=n, label=n, fill=SKA_Subspecies))+
  geom_col()+
  scale_x_continuous(n.breaks = max(df$viral_count), name = "# viral species")+
  scale_y_continuous(breaks=seq(0, 125, 25), name = "# mosquitoes")+
  scale_fill_viridis_d(begin=0.4, end=0.9, name="")+
  #geom_text(position =position_stack(vjust=.5), size=3)+
  theme_bw()+
  theme(legend.text = element_text(face="italic"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p

#+ eval=FALSE,echo=FALSE
ggarrange(NULL, ggarrange(NULL, p+theme(legend.position = "none"), ncol = 1, labels = c("B", "C")), ncol=2, 
          labels = c("A", NULL), widths = c(.6,.4))
ggsave("figures/figure1.pdf", dpi=300)

# library(grImport)
# PostScriptTrace("figures/species_assignment_clustermap.pdf", charpath = F)
# my_shape <- readPicture("species_assignment_clustermap.pdf.xml")
# 
# my_shape_grob <- pictureGrob(my_shape)
# #my_shape_grob
# ggsave("test.pdf", cowplot::plot_grid(my_shape_grob))
# 
# PostScriptTrace("figures/Belgium_subspecies.pdf", charpath = T)
# my_shape <- readPicture("Belgium_subspecies.pdf.xml")
# #grid.picture(my_shape)
# 
# my_shape_grob <- pictureGrob(my_shape)
# cowplot::plot_grid(my_shape_grob)
# ggsave("test.pdf", cowplot::plot_grid(my_shape_grob))

#' **Only keep taxa with more than 0 reads:**
BMV_species <- prune_taxa(taxa_sums(BMV_species) > 0, BMV_species)
BMV_species

#' **Only keep samples with more than 0 viral reads:**
BMV_final <- prune_samples(sample_sums(BMV_species) > 0, BMV_species)
BMV_final

#' ### Prepare relative abundance data
#' **Agglomerate taxa on viral Family level:**
BMV_family <- BMV_final %>%
  tax_glom(taxrank = "Family")
BMV_family

#' **Calculate relative abundance of viral families:**
BMV_family_re.abund <- BMV_family %>%                   
  transform_sample_counts(function(x) {round(100*(x/sum(x)),5)} )
BMV_family_re.abund

#' **Count mosquito species per location:**
sp_count <- psmelt(BMV_family_re.abund) %>% 
  select(Sample,SKA_Subspecies,Municipality) %>% 
  distinct() %>% 
  group_by(Municipality) %>% 
  count(SKA_Subspecies) %>% 
  pivot_wider(names_from=Municipality, values_from = n) %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
  rename("Villers_Le_Bouillet" = `Villers-Le-Bouillet`)
#+ echo=FALSE
datatable(sp_count)

#' **Create list of locations:**
loclist<-colnames(sp_count)[2:ncol(sp_count)]
rev(loclist)

#' **Create list of barplots with number of mosquitoes per location:**
plist <- list(ggplot()+theme_void()+#ggtitle("Mosquito species")+
                theme(title = element_text(size = 6)))
for (i in rev(loclist)){
  plist[[i]] <- ggplot(sp_count, aes_string(x = "SKA_Subspecies", y=i, fill="SKA_Subspecies"))+
    geom_col()+
    #scale_color_viridis_d(begin=0.4, end=0.9)+
    scale_fill_viridis_d(begin=0.4, end=0.9, name="")+
    theme_classic()+
    theme(legend.position = "none",
          axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size=10),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(colour = 'black', size = 0.25),
          axis.ticks = element_line(colour = "black", size = 0.25))+
    #ggtitle(gsub("_","-",i))+
    scale_y_continuous(limits = c(0,20), expand = c(0,0), position = "right")
}
plist[10] <- list(ggplot()+theme_void())
ggarrange(plotlist = plist[2:9], labels = gsub("_","-",rev(loclist)), 
          font.label = list(size=10, face="plain"), hjust = 0, vjust = 1)

#' **Create colorpalette for viral families:**
BMV_smelt <- psmelt(BMV_final)
FamLevel <- levels(as.factor(BMV_smelt[(BMV_smelt$Family!="unclassified"),]$Family))
FamLevel <- c(FamLevel, "unclassified")
BMV_smelt$Family <- factor(BMV_smelt$Family, levels = FamLevel)

myFamCol <- c(stepped(n=20), "#E6550D", "#FD8D3C", "#FDAE6B","lightgrey")
names(myFamCol) <- levels(as.factor(BMV_smelt$Family))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(myFamCol, main = "Viral families")

#+ echo=TRUE
#' **Create colors for each mosquito species:**
myColors <- viridis::viridis(4,1, begin=0.4, end = 0.9, direction = 1)
names(myColors) <- levels(as.factor(BMV_smelt$SKA_Subspecies))

#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(myColors, main="Mosquito species")

#' ### Prepare heatmap data
#' **Make table with fungi reads**:
BMV_fungi <- subset_taxa(BMV, Phylum %in% c("Ascomycota", "Basidiomycota", 
                                            "Chytridiomycota", "Glomeromycota",
                                            "Zygomycota","Neocallimastigomycota",
                                            "Microsporidia"))
BMV_fungi

BMV_fungi <- tax_glom(BMV_fungi, taxrank = 'Phylum')

sample_variables(BMV_fungi)
fungi <- as.data.frame(otu_table(BMV_fungi))
fungi2 <- fungi %>% 
  mutate(Phylum=rownames(.), .before=1) %>% 
  pivot_longer(cols = 2:198, names_to="Sample", values_to = "fungi_reads") %>% 
  pivot_wider(names_from = Phylum, values_from = fungi_reads)
fungi2

tax_table(BMV_fungi)

# fungi.df <- pivot_longer(fungi, cols = 1:197, names_to="Sample", values_to = "fungi_reads")
# fungi.df
# meta <- left_join(meta, fungi.df, by="Sample")
# rownames(meta) <- meta$Sample
# sample_data(BMV_final) <- data.frame(meta)

#' **Convert phyloseq object to metagenomeseq object to make heatmap:**
BMV_metaseq <- phyloseq_to_metagenomeSeq(BMV_final)

#' **Aggregate by species:**
BMV_metaseq_species <- aggregateByTaxonomy(BMV_metaseq, lvl = "Species", norm = F, 
                                           aggfun = colSums, out = "MRexperiment", alternate = T)

#' **Count number of unique viral species:**
n_species <- length(unique(featureData(BMV_metaseq_species)$Species))

#' **Assign mosquito species for each sample:**
mosquito_species <- pData(BMV_metaseq_species)$SKA_Subspecies

#' **Assign colors to location:**
location <- pData(BMV_metaseq_species)$Municipality
unique_locations <- length(unique(pData(BMV_metaseq_species)$Municipality))
locColors <- viridis::plasma(unique_locations, 1, begin=0, end = 0.9, direction = 1)
names(locColors) <- levels(as.factor(pData(BMV_metaseq_species)$Municipality))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(locColors, main="Location colors")

#' **Set heatmap colors:**
heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(200)
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(heatmapCols, main="Heatmap colors for abundance")

#' **Calculate average BLASTx values:**
blastx <- read.table("data/BEmosquitoes.tab", header=T, row.names=1, dec=".", sep="\t")
blastx <- dplyr::select(blastx, 2)
blastx.UF <- otu_table(as.matrix(blastx), taxa_are_rows=T)
blastx.ps <- phyloseq(blastx.UF, tax_table(prune_taxa(taxa_sums(BMV.V2) >= 1, BMV.V2)))
blastx_metaseq <- phyloseq_to_metagenomeSeq(blastx.ps)
blastx_metaseq

#Aggregate by species
blastx_mean = aggregateByTaxonomy(blastx_metaseq, lvl = "Species", norm = F, aggfun = mean, out = "MRexperiment", alternate = T)
blastx_mean

#' **Store average BLASTx identities in dataframe:**
rowanno <- as.data.frame(returnAppropriateObj(blastx_mean, log=F, norm=F))
colnames(rowanno)[1] <- "blastx"
#+ echo=FALSE
datatable(rowanno)

#' **Create color function for BLASTx values:**
col_fun=colorRamp2(c(0,100), c("white","deepskyblue4"))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(col_fun(0:100), main="Blastx % identity colors")

#' **Get Baltimore classification/genome composition for each viral species from ICTV:**
taxV <- as.data.frame(tax_table(BMV_final))
baltimore <- read.table("data/ICTV_classification_family.tsv", header=T, sep="\t")

balt_gc <- baltimore %>% 
  select(Family, Genome.Composition) %>% 
  distinct() %>% 
  add_row(Family="Picorna-like", Genome.Composition="ssRNA(+)") %>% 
  add_row(Family="Orthomyxo-like", Genome.Composition="ssRNA(-)") %>% 
  add_row(Family="Negeviridae", Genome.Composition="ssRNA(+)") %>% 
  add_row(Family="Luteoviridae", Genome.Composition="ssRNA(+)")

tax_gc <- left_join(taxV, balt_gc, by="Family") %>%
  filter(!Genome.Composition %in% c("ssRNA(+/-)", "ssDNA(+/-)", "ssDNA(+)", "ssRNA")) %>% # remove double genome compositions
  mutate(Genome.Composition = case_when(Species == 'Sclerotinia sclerotiorum dsRNA virus 3' ~ 'dsRNA',
                                        Family == 'Reo-like'~'dsRNA',
                                        Species == 'unclassified' ~'unknown',
                                        Family == 'Totiviridae'~'dsRNA',
                                        Species == 'Hubei orthoptera virus 4' ~ 'ssRNA(+)',
                                        Species == 'Hubei odonate virus 15' ~ 'dsRNA',
                                        Species == 'Wuhan insect virus 26' ~ 'dsRNA',
                                        Species == 'Wuhan insect virus 13' ~ 'ssRNA(+)',
                                        Species == "Linepithema humile qinvirus-like virus 1" ~'ssRNA(-)',
                                        Species =="Hubei virga-like virus 23" ~ 'ssRNA(+)',
                                        Species == "Culex mononega like virus 2" ~'ssRNA(-)',
                                        TRUE ~ `Genome.Composition`))

gcanno<-tax_gc %>% select(Genome.Composition, Family)
rownames(gcanno)<-tax_gc$Species
colnames(gcanno) <- c('GC', 'Family')
gcanno<-gcanno[order(rownames(gcanno)) , ,drop=F]
gcanno$Family <- factor(gcanno$Family, levels = FamLevel)
#+ echo=FALSE
datatable(gcanno)

#' ## Alpha diversity
alpha <- plot_richness(BMV_final, measures=c("Observed", "Shannon", "Simpson"), x="SKA_Subspecies", color = "SKA_Subspecies")+
  geom_boxplot()+
  geom_jitter(width=0)+
  theme_bw()+
  scale_color_viridis_d(begin=0.4, end =0.9, name="", labels=c(expression(paste(italic("Aedes japonicus"), " (n=8)")),
                                                               expression(paste(italic("Culex pipiens molestus"), " (n=23)")),
                                                               expression(paste(italic("Culex pipiens pipiens"), " (n=48)")),
                                                               expression(paste(italic("Culex torrentium"), " (n=3)"))))+
  theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0)+
  scale_y_continuous(limits = c(0, NA))+
  stat_compare_means(method = "wilcox.test", aes_string(x="SKA_Subspecies", y="value",
                                                        group="SKA_Subspecies"), 
                     label = "p.format", hide.ns=TRUE, size=3, show.legend = F, tip.length = 0.01,
                     comparisons = list(c(1, 2),c(1, 3)))
alpha
ggsave("figures/alpha-diversity-species.pdf", alpha, dpi=300)

#' ## Ordination
#' ### PCoA
v.ord <- ordinate(BMV_final, method = "PCoA")
pcoa <- plot_ordination(BMV_final, v.ord, type="samples", color="SKA_Subspecies", shape = "Municipality")+
  scale_shape_manual(values = c(15,16,17,3:8), guide =guide_legend(label.theme = element_text(size=10)), name="Location")+
  scale_color_viridis_d(begin=0, end =1, name="Mosquito species")+
  #stat_ellipse(type = "norm", linetype = 2, aes_string(group="SKA_Subspecies"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  guides(col = guide_legend(override.aes = list(shape = 15, size = 5),
                            label.theme = element_text(size=10, face="italic")))
pcoa
ggsave("figures/PCoA-eukaryotic-virome.pdf", pcoa, dpi=300)
?ggsave
#' ### NMDS full dataset
v.ord <- ordinate(BMV_final, method = "NMDS", k=2)
nmds <- plot_ordination(BMV_final, v.ord, type="samples", color="SKA_Subspecies", shape = "Municipality")+
  scale_shape_manual(values = c(15,16,17,3:8), guide =guide_legend(label.theme = element_text(size=10)), name="Location")+
  scale_color_viridis_d(begin=0, end =1, name="Mosquito species")+
  #stat_ellipse(type = "t", linetype = 2, aes_string(color="SKA_Subspecies"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  guides(col = guide_legend(override.aes = list(shape = 15, size = 5),
                            label.theme = element_text(size=10, face="italic")))
nmds
ggsave("figures/NMDS-eukaryotic-virome.pdf", nmds, dpi=300)

#' ### NMDS singletons removed
BMV_no_singletons <- BMV_final
for (i in c("MEMO014","MEMO078","MEMO145","MEMO146","MEMO149","MEMO084","MEMO018","MEMO006",
            "NEMO27","NEMO30","NEMO07","NEMO46","NEMO51","MEMO127","MEMO008")) {
  BMV_no_singletons <- subset_samples(BMV_no_singletons, Sample != i)
}
BMV_no_singletons

v.ord <- ordinate(BMV_no_singletons, method = "NMDS", k=2)
nmds_ns <- plot_ordination(BMV_no_singletons, v.ord, type="samples", color="SKA_Subspecies", shape = "Municipality")+
  scale_shape_manual(values = c(15,16,17,3:8), guide =guide_legend(label.theme = element_text(size=10)), name="Location")+
  scale_color_viridis_d(begin=0, end =1, name="Mosquito species")+
  #stat_ellipse(type = "t", linetype = 2, aes_string(color="SKA_Subspecies"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  guides(col = guide_legend(override.aes = list(shape = 15, size = 5),
                            label.theme = element_text(size=10, face="italic")))
nmds_ns
ggsave("figures/NMDS-eukaryotic-virome-no-singletons.pdf", nmds_ns, dpi=300)

ord.legend <- cowplot::get_legend(pcoa)
alpha.legend <- cowplot::get_legend(alpha)
plots <- cowplot::align_plots(alpha+theme(legend.position = "none"), 
                              pcoa+theme(legend.position = "none"), align = 'v', axis = 'l')
bottom<-cowplot::plot_grid(plots[[2]], 
                           nmds+theme(legend.position = "none"), labels=c("B","C"))
cowplot::plot_grid(plots[[1]], alpha.legend, bottom, ord.legend, labels = c('A', ''), ncol=2, rel_widths = c(1, .3))

#ggarrange(alpha, ggarrange(pcoa, nmds, labels=c("B","C"), common.legend = T, legend = "right"), labels="A", ncol=1)
ggsave("figures/combined-alpha-ordination.pdf", dpi=300)

#' ## Relative abundance
#' ### Relative abundance per location
phylobar_location <- merge_samples(BMV_family, "Municipality")
phylobar_location_re.abund <- phylobar_location %>%                   
  transform_sample_counts(function(x) {round(100*(x/sum(x)),5)} )

lc <- plot_relabund(phylobar_location_re.abund, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack")+
  ylab("Relative abundance (%)")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0, size = 10, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_text(size = 12,hjust = 0.5),
        axis.ticks.y = element_blank(),
        legend.text = element_text(face="italic"))+
  scale_fill_manual(values = myFamCol, name="Viral family")+
  scale_y_continuous(expand = c(0,0))+
  coord_flip()

#+ echo=FALSE, fig.width=10
lc

#+ echo=TRUE
#' **Combine viral family barplot with mosquito species barplots from each location:**
lc2<-ggarrange(lc+theme(legend.position = "none"), ggarrange(plotlist=plist, nrow=10, heights = c(0.2, rep(1,8),0.6)), ncol=2, widths = c(1,0.1))
ggsave("figures/viralfamily_barplot_loc.pdf", lc2, height = 6.92, width = 11.2, dpi=300)

#+ echo=FALSE
lc2

#' ### Relative abundance per mosquito species
phylobar_species <- merge_samples(BMV_family, "SKA_Subspecies")
phylobar_species_re.abund <- phylobar_species %>%                   
  transform_sample_counts(function(x) {round(100*(x/sum(x)),5)} )

sp <- plot_relabund(phylobar_species_re.abund, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack")+
  ylab("Relative abundance (%)")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0, size = 10, vjust = 0.5, hjust = 0.5, face="italic"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, vjust=0.5),
        axis.title.x = element_text(size = 12,hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face="italic"),
        legend.title = element_blank())+
  scale_fill_manual(values = myFamCol, name="Viral family")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(guide = guide_axis(n.dodge=2))
ggsave("figures/viralfamily_barplot_species.pdf", sp , height = 6.92, width = 11.2, dpi=300)

#+ echo=FALSE
sp

#' **Create combined legend for viral families and mosquito species:**
#+ echo=TRUE
sp2 <- plot_relabund(phylobar_species_re.abund, fill="Family", color="Sample")+
  theme_bw()+
  theme(legend.text = element_text(face="italic"),
        legend.title = element_blank(),
        legend.position = "top")+
  scale_fill_manual(values = myFamCol, name="a")+
  scale_color_manual(values = myColors, name="b")+
  guides(fill=guide_legend(nrow = 4),
         col=guide_legend(override.aes =list(fill = myColors), direction = "vertical"))
sp_leg <- get_legend(sp2)

#' ### Combine relative abundances
#+ fig.width=15, fig.height=9
ggarrange(sp, lc2, common.legend = T, legend = "top", labels ='AUTO',
          font.label = list(size = 16), widths = c(0.5,1), legend.grob = sp_leg)
ggsave("figures/viralfamily_barplot.pdf", width = 15, height=9, dpi=300)

#' ## Eukaryotic virome heatmap
#p.fungi<-log10(pData(BMV_metaseq_species)$fungi_reads+1)
#p.fungi<-pData(BMV_metaseq_species)$fungi_reads
#p.fungi<-as.matrix(fungi2[fungi2$Sample %in% pData(BMV_metaseq_species)$Sample,2:5])
p.fungi<-log10(as.matrix(fungi2[fungi2$Sample %in% pData(BMV_metaseq_species)$Sample,2:5]))
p.fungi[p.fungi==-Inf]<-NA
colnames(p.fungi)<-as.data.frame(tax_table(BMV_fungi))$Phylum

left_ra=rowAnnotation('Family'=gcanno$Family,
                      "Blastx"= rowanno$blastx,
                      col=list('Family'=myFamCol,
                               "Blastx"=col_fun),
                      show_annotation_name=T,
                      annotation_name_side = "top",
                      annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                      annotation_name_rot = 270,
                      annotation_legend_param=list("Blastx" = list(title ="Blastx % Identity", 
                                                                   direction="horizontal",
                                                                   at=c(0,25,50,75,100)),
                                                   "Family"=list(title="Viral family", labels_gp = gpar(fontface="italic"))))

#Draw heatmap
column_ha = HeatmapAnnotation(Location=location,
                              'Mosquito species' = mosquito_species,
                              Fungi = anno_points(p.fungi, ylim = c(0, 6),
                                                   axis_param = list(side = "right"),
                                                  gp = gpar(fill = c(2:4,7), col = c(2:4,7)),
                                                  pch=1:length(colnames(p.fungi))),
                              show_annotation_name = T,
                              annotation_label = gt_render(c("Location", "Species", "log<sub>10</sub>(Fungi)")),
                              annotation_name_side = "right",
                              annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                              annotation_name_rot = 0,
                              annotation_legend_param=list('Mosquito species'= list(title ='Mosquito species', 
                                                                                    labels_gp = gpar(fontface="italic")),
                                                           'Location' = list(title='Location')),
                              col = list('Mosquito species'=c(myColors), Location=c(locColors)))

#n depends on the amount of taxa in 'BMV_metaseq_species'
hm <- plot_abundance(BMV_metaseq_species, n = n_species, log = T, norm = F, colclust = "bray",
                     col = heatmapCols, name = "Log2 Read Counts", 
                     top_annotation=column_ha,
                     row_names_gp = gpar(fontsize = 8),
                     show_column_names = FALSE,
                     heatmap_legend_param = list(direction = "horizontal"),
                     cluster_rows = FALSE,
                     cluster_row_slices = FALSE,
                     row_split=factor(gcanno$GC, 
                                      levels = c("dsDNA","ssDNA","dsRNA","ssRNA(+)","ssRNA(-)","unknown")),
                     row_order=tax_gc[with(tax_gc, order(tax_gc$Genome.Composition, Family)),"Species"],
                     left_annotation=left_ra,
                     row_title_gp = gpar(fontsize=10),
                     row_title_rot=0,
                     border=T)

#' Legends
lgd_list<-list(Legend(labels = colnames(p.fungi), title = "Fungi", type = "points", pch = 1:length(colnames(p.fungi)), 
            legend_gp = gpar(col = c(2:4,7)), background ="transparent"))
#+ echo=TRUE, fig.width=15, fig.height=9
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left", 
     merge_legend = T, legend_grouping="original",
     annotation_legend_list = lgd_list)
#+ include=FALSE
pdf("figures/BEmosquitoes_hm_baltimore.pdf", width = 16.6, height =10)
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left", 
     merge_legend = T, legend_grouping="original", 
     annotation_legend_list = lgd_list)
dev.off()

#' ***
#' # RT-qPCR analysis
#' ## Prepare the data
library(ggh4x)
qPCR <- read_csv("data/BMV-RT-qPCR-analysis.csv", col_names = T)

#' **Recalculate quantity based on mosquito dilutions:**
qPCR <- qPCR %>%  
  mutate(Quantity = case_when(str_detect(Sample, "MEMO") ~ qPCR$"Quantity Mean" * 40,
                              str_detect(Sample, "NEMO") ~ qPCR$"Quantity Mean" * 100))
#+ echo=FALSE
datatable(qPCR)

#' **Divide quantity of *Xanthi chryso-like virus* by 2 because it is a dsRNA virus:**
qPCR <- qPCR %>%  
  mutate(Quantity = case_when(str_detect(Target, "XCV") ~ qPCR$Quantity/2,
                              TRUE ~ Quantity))
#+ echo=FALSE
datatable(qPCR)

#' **Replace NA to 0:**
qPCR_noNA <- qPCR %>% mutate_all(~replace(., is.na(.), 0))

#' **Replace abbreviations for full names:**
qPCR_noNA <- qPCR_noNA %>% mutate(Target, Target= case_when(Target=="CPV"~"Culex phasma virus", 
                                                            Target=="XCV"~"Xanthi chryso-like virus", 
                                                            Target=="WMV6"~"Wuhan Mosquito Virus 6",
                                                            Target=="WMV4"~"Wuhan Mosquito Virus 4",
                                                            Target=="Daeseongdong"~"Daeseongdong virus",
                                                            Target=="HMV4"~"Hubei mosquito virus 4",
                                                            TRUE ~ Target))

#' **Merge metadata with qPCR data:**
#' Abrreviate locations to fit the names on the plots later on.
metadata <- meta %>% mutate(Municipality = case_when(Municipality == 'Frameries' ~ 'Fr',
                                                         Municipality == 'Dilsen-Stokkem' ~ 'DS',
                                                         Municipality == 'Kallo' ~ 'K',
                                                         TRUE ~ Municipality))
metaqPCR <- merge.data.frame(qPCR_noNA, metadata, by="Sample")
#+ echo=FALSE
datatable(metaqPCR)

#' **Count number of tested samples per mosquito species:**
metaqPCR %>%
  select(Sample, SKA_Subspecies) %>% 
  distinct() %>% 
  count(SKA_Subspecies)

#' ## Plot data
legendlabels <- c(expression(paste(italic("Aedes japonicus"), " (n=8)")),
  expression(paste(italic("Culex pipiens molestus"), " (n=47)")),
  expression(paste(italic("Culex pipiens pipiens"), " (n=127)")),
  expression(paste(italic("Culex torrentium"), " (n=16)")))

qPCRviolin <- ggplot(metaqPCR, aes(x=SKA_Subspecies, y=Quantity))+
  geom_violin(aes(fill=SKA_Subspecies, color=SKA_Subspecies), position="identity")+
  geom_jitter(aes(color=SKA_Subspecies), size=1, width=0.1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.text.align = 0)+
  scale_y_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.1)), 
                     trans=scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 10^seq(2,8,2)),
                     labels = scientific_10(c(0, 10^seq(2,8,2))))+
  xlab(label = "")+
  ylab(label = "Viral copies per mosquito")+
  scale_fill_viridis_d(begin=0.4, end =0.9, name="", alpha = 0.5,
                       labels=legendlabels)+
  scale_color_viridis_d(begin=0.4, end =0.9, name="",
                        labels=legendlabels)+ 
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  facet_wrap(~Target, nrow=2, ncol=3)
ggsave("figures/qPCR-violin.pdf", dpi=300)
#+ echo=FALSE
qPCRviolin

#' **Plot qPCR overview per sample:**
qPCRoverview <- ggplot(metaqPCR, 
             aes(x=Sample, y=Quantity))+
  facet_nested(Target~Municipality,
               space="free",
               scales = "free_x",
               labeller = labeller(Target = label_wrap_gen(10), Municipality=label_wrap_gen(5)))+
  geom_col(aes(fill=SKA_Subspecies), position="identity", alpha=1)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank())+
  scale_y_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.1)),
                     trans = scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 10^seq(2,8,2)),
                     labels = scientific_10(c(0, 10^seq(2,8,2))))+
  xlab(label = "Sample")+
  ylab(label = "Viral copies per mosquito")+
  scale_fill_viridis_d(begin=0.4, end =0.9, name="")+ 
  guides(fill = guide_legend(override.aes = list(alpha = 1), 
                             label.theme = element_text(size=10, angle = 0, face = "italic")))
ggsave("figures/qPCR_overview.pdf", width = 18.7, height = 10.3, units = "in", dpi=300)
#+ echo=FALSE, fig.width=18.7, fig.height=10.3, fig.cap="DS=Dilsen-Stokkem, Fr=Frameries, K=Kallo"
qPCRoverview

#' ***
#' # Phylogenetics
#+ results="hide"
phylo_path <- here("phylogenetics")
phylo_packages <- c('phytools', 'ggtree', 'treeio', 'rphylopic',
                    'ggimage', 'aplot')
lapply(phylo_packages, library, character.only = TRUE)

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

bunya@phylo$tip.label <- gsub(bunya@phylo$tip.label, pattern = "_", replacement=" ")
bunya@phylo$tip.label <- gsub(bunya@phylo$tip.label, pattern = "P ", replacement="P_")
bunya@phylo$tip.label <- gsub(bunya@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

bunya_clean <- clean_phylometa(bunya, metadata = meta.bunya)
bunya_ictv_clean <- clean_phylometa(bunya, metadata = bunya.ictv)

bunya_tree <- group_phylotaxa(bunya, bunya_clean, group = c("Family","Genus", "Host", "Geo_Location"))
bunya_tree

#' Give clean metadata ictv!
bunyacladedf <- mrca_ictv(tree = bunya_tree, meta = bunya_ictv_clean, 
                          group="Genus")

famvec <- sort(unique(bunya_clean$Family))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.4, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

bunya_tree@phylo$tip.label <- gsub(bunya_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)
                   
p <- plot_phylotree(bunya_tree, col=col, shape=shape, cladedf = bunyacladedf, 
                    labels=c(famvec, "Belgian mosquitoes"),
               plainlabels = c("unclassified", "Belgian mosquitoes"), align=F)+
  geom_nodepoint(aes(subset= node %in% bunyacladedf$node), size=1, shape=18)+
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

endorna@phylo$tip.label <- gsub(endorna@phylo$tip.label, pattern = "_", replacement=" ")
endorna@phylo$tip.label <- gsub(endorna@phylo$tip.label, pattern = "P ", replacement="P_")
endorna@phylo$tip.label <- gsub(endorna@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

endorna_clean <- clean_phylometa(endorna, metadata = meta.endorna)
endorna_ictv_clean <- clean_phylometa(endorna, metadata = endorna.ictv)

endorna_tree <- group_phylotaxa(endorna, endorna_clean, group = c("Family","Genus", "Host", "Geo_Location"))
endorna_tree

#' Give clean metadata ictv!
endornacladedf <- mrca_ictv(tree= endorna_tree, meta = endorna_ictv_clean,
                            group = "Genus")

famvec <- sort(unique(endorna_clean$Genus))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.7, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

endorna_tree@phylo$tip.label <- gsub(endorna_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(endorna_tree, col=col, shape=shape, cladedf = endornacladedf, labels=c(famvec, "Belgian mosquitoes"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F, group="Genus")+
  geom_nodepoint(aes(subset= node %in% endornacladedf$node), size=1, shape=18)+
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

ghabri@phylo$tip.label <- gsub(ghabri@phylo$tip.label, pattern = "_", replacement=" ")
ghabri@phylo$tip.label <- gsub(ghabri@phylo$tip.label, pattern = "P ", replacement="P_")
ghabri@phylo$tip.label <- gsub(ghabri@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

ghabri_clean <- clean_phylometa(ghabri, metadata = meta.ghabri)
ghabri_ictv_clean <- clean_phylometa(ghabri, metadata = ghabri.ictv)

ghabri_tree <- group_phylotaxa(ghabri, ghabri_clean, group = c("Family","Genus", "Host", "Geo_Location"))
ghabri_tree

#' Give clean metadata ictv!
ghabricladedf <- mrca_ictv(tree=ghabri_tree, meta =ghabri_ictv_clean, 
                          group="Genus", subset = c("unclassified", "Chrysovirus"))

famvec <- sort(unique(ghabri_clean$Family))
famvec
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.2, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

ghabri_tree@phylo$tip.label <- gsub(ghabri_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(ghabri_tree, col=col, shape=shape, cladedf = ghabricladedf, labels=c(famvec, "Belgian mosquitoes"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F)+
  geom_nodepoint(aes(subset= node %in% ghabricladedf$node), size=1, shape=18)+
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
#' Rerun tree with hmv4!
meta.noda <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Nodaviridae/"), pattern = "*.csv", full.names = T), read_csv))
noda.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Nodaviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
noda <- read.iqtree(here(phylo_path, "Nodaviridae/nodamuvirales.treefile"))
noda@phylo <- phytools::midpoint.root(noda@phylo)

noda@phylo$tip.label <- gsub(noda@phylo$tip.label, pattern = "_", replacement=" ")
noda@phylo$tip.label <- gsub(noda@phylo$tip.label, pattern = "P ", replacement="P_")
noda@phylo$tip.label <- gsub(noda@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

noda_clean <- clean_phylometa(noda, metadata = meta.noda)
noda_ictv_clean <- clean_phylometa(noda, metadata = noda.ictv)

noda_tree <- group_phylotaxa(noda, noda_clean, group = c("Family","Genus", "Host", "Geo_Location"))
noda_tree

#' Give clean metadata ictv!
nodacladedf <- mrca_ictv(tree = noda_tree, meta = noda_ictv_clean, group = "Genus")

famvec <- sort(unique(noda_clean$Family))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.4, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

noda_tree@phylo$tip.label <- gsub(noda_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(noda_tree, col=col, shape=shape, cladedf = nodacladedf, labels=c(famvec, "Belgian mosquitoes"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F)+
  geom_nodepoint(aes(subset= node %in% nodacladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,4)

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

orthomyxo@phylo$tip.label <- gsub(orthomyxo@phylo$tip.label, pattern = "_", replacement=" ")
orthomyxo@phylo$tip.label <- gsub(orthomyxo@phylo$tip.label, pattern = "P ", replacement="P_")
orthomyxo@phylo$tip.label <- gsub(orthomyxo@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

orthomyxo_clean <- clean_phylometa(orthomyxo, metadata = meta.orthomyxo)
orthomyxo_ictv_clean <- clean_phylometa(orthomyxo, metadata = orthomyxo.ictv)

orthomyxo_tree <- group_phylotaxa(orthomyxo, orthomyxo_clean, group = c("Family","Genus", "Host", "Geo_Location"))
orthomyxo_tree

#' Give clean metadata ictv!
orthocladedf <- mrca_ictv(tree=orthomyxo_tree, meta = orthomyxo_ictv_clean, 
                          group="Genus")

famvec <- sort(unique(orthomyxo_clean$Genus))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

orthomyxo_tree@phylo$tip.label <- gsub(orthomyxo_tree@phylo$tip.label, pattern = "lcl [ORF]*[0-9]* *", replacement="", perl = T)

p <- plot_phylotree(orthomyxo_tree, col=col, shape=shape, cladedf = orthocladedf, labels=c(famvec, "Belgian mosquitoes"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F, group="Genus")+
  geom_nodepoint(aes(subset= node %in% orthocladedf$node), size=1, shape=18)+
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

picorna@phylo$tip.label <- gsub(picorna@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

picorna_clean <- clean_phylometa(picorna, metadata = meta.picorna)
picorna_ictv_clean <- clean_phylometa(picorna, metadata = picorna.ictv)

picorna_tree <- group_phylotaxa(picorna, picorna_clean, 
                                group = c("Family","Genus", "Host", "Geo_Location"))

#' Give clean metadata ictv!
picornacladedf <- mrca_ictv(tree=picorna_tree, meta = picorna_ictv_clean, 
                          group="Genus", subset="unclassified")#, subset = "Cripavirus")

famvec <- sort(unique(picorna_clean$Family))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.2, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

picorna_tree@phylo$tip.label <- gsub(picorna_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)


p <- plot_phylotree(picorna_tree, col=col, shape=shape, cladedf = picornacladedf, labels=c(famvec, "Belgian mosquitoes"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F)+
  geom_nodepoint(aes(subset= node %in% picornacladedf$node), size=1, shape=18)+
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

#' ## Reoviridae

meta.reo <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Reoviridae/"), pattern = "*.csv", full.names = T), read_csv))
reo.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Reoviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
reo <- read.iqtree(here(phylo_path, "Reoviridae/reoviridae.treefile"))
reo@phylo <- phytools::midpoint.root(reo@phylo)

reo@phylo$tip.label <- gsub(reo@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

reo_clean <- clean_phylometa(reo, metadata = meta.reo)
reo_ictv_clean <- clean_phylometa(reo, metadata = reo.ictv)

reo_tree <- group_phylotaxa(reo, reo_clean, group = c("Family","Genus", "Host", "Geo_Location"))
reo_tree@phylo$tip.label

#' Give clean metadata ictv!
reocladedf <- mrca_ictv(tree=reo_tree, meta =reo_ictv_clean, 
                        group="Genus")

famvec <- sort(unique(reo_clean$Genus))
col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

reo_tree@phylo$tip.label <- gsub(reo_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)
#+echo=FALSE
reo_tree@phylo$tip.label <- gsub(reo_tree@phylo$tip.label, pattern = " pasted", replacement="", perl = T)
#+echo=TRUE
p <- plot_phylotree(reo_tree, col=col, shape=shape, cladedf = reocladedf, labels=c(famvec, "Belgian mosquitoes"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F, group="Genus")+
  geom_nodepoint(aes(subset = node %in% reocladedf$node), size=1, shape=18)+
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

meta.tombus <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Tombusviridae/"), pattern = "*.csv", full.names = T), read_csv))
tombus.ictv <- do.call(rbind, lapply(list.files(path=here(phylo_path, "Tombusviridae/"), pattern = "*_ictv.csv", full.names = T), read_csv))
tombus <- read.iqtree(here(phylo_path, "Tombusviridae/tombusviridae2.treefile"))
tombus@phylo <- phytools::midpoint.root(tombus@phylo)

tombus@phylo$tip.label <- gsub(tombus@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

tombus_clean <- clean_phylometa(tombus, metadata = meta.tombus)
tombus_ictv_clean <- clean_phylometa(tombus, metadata = tombus.ictv)

tombus_tree <- group_phylotaxa(tombus, tombus_clean, 
                               group = c("Family","Genus", "Host", "Geo_Location"))

#' Give clean metadata ictv!
tombuscladedf <- mrca_ictv(tree = tombus_tree, meta = tombus_ictv_clean, 
                          group="Genus")

famvec <- sort(unique(tombus_clean$Genus))

col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

tombus_tree@phylo$tip.label <- gsub(tombus_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)


p <- plot_phylotree(tombus_tree, col=col, shape=shape, cladedf = tombuscladedf,
                    labels=c(famvec, "Belgian mosquitoes"), 
                    plainlabels = c("unclassified", "Belgian mosquitoes"), align=F,
                    group="Genus")+
  geom_nodepoint(aes(subset= node %in% tombuscladedf$node), size=1, shape=18)+
  scale_y_reverse()+
  xlim(-.05,3.5)
tombusplot <- addSmallLegend(p, pointSize = 2.5, textSize = 5, spaceLegend = .5)+
  theme(legend.text = element_text(vjust = .4),
        legend.position = c(0.1,.2))
tombusplot

ggsave(filename = here(phylo_path, "Tombusviridae/tombusviridae.pdf"), tombusplot, dpi=1200, width=6.875)

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

negev@phylo$tip.label <- gsub(negev@phylo$tip.label, pattern = "lcl\\|", replacement="lcl_")

negev_clean <- clean_phylometa_negev(negev, metadata = meta.negev)

#negev_ictv_clean <- clean_phylometa(negev, metadata = negev.ictv)

negev_tree <- group_phylotaxa(negev, negev_clean, 
                               group = c("Family","Genus", "Host", "Geo_Location"))

#' Give clean metadata ictv!
# negevcladedf <- mrca_ictv(tree = negev_tree, meta = negev_ictv_clean, 
#                            group="Genus")

famvec <- sort(unique(negev_clean$Genus))

col <- c(viridis::plasma(n = length(famvec)-1, begin=0.1, end=0.9))
names(col) <- famvec[famvec!="unclassified"]
col <- c(col, unclassified="grey", NODE="#43BF71FF")
shape <- rep(16, length(famvec))
names(shape) <- famvec
shape <- c(shape, NODE=17)

negev_tree@phylo$tip.label <- gsub(negev_tree@phylo$tip.label, pattern = "lcl ORF[0-9]* ", replacement="", perl = T)

p <- plot_phylotree(negev_tree, col=col, shape=shape,
                    labels=c(famvec, "Belgian mosquitoes"), 
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

ggsave(filename = here(phylo_path, "Negevirus/negeviruses_minimal.pdf"), tombusplot_minimal, dpi=1200, width=6.875)


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
                picornaplot, bunyaplot, 
                ghabriplot, #reoplot,
                #tombusplot, negevplot,
                ncol=3, nrow=2, labels="AUTO",
                heights = c(.4,.6))

#+ echo=TRUE, fig.width=20, fig.height=20
# ggarrange(p1, 
#           ggarrange(NULL, tombusplot+theme(legend.position = "right"), NULL, 
#                     widths=c(0.1,0.8,0.1), labels=c("","G",""), ncol=3), 
#           ncol=1, heights = c(.75, .25))
ggarrange(p1, 
          ggarrange(reoplot, tombusplot,#+theme(legend.position = "right"), 
                    widths=c(0.5,0.5), labels=c("G","H"), ncol=2), 
          ncol=1, heights = c(.7, .3))

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
ggarrange(p2, 
          ggarrange(reoplot_minimal, tombusplot_minimal, 
                    widths=c(0.5,0.5), labels=c("G","H"), ncol=2), 
          ncol=1, heights = c(.7, .3))

#ggsave("figures/treeplots_minimal.pdf", width=17, height = 28, units = "in", dpi=300)
ggsave("figures/treeplots_minimal.pdf", width=20, height = 20, units = "in", dpi=600)
ggsave("figures/treeplots_minimal_small.pdf", width=6.875, height = 9.0625, units = "in", dpi=600)

#' ***
#' # Phageome & Wolbachia bacteria
library(ggbeeswarm)
phage <- read.table("data/wolbachia_phage_ratio.tsv", sep="\t", header = T)
phage <- merge(phage, meta, by="Sample")
ggplot(phage, aes(x = SKA_Subspecies, y= ratio, color=ratio>4))+
  geom_beeswarm(cex=2)+
  scale_color_brewer(palette = "Set2", direction = -1)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face="italic"))

phageWO_real <- phage[phage$ratio>4,]$Sample


#' ## Wolbachia
#' Wolbachia bwa-mem2 mapped
wMel <- read.table("data/wMel.txt", header=F, dec=".", sep=" ")
wPip <- read.table("data/wPip.txt", header=F, dec=".", sep=" ")
phageWO <- read.table("data/phageWO.txt", header=F, dec=".", sep=" ")
colnames(wMel) <- c("Sample", "numreads", "covbases", "coverage", "meandepth")
colnames(wPip) <- c("Sample", "numreads", "covbases", "coverage", "meandepth")
colnames(phageWO) <- c("Sample", "supergroup", "numreads", "covbases", "coverage", "meandepth")
wPip$supergroup <- "wPip"
wMel$supergroup <- "wMel"

wolb <- rbind(wMel, wPip, phageWO)

wolb <- wolb %>% 
  plyr::mutate(supergroup=str_replace(string=supergroup, pattern = "NODE_12_length_11674_cov_94.907217_MEMO129", replacement = "phageWO2")) %>% 
  plyr::mutate(supergroup=str_replace(string=supergroup, pattern = "NODE_1_length_31159_cov_25.624799_MEMO050", replacement = "phageWO1"))
wolb$numreads[wolb$coverage<5] <- 0
wolb$meandepth[wolb$coverage<5] <- 0

wviolin <- pivot_wider(wolb, id_cols = "Sample", names_from = "supergroup", values_from = "numreads") %>% 
  ggplot(aes(x = meta[meta$Sample!="NEMO39",]$SKA_Subspecies, y=wPip))+
  geom_boxplot(aes(color=meta[meta$Sample!="NEMO39",]$SKA_Subspecies), width=0.1, outlier.shape = NA)+
  geom_violin(aes(fill=meta[meta$Sample!="NEMO39",]$SKA_Subspecies), color=NA, show.legend = F, trim = T)+
  geom_jitter(aes(color=meta[meta$Sample!="NEMO39",]$SKA_Subspecies), size=1, alpha=0.6,
              position=position_jitter(w=0.1,h=0)) +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank())+
  xlab(label = "")+
  ylab(label = "Wolbachia read count")+
  scale_fill_manual(values=alpha(myColors, 0.1), name='')+
  scale_color_viridis_d(begin=0.4, end =0.9, name="")+ 
  guides(fill = guide_legend(#override.aes = list(alpha = 1), 
    label.theme = element_text(size=10, angle = 0, face = "italic")))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 10^seq(2,8,2)))
wviolin
ggsave("figures/wolbachia_readcount.pdf", dpi=300)


#' barplot
wb <-  pivot_wider(wolb, id_cols = "Sample", names_from = "supergroup", values_from = "numreads") %>% 
  ggplot(aes(x = Sample, y=wPip))+
  geom_col(aes(fill=meta[meta$Sample!="NEMO39",]$SKA_Subspecies))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank())+
  scale_y_continuous(expand=expansion(mult = c(0, 0.01), add = c(0, 0.1)))+
  xlab(label = "Sample")+
  ylab(label = "Wolbachia read count")+
  labs(shape="Virus", color="Virus")+
  scale_fill_viridis_d(begin=0.4, end =0.9, name="")+ 
  guides(fill = guide_legend(override.aes = list(alpha = 1), 
                             label.theme = element_text(size=10, angle = 0, face = "italic")))+
  facet_grid(~meta[meta$Sample!="NEMO39",]$Municipality, scales = "free", space = "free")
wb


#' ## Phageome
#' 
phageotu<-read.table("data/BEmosq_classification-85-taxfile-1000nt.tsv", sep="\t", header = T)
phage<-phageotu[phageotu$Phylum == "Phage",]
rownames(phage)<- phage$Contig
phage <- subset(phage, select=-Contig)

newtax<-tax[!rownames(tax) %in% rownames(phage),]

newtax <- rbind(newtax, phage)

tax_table(BMV) <- as.matrix(newtax) #Maybe better solution

#newtax.UF <- tax_table(as.matrix(newtax))

#BMV <- phyloseq(OTU.UF, newtax.UF, meta.UF)

BMV.V <- subset_taxa(BMV, Kingdom=="Viruses")

phageome <- subset_taxa(BMV.V, Phylum=="Phage" | Order=="Caudovirales" | Family=="Inoviridae")
phageome

otu_table(phageome) <- as.data.frame(otu_table(phageome)) %>% 
  mutate(otu=rownames(.)) %>% 
  pivot_longer(cols=1:197, names_to = "Sample", values_to = "abundance") %>% 
  mutate(abundance=if_else(otu %in% c("NODE_1_length_31159_cov_25.624799_MEMO050|full", 
                                            "NODE_12_length_11674_cov_94.907217_MEMO129|full") 
                                 & Sample != phageWO_real, 0, as.double(abundance))) %>% 
  pivot_wider(names_from = Sample, values_from = abundance) %>% 
  column_to_rownames('otu') %>% 
  otu_table(taxa_are_rows = T)

#sample_data(phageome) <- data.frame(wv.meta, row.names = 1) 
#sample_data(phageome)

phage_pruned <- prune_samples(sample_sums(phageome)>0, phageome)
phage_pruned <- prune_taxa(taxa_sums(phage_pruned)>0, phage_pruned)
phage_pruned

pp <- as.data.frame(otu_table(phage_pruned))

p.merge <- merge_taxa(phage_pruned, 1:2, 2)
p.merge <- merge_taxa(p.merge, 2:3, 2)
p.merge <- merge_taxa(p.merge, 3:25, 2)
p.merge <- merge_taxa(p.merge, 4:5, 2)
tax_table(p.merge)
otu_table(p.merge)
plot_heatmap(p.merge, taxa.label = "Species", low ="#FFF7B8", high = "#800026", na.value = "#FFFFCC")+
  theme_bw()+
  theme(axis.title = element_blank())


phage_tax <- as.data.frame(tax_table(p.merge))
phage_tax <- replace_na(phage_tax, list(Phylum="Uroviricota", Class="Caudoviricetes", Order="Caudovirales", 
                                        Family="unclassified", Genus="unclassified", Species="unclassified"))
tax_table(p.merge) <- as.matrix(phage_tax)

#+ include=FALSE
spin("BelgianMosquitoVirome.R", knit=F, format = "Rmd")
#' <div class="tocify-extend-page" data-unique="tocify-extend-page" style="height: 0;"></div>
