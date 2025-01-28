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
knitr::opts_chunk$set(
  warning = FALSE, message = FALSE,
  fig.width = 11.2, fig.height = 6.92, cache = F
) # , fig.path='figures/')
#+ include=TRUE, results="hide"
necessary_packages <- c(
  "tidyverse", "data.table", "scales", "knitr", "grid", "magrittr", "DT", "here", "gridtext",
  "ggpubr", "ggthemes", "reshape2", "viridis", "pals", "circlize", "scatterpie", "ggh4x",
  "vegan", "phyloseq", "metagenomeSeq", "ComplexHeatmap", "decontam"
)
lapply(necessary_packages, library, character.only = TRUE)
i_am("BelgianMosquitoVirome.R")
setwd(here::here())
source("BMV_functions.R")
#+ echo=FALSE
si <- sessioninfo::session_info()
pckgs <- map2(
  si$packages$package[si$packages$attached == T],
  si$packages$loadedversion[si$packages$attached == T],
  ~ paste0(.x, " ", .y)
) %>%
  simplify()
print(si$platform$version)
print(paste0("Running under: ", si$platform$os))
print(pckgs)

#' ***
#' # Eukaryotic virome analysis
#' ## General sequencing info
seqnum_raw <- read.delim("data/sequencing-info.tsv", sep = "\t", header = T)
seqnum_raw_summarized <- seqnum_raw %>%
  group_by(sample) %>%
  summarise(sumseqs = sum(num_seqs))


seqnum <- read.delim("data/sequencing-info-hostout.tsv", sep = "\t", header = T)
seqnum_summarized <- seqnum %>%
  group_by(sample) %>%
  summarise(sumseqs_nonhost = sum(num_seqs))

seqnum_final <- merge(seqnum_raw_summarized, seqnum_summarized)
seqnum_final <- seqnum_final %>%
  mutate(proportion_nonhost = sumseqs_nonhost / sumseqs * 100)

datatable(seqnum_final)

max(seqnum_final$sumseqs_nonhost)
min(seqnum_final$sumseqs_nonhost)
mean(seqnum_final$sumseqs_nonhost)

max(seqnum_final$proportion_nonhost)
min(seqnum_final$proportion_nonhost)
mean(seqnum_final$proportion_nonhost)

sum(seqnum_final$sumseqs)
sum(seqnum_final$sumseqs_nonhost)

#' ## Prepare the data
#' ### Load the OTU table, taxonomy file and metadata into R
OTU <- read.table("data/abundance-table.txt", header = TRUE, row.names = 1, sep = "\t", dec = ".")
names(OTU) <- gsub(x = names(OTU), pattern = "\\.", replacement = "-")
summary(rowSums(OTU))
summary(colSums(OTU))

tax <- read.table("data/BEmosq_classification-1000nt.tsv", header = TRUE, row.names = 1, sep = "\t", dec = ".")
meta <- read.table("data/BEmosq_metadata.csv", header = TRUE, row.names = 1, sep = ";", dec = ".")
meta <- cbind(Sample = rownames(meta), meta)

meta <- meta %>%
  mutate(Trap = case_when(
    substr(Original_name, 4, 4) == "G" ~ "Gravid",
    substr(Original_name, 4, 4) == "M" ~ "Mosquito Magnet",
    substr(Original_name, 4, 4) == "B" ~ "BG Sentinel",
    Location == "control" ~ NA,
    Municipality %in% c("Leuven", "Bertem") ~ "BG Sentinel"
  ))

#+ echo=FALSE
datatable(meta)

#' ### Make a phyloseq object
OTU.UF <- otu_table(as.matrix(OTU), taxa_are_rows = T)
tax.UF <- tax_table(as.matrix(tax))
meta.UF <- sample_data(meta)

BMV <- phyloseq(OTU.UF, tax.UF, meta.UF)
BMV

#' ### Remove contamination
#' #### Visualize library sizes of samples and negative controls
decontam <- as.data.frame(sample_data(BMV))
decontam$LibrarySize <- sample_sums(BMV)
decontam <- decontam[order(decontam$LibrarySize), ]
decontam$Index <- seq(nrow(decontam))
ggplot(data = decontam, aes(x = Index, y = LibrarySize, color = Control)) +
  geom_point() +
  ggtitle("Library sizes")

#' #### Detect contaminants
sample_data(BMV)$is.neg <- sample_data(BMV)$Control == "Yes"
contamdf.prev <- isContaminant(BMV, method = "prevalence", neg = "is.neg")

#' **Number of contaminating contigs:**
table(contamdf.prev$contaminant)

#' #### Visualize prevalence of contaminants in samples and negative controls
ps.pa <- transform_sample_counts(BMV, function(abund) 1 * (abund > 0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "Yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "No", ps.pa)

df.pa <- data.frame(
  pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg),
  contaminant = contamdf.prev$contaminant
)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (True Samples)") +
  ggtitle("Prevalence of contaminants in samples vs NCs")

#' #### Remove negative controls and contaminants from phyloseq object
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, BMV)
ps.noncontam
ps.noncontam <- prune_samples(sample_data(BMV)$Control != "Yes", ps.noncontam)
BMV <- ps.noncontam
BMV

#' Remove negative controls from metadata table
meta <- meta[meta$SKA_Subspecies != "control", ]

#' ### Subset virome data
#' #### Subset only viruses and remove prokaryotic viruses and EVEs
BMV.V <- subset_taxa(BMV, Kingdom == "Viruses")
BMV.V

EVE_phage <- c(
  "Atrato Retro-like virus", "Gurupi chuvirus-like 1", "Aedes aegypti To virus 1",
  "Aedes aegypti To virus 2", "Guato virus", "Kaiowa virus", "Atrato Chu-like virus 1",
  "Chuvirus Mos8Chu0", "Chibugado virus", "Prokaryotic dsDNA virus sp."
)

BMV.V2 <- subset_taxa(BMV.V, !is.element(Species, EVE_phage))
BMV.V2 <- subset_taxa(BMV.V2, Order != "Caudovirales")
BMV.V2 <- subset_taxa(BMV.V2, Phylum != "Phage")
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
  summarise(n = n()) %>%
  ggplot(aes(x = viral_count, y = n)) +
  geom_col() +
  theme_bw()

df <- merge(meta, vcount, by = 0) %>%
  group_by(viral_count, SKA_Subspecies) %>%
  summarise(n = n())

p <- ggplot(df, aes(x = viral_count, y = n, label = n, fill = SKA_Subspecies)) +
  geom_col() +
  scale_x_continuous(n.breaks = max(df$viral_count), name = "# viral species") +
  scale_y_continuous(breaks = seq(0, 125, 25), name = "# mosquitoes") +
  scale_fill_viridis_d(begin = 0.4, end = 0.9, name = "") +
  # geom_text(position =position_stack(vjust=.5), size=3)+
  theme_bw() +
  theme(
    legend.text = element_text(face = "italic"),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p

#+ eval=FALSE,echo=FALSE
ggarrange(ggarrange(NULL, p + theme(legend.position = "none"), ncol = 1, labels = c("A", "C")), NULL,
  ncol = 2,
  labels = c("", "B"), widths = c(.4, .6)
)

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
  transform_sample_counts(function(x) {
    round(100 * (x / sum(x)), 5)
  })
BMV_family_re.abund

#' **Count mosquito species per location:**
sp_count <- psmelt(BMV_family_re.abund) %>%
  select(Sample, SKA_Subspecies, Municipality) %>%
  distinct() %>%
  group_by(Municipality) %>%
  count(SKA_Subspecies) %>%
  pivot_wider(names_from = Municipality, values_from = n) %>%
  replace(is.na(.), 0) %>%
  rename("Villers_Le_Bouillet" = `Villers-Le-Bouillet`)
#+ echo=FALSE
datatable(sp_count)

#' **Create list of locations:**
loc_order <- c(
  "Natoye", "Vrasene", "Villers_Le_Bouillet",
  "Frameries", "Maasmechelen", "Bertem",
  "Leuven", "Eupen"
)

#' **Create list of barplots with number of mosquitoes per location:**
plist <- list(ggplot() +
  theme_void() + # ggtitle("Mosquito species")+
  theme(title = element_text(size = 6)))
for (i in loc_order) {
  plist[[i]] <- ggplot(sp_count, aes_string(x = "SKA_Subspecies", y = i, fill = "SKA_Subspecies")) +
    geom_col() +
    # scale_color_viridis_d(begin=0.4, end=0.9)+
    scale_fill_viridis_d(begin = 0.4, end = 0.9, name = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      axis.line = element_line(colour = "black", size = 0.25),
      axis.ticks = element_line(colour = "black", size = 0.25)
    ) +
    # ggtitle(gsub("_","-",i))+
    scale_y_continuous(limits = c(0, 20), expand = c(0, 0), position = "right")
}
plist[10] <- list(ggplot() +
  theme_void())
ggarrange(
  plotlist = plist[2:9], labels = gsub("_", "-", loc_order),
  font.label = list(size = 10, face = "plain"), hjust = 0, vjust = 1
)

#' **Create colorpalette for viral families:**
BMV_smelt <- psmelt(BMV_final)
FamLevel <- levels(as.factor(BMV_smelt[(BMV_smelt$Family != "unclassified"), ]$Family))
FamLevel <- c(FamLevel, "unclassified")
BMV_smelt$Family <- factor(BMV_smelt$Family, levels = FamLevel)

myFamCol <- c(stepped(n = 20), "#E6550D", "#FD8D3C", "#FDAE6B", "lightgrey")
names(myFamCol) <- levels(as.factor(BMV_smelt$Family))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(myFamCol, main = "Viral families")

#+ echo=TRUE
#' **Create colors for each mosquito species:**
# myColors <- viridis::viridis(4,1, begin=0.4, end = 0.9, direction = 1)
myColors <- viridis::viridis(4, 1, begin = 0, end = .9, direction = 1)
names(myColors) <- levels(as.factor(BMV_smelt$SKA_Subspecies))

#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(myColors, main = "Mosquito species")

#' ### Prepare heatmap data
#' **Make table with fungi reads**:
BMV_fungi <- subset_taxa(BMV, Phylum %in% c(
  "Ascomycota", "Basidiomycota",
  "Chytridiomycota", "Glomeromycota",
  "Zygomycota", "Neocallimastigomycota",
  "Microsporidia"
))
BMV_fungi

BMV_fungi <- tax_glom(BMV_fungi, taxrank = "Phylum")

sample_variables(BMV_fungi)
fungi <- as.data.frame(otu_table(BMV_fungi))
fungi2 <- fungi %>%
  mutate(Phylum = rownames(.), .before = 1) %>%
  pivot_longer(cols = 2:198, names_to = "Sample", values_to = "fungi_reads") %>%
  pivot_wider(names_from = Phylum, values_from = fungi_reads)
fungi2

tax_table(BMV_fungi)

#' **Convert phyloseq object to metagenomeseq object to make heatmap:**
BMV_metaseq <- phyloseq_to_metagenomeSeq(BMV_final)

#' **Aggregate by species:**
BMV_metaseq_species <- aggregateByTaxonomy(BMV_metaseq,
  lvl = "Species", norm = F,
  aggfun = colSums, out = "MRexperiment", alternate = T
)

#' **Count number of unique viral species:**
n_species <- length(unique(featureData(BMV_metaseq_species)$Species))

#' **Assign mosquito species for each sample:**
mosquito_species <- pData(BMV_metaseq_species)$SKA_Subspecies

#' **Assign colors to location:**
location <- pData(BMV_metaseq_species)$Municipality
unique_locations <- length(unique(pData(BMV_metaseq_species)$Municipality))
locColors <- viridis::plasma(unique_locations, 1, begin = 0, end = 0.9, direction = 1)
names(locColors) <- levels(as.factor(pData(BMV_metaseq_species)$Municipality))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(locColors, main = "Location colors")

#' **Set heatmap colors:**
heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(200)
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(heatmapCols, main = "Heatmap colors for abundance")

#' **Calculate average BLASTx values:**
# blastx <- read.table("data/BEmosquitoes.tsv", header=T, row.names=1, dec=".", sep="\t")
blastx <- read.table("data/BMV_2023.pid.tsv", header = T, row.names = 1, dec = ".", sep = "\t")
blastx <- dplyr::select(blastx, 2)
blastx.UF <- phyloseq::otu_table(as.matrix(blastx), taxa_are_rows = T)
blastx.ps <- phyloseq(blastx.UF, tax_table(prune_taxa(taxa_sums(BMV.V2) >= 1, BMV.V2)))
blastx_metaseq <- phyloseq_to_metagenomeSeq(blastx.ps)
blastx_metaseq

# Aggregate by species
blastx_mean <- aggregateByTaxonomy(blastx_metaseq, lvl = "Species", norm = F, aggfun = mean, out = "MRexperiment", alternate = T)
blastx_mean

#' **Store average BLASTx identities in dataframe:**
rowanno <- as.data.frame(returnAppropriateObj(blastx_mean, log = F, norm = F))
colnames(rowanno)[1] <- "blastx"
#+ echo=FALSE
datatable(rowanno)
table(rowanno$blastx > 95)

#' **Create color function for BLASTx values:**
col_fun <- colorRamp2(c(0, 100), c("white", "deepskyblue4"))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(col_fun(0:100), main = "Blastx % identity colors")

#' **Get Baltimore classification/genome composition for each viral species from ICTV:**
taxV <- as.data.frame(tax_table(BMV_final))
baltimore <- read.table("data/ICTV_classification_family.tsv", header = T, sep = "\t")

balt_gc <- baltimore %>%
  select(Family, Genome.Composition) %>%
  distinct() %>%
  add_row(Family = "Picorna-like", Genome.Composition = "ssRNA(+)") %>%
  add_row(Family = "Orthomyxo-like", Genome.Composition = "ssRNA(-)") %>%
  add_row(Family = "Negeviridae", Genome.Composition = "ssRNA(+)") %>%
  add_row(Family = "Luteoviridae", Genome.Composition = "ssRNA(+)")

tax_gc <- left_join(taxV, balt_gc, by = "Family") %>%
  filter(!Genome.Composition %in% c("ssRNA(+/-)", "ssDNA(+/-)", "ssDNA(+)", "ssRNA")) %>% # remove double genome compositions
  mutate(Genome.Composition = case_when(
    Species == "Sclerotinia sclerotiorum dsRNA virus 3" ~ "dsRNA",
    Family == "Reo-like" ~ "dsRNA",
    Species == "unclassified" ~ "unknown",
    Family == "Totiviridae" ~ "dsRNA",
    Species == "Hubei orthoptera virus 4" ~ "ssRNA(+)",
    Species == "Hubei odonate virus 15" ~ "dsRNA",
    Species == "Wuhan insect virus 26" ~ "dsRNA",
    Species == "Wuhan insect virus 13" ~ "ssRNA(+)",
    Species == "Linepithema humile qinvirus-like virus 1" ~ "ssRNA(-)",
    Species == "Hubei virga-like virus 23" ~ "ssRNA(+)",
    Species == "Culex mononega like virus 2" ~ "ssRNA(-)",
    TRUE ~ `Genome.Composition`
  ))

gcanno <- tax_gc %>% select(Genome.Composition, Family)
rownames(gcanno) <- tax_gc$Species
colnames(gcanno) <- c("GC", "Family")
gcanno <- gcanno[order(rownames(gcanno)), , drop = F]
gcanno$Family <- factor(gcanno$Family, levels = FamLevel)
#+ echo=FALSE
datatable(gcanno)

#' ## Alpha diversity
df <- as_tibble(sample_data(BMV_final))
df$LibrarySize <- sample_sums(BMV_final)
df <- df[order(df$LibrarySize), ]

threshold <- quantile(df$LibrarySize, 0.05)

min_depth <- df %>%
  dplyr::filter(LibrarySize >= 100) %>%
  select(LibrarySize) %>%
  min()

ggplot(data = df, aes(x = LibrarySize, y = "samples", color = SKA_Subspecies)) +
  geom_jitter(height = 0.02) +
  ggtitle("Library sizes") +
  geom_vline(xintercept = min_depth, linetype = "dashed", color = "grey") +
  labs(y = "") +
  scale_x_log10() +
  scale_color_manual(values = myColors, name = "Mosquito species") +
  theme_bw()

set.seed(1234)
alpha_df_list <- purrr::map(1:1000, ~ alpha_rarefied(ab_table = as_tibble(otu_table(BMV_final)), sequencing_depth = min_depth))
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df %<>%
  mutate(Observed = round(Observed))

labels <- c(
  expression(paste(italic("Aedes japonicus"), " (n=8)")),
  expression(paste(italic("Culex pipiens molestus"), " (n=23)")),
  expression(paste(italic("Culex pipiens pipiens"), " (n=48)")),
  expression(paste(italic("Culex torrentium"), " (n=3)"))
)

alpha <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  left_join(meta %>% select(Sample, SKA_Subspecies)) %>%
  pivot_longer(c(-SKA_Subspecies, -Sample), names_to = "Metric", values_to = "Diversity") %>%
  ggplot(aes(x = SKA_Subspecies, y = Diversity, color = SKA_Subspecies)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  theme_bw() +
  # scale_color_viridis_d(begin=0.4, end =0.9, name="", labels=labels)+
  scale_color_viridis_d(begin = 0, end = .9, name = "", labels = labels) +
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text.align = 0
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  stat_pwc(
    method = "wilcox.test", aes(
      x = SKA_Subspecies, y = Diversity,
      group = SKA_Subspecies
    ), p.adjust.method = "BH",
    label = "p.adj", hide.ns = "p.adj", show.legend = F, tip.length = 0.01
  )
alpha

# Without rarefaction
# alpha <- plot_richness(BMV_final, measures=c("Observed", "Shannon", "Simpson"), x="SKA_Subspecies", color = "SKA_Subspecies")+
#  geom_boxplot()+
#  geom_jitter(width=0)+
#  theme_bw()+
#  scale_color_viridis_d(begin=0.4, end =0.9, name="", labels=c(expression(paste(italic("Aedes japonicus"), " (n=8)")),
#                                                               expression(paste(italic("Culex pipiens molestus"), " (n=23)")),
#                                                               expression(paste(italic("Culex pipiens pipiens"), " (n=48)")),
#                                                               expression(paste(italic("Culex torrentium"), " (n=3)"))))+
#  theme(strip.text.x = element_text(size = 10),
#        axis.title.x = element_blank(),
#        axis.text.x = element_blank(),
#        axis.ticks.x = element_blank(),
#        legend.text.align = 0)+
#  scale_y_continuous(limits = c(0, NA))+
#  stat_pwc(method = "wilcox.test", aes_string(x="SKA_Subspecies", y="value",
#                                              group="SKA_Subspecies"), p.adjust.method="BH",
#           label = "p.adj", hide.ns="p.adj", show.legend = F, tip.length = 0.01)
# alpha

ggsave("figures/alpha-diversity-species.pdf", alpha, dpi = 300)

#' ## Ordination
#' ### PCoA

# Without rarefaction
# v.ord <- ordinate(BMV_final, method = "PCoA")
# pcoa <- plot_ordination(BMV_final, v.ord, type="samples", color="SKA_Subspecies", shape = "Municipality")+
#  scale_shape_manual(values = c(15,16,17,3:8), guide =guide_legend(label.theme = element_text(size=10)), name="Location")+
#  scale_color_viridis_d(begin=0, end =1, name="Mosquito species")+
#  stat_ellipse(type = "norm", linetype = 2, aes_string(group="SKA_Subspecies"), show.legend = F)+
#  theme_bw()+
#  theme(panel.grid.minor = element_blank())+
#  guides(col = guide_legend(override.aes = list(shape = 15, size = 5),
#                            label.theme = element_text(size=10, face="italic")))
#
# pcoa

vegan_avgdist <- avgdist(as.data.frame(t(otu_table(BMV_final))),
  sample = min_depth, iterations = 1000
)

v.ord <- ape::pcoa(vegan_avgdist)

pcoa <- plot_ordination(BMV_final, v.ord, type = "samples", color = "SKA_Subspecies", shape = "Municipality", axes = c(1, 2)) +
  scale_shape_manual(values = c(15, 16, 17, 3:8), guide = guide_legend(label.theme = element_text(size = 10)), name = "Location") +
  scale_color_viridis_d(begin = 0, end = .9, name = "Mosquito species") +
  stat_ellipse(type = "norm", linetype = 2, aes_string(group = "SKA_Subspecies"), show.legend = F) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  guides(col = guide_legend(
    override.aes = list(shape = 15, size = 5),
    label.theme = element_text(size = 10, face = "italic")
  ))
pcoa
#' Permanova test
samples_vegan <- row.names(as.matrix(vegan_avgdist))

metadata_vegan <- meta %>%
  filter(Sample %in% samples_vegan)

perm <- adonis2(vegan_avgdist ~ SKA_Subspecies * Municipality * Trap,
  data = metadata_vegan
)
perm
pval <- perm$`Pr(>F)`[1]
Rsq <- perm$R2[1]

pcoa <- pcoa +
  ggplot2::annotate(geom = "text", x = -.9, y = c(.7, .63, .57), size = c(4, 3, 3), hjust = 0, label = c(
    "Adonis test",
    as.expression(bquote(italic(R^2) ~ "=" ~ .(round(Rsq, 2)))),
    as.expression(bquote(italic(p) ~ "=" ~ .(pval)))
  ))
pcoa
ggsave("figures/PCoA-eukaryotic-virome.pdf", dpi = 300)

#' Beta dispersion
bd <- betadisper(as.dist(vegan_avgdist), metadata_vegan$SKA_Subspecies)
permutest(bd)

#' ### NMDS full dataset

# v.ord <- ordinate(BMV_final, method = "NMDS", k=2)

v.ord <- metaMDS(vegan_avgdist, k = 2)

nmds <- plot_ordination(BMV_final, v.ord, type = "samples", color = "SKA_Subspecies", shape = "Municipality") +
  scale_shape_manual(values = c(15, 16, 17, 3:8), guide = guide_legend(label.theme = element_text(size = 10)), name = "Location") +
  scale_color_viridis_d(begin = 0, end = 1, name = "Mosquito species") +
  # stat_ellipse(type = "norm", linetype = 2, aes_string(group="SKA_Subspecies"), show.legend = F)+
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  guides(col = guide_legend(
    override.aes = list(shape = 15, size = 5),
    label.theme = element_text(size = 10, face = "italic")
  ))
nmds
ggsave("figures/NMDS-eukaryotic-virome.pdf", nmds, dpi = 300)

#' ### NMDS singletons removed
BMV_no_singletons <- BMV_final
for (i in c(
  "MEMO014", "MEMO078", "MEMO145", "MEMO146", "MEMO149", "MEMO084", "MEMO018", "MEMO006",
  "NEMO27", "NEMO30", "NEMO07", "NEMO46", "NEMO51", "MEMO127", "MEMO008"
)) {
  BMV_no_singletons <- subset_samples(BMV_no_singletons, Sample != i)
}
BMV_no_singletons

# v.ord <- ordinate(BMV_no_singletons, method = "NMDS", k=2)
vegan_avgdist_no_single <- avgdist(as.data.frame(t(otu_table(BMV_no_singletons))),
  sample = min_depth, iterations = 1000
)


v.ord_no_single <- metaMDS(vegan_avgdist_no_single, k = 2)

nmds_ns <- plot_ordination(BMV_no_singletons, v.ord_no_single, type = "samples", color = "SKA_Subspecies", shape = "Municipality") +
  scale_shape_manual(values = c(15, 16, 17, 3:8), guide = guide_legend(label.theme = element_text(size = 10)), name = "Location") +
  scale_color_viridis_d(begin = 0, end = .9, name = "Mosquito species") +
  # stat_ellipse(type = "norm", linetype = 2, aes_string(group="SKA_Subspecies"), show.legend = F)+
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  guides(col = guide_legend(
    override.aes = list(shape = 15, size = 5),
    label.theme = element_text(size = 10, face = "italic")
  ))
nmds_ns
ggsave("figures/NMDS-eukaryotic-virome-no-singletons.pdf", nmds_ns, dpi = 300)

v.ord <- ordinate(BMV_no_singletons, method = "PCoA")
pcoa_ns <- plot_ordination(BMV_no_singletons, v.ord, type = "samples", color = "SKA_Subspecies", shape = "Municipality") +
  scale_shape_manual(values = c(15, 16, 17, 3:8), guide = guide_legend(label.theme = element_text(size = 10)), name = "Location") +
  scale_color_viridis_d(begin = 0, end = 1, name = "Mosquito species") +
  stat_ellipse(type = "norm", linetype = 2, aes_string(group = "SKA_Subspecies"), show.legend = F) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  guides(col = guide_legend(
    override.aes = list(shape = 15, size = 5),
    label.theme = element_text(size = 10, face = "italic")
  ))
pcoa_ns

#' Permanova test

samples_vegan_ns <- row.names(as.matrix(vegan_avgdist_no_single))

metadata_vegan_ns <- meta %>%
  filter(Sample %in% samples_vegan_ns)

perm_ns <- adonis2(vegan_avgdist_no_single ~ SKA_Subspecies * Municipality,
  data = metadata_vegan_ns
)
perm_ns
pval_ns <- perm_ns$`Pr(>F)`[1]
Rsq_ns <- perm_ns$R2[1]

ord.legend <- cowplot::get_legend(pcoa)
alpha.legend <- cowplot::get_legend(alpha)
plots <- cowplot::align_plots(alpha + theme(legend.position = "none"),
  pcoa + theme(legend.position = "none"),
  align = "v", axis = "l"
)
bottom <- cowplot::plot_grid(plots[[2]],
  nmds + theme(legend.position = "none"),
  labels = c("B", "C")
)
cowplot::plot_grid(plots[[1]], alpha.legend, bottom, ord.legend, labels = c("A", ""), ncol = 2, rel_widths = c(1, .3))

ggsave("figures/combined-alpha-ordination.pdf", dpi = 300, width = 11.2, height = 8)

#' ## Relative abundance
#' ### Relative abundance per location
phylobar_location <- merge_samples(BMV_family, "Municipality")
phylobar_location_re.abund <- phylobar_location %>%
  transform_sample_counts(function(x) {
    round(100 * (x / sum(x)), 5)
  })

lc <- plot_relabund(phylobar_location_re.abund, fill = "Family") +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack") +
  ylab("Relative abundance (%)") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, vjust = 0.5),
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.ticks.y = element_blank(),
    legend.text = element_text(face = "italic")
  ) +
  scale_fill_manual(values = myFamCol, name = "Viral family") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = rev(gsub("_", "-", loc_order))) +
  coord_flip()

#+ echo=FALSE, fig.width=10
lc

#+ echo=TRUE
#' **Combine viral family barplot with mosquito species barplots from each location:**
lc2 <- ggarrange(lc + theme(legend.position = "none"), ggarrange(plotlist = plist, nrow = 10, heights = c(0.2, rep(1, 8), 0.6)), ncol = 2, widths = c(1, 0.1))
ggsave("figures/viralfamily_barplot_loc.pdf", lc2, height = 6.92, width = 11.2, dpi = 300)

#+ echo=FALSE
lc2

#' ### Relative abundance per mosquito species
phylobar_species <- merge_samples(BMV_family, "SKA_Subspecies")
phylobar_species_re.abund <- phylobar_species %>%
  transform_sample_counts(function(x) {
    round(100 * (x / sum(x)), 5)
  })

sp <- plot_relabund(phylobar_species_re.abund, fill = "Family") +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack") +
  ylab("Relative abundance (%)") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5, hjust = 0.5, face = "italic"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, vjust = 0.5),
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.text = element_text(face = "italic"),
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = myFamCol, name = "Viral family") +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_x_discrete(guide = guide_axis(n.dodge=2))
  scale_x_discrete(labels = c("Ae. japonicus", "Cx. p. molestus", "Cx. p. pipiens", "Cx. torrentium"))
ggsave("figures/viralfamily_barplot_species.pdf", sp, height = 6.92, width = 11.2, dpi = 300)

#+ echo=FALSE
sp

#' **Create combined legend for viral families and mosquito species:**
#+ echo=TRUE
sp2 <- plot_relabund(phylobar_species_re.abund, fill = "Family", color = "Sample") +
  theme_bw() +
  theme(
    legend.text = element_text(face = "italic"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.box.just = "bottom"
  ) +
  scale_fill_manual(values = myFamCol, name = "a") +
  scale_color_manual(values = myColors, name = "b") +
  guides(
    fill = guide_legend(nrow = 4),
    col = guide_legend(override.aes = list(fill = myColors), direction = "vertical")
  )
sp_leg <- get_legend(sp2)
sp2

#' ### Combine relative abundances
#+ fig.width=15, fig.height=9
ggarrange(sp, lc2,
  common.legend = T, legend = "top", labels = "AUTO",
  font.label = list(size = 16), widths = c(0.5, 1), legend.grob = sp_leg
)
ggsave("figures/viralfamily_barplot.pdf", width = 15, height = 9, dpi = 300)

#' ## Eukaryotic virome heatmap
# p.fungi<-log10(pData(BMV_metaseq_species)$fungi_reads+1)
# p.fungi<-pData(BMV_metaseq_species)$fungi_reads
# p.fungi<-as.matrix(fungi2[fungi2$Sample %in% pData(BMV_metaseq_species)$Sample,2:5])
p.fungi <- log10(as.matrix(fungi2[fungi2$Sample %in% pData(BMV_metaseq_species)$Sample, 2:5]))
p.fungi[p.fungi == -Inf] <- NA
colnames(p.fungi) <- as.data.frame(tax_table(BMV_fungi))$Phylum

left_ra <- rowAnnotation(
  "Family" = gcanno$Family,
  "Blastx" = rowanno$blastx,
  col = list(
    "Family" = myFamCol,
    "Blastx" = col_fun
  ),
  show_annotation_name = T,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_name_rot = 270,
  annotation_legend_param = list(
    "Blastx" = list(
      title = "Blastx % Identity",
      direction = "horizontal",
      at = c(0, 25, 50, 75, 100)
    ),
    "Family" = list(title = "Viral family", labels_gp = gpar(fontface = "italic"))
  )
)

# Draw heatmap
column_ha <- HeatmapAnnotation(
  Location = location,
  "Mosquito species" = mosquito_species,
  Fungi = anno_points(p.fungi,
    ylim = c(0, 6),
    axis_param = list(side = "right"),
    gp = gpar(fill = c(2:4, 7), col = c(2:4, 7)),
    pch = 1:length(colnames(p.fungi))
  ),
  show_annotation_name = T,
  annotation_label = gt_render(c("Location", "Species", "log<sub>10</sub>(Fungi)")),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_name_rot = 0,
  annotation_legend_param = list(
    "Mosquito species" = list(
      title = "Mosquito species",
      labels_gp = gpar(fontface = "italic")
    ),
    "Location" = list(title = "Location")
  ),
  col = list("Mosquito species" = c(myColors), Location = c(locColors))
)

# n depends on the amount of taxa in 'BMV_metaseq_species'
hm <- plot_abundance(BMV_metaseq_species,
  n = n_species, log = T, norm = F, colclust = "bray",
  col = heatmapCols, name = "Log2 Read Counts",
  top_annotation = column_ha,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = FALSE,
  heatmap_legend_param = list(direction = "horizontal"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  row_split = factor(gcanno$GC,
    levels = c("dsDNA", "ssDNA", "dsRNA", "ssRNA(+)", "ssRNA(-)", "unknown")
  ),
  row_order = tax_gc[with(tax_gc, order(tax_gc$Genome.Composition, Family)), "Species"],
  left_annotation = left_ra,
  row_title_gp = gpar(fontsize = 10),
  row_title_rot = 0,
  border = T
)

#' Legends
lgd_list <- list(Legend(
  labels = colnames(p.fungi), title = "Fungi", type = "points", pch = 1:length(colnames(p.fungi)),
  legend_gp = gpar(col = c(2:4, 7)), background = "transparent"
))
#+ echo=TRUE, fig.width=15, fig.height=9
draw(hm,
  heatmap_legend_side = "left", annotation_legend_side = "left",
  merge_legend = T, legend_grouping = "original",
  annotation_legend_list = lgd_list
)
#+ include=FALSE
pdf("figures/BEmosquitoes_hm_baltimore.pdf", width = 16.6, height = 10)
draw(hm,
  heatmap_legend_side = "left", annotation_legend_side = "left",
  merge_legend = T, legend_grouping = "original",
  annotation_legend_list = lgd_list
)
dev.off()

#' ***
#' # RT-qPCR analysis
#' ## Prepare the data
library(ggh4x)
qPCR <- read_csv("data/BMV-RT-qPCR-analysis.csv", col_names = T)

#' **Recalculate quantity based on mosquito dilutions:**
qPCR <- qPCR %>%
  mutate(Quantity = case_when(
    str_detect(Sample, "MEMO") ~ qPCR$"Quantity Mean" * 40,
    str_detect(Sample, "NEMO") ~ qPCR$"Quantity Mean" * 100
  ))
#+ echo=FALSE
datatable(qPCR)

#' **Divide quantity of *Xanthi chryso-like virus* by 2 because it is a dsRNA virus:**
qPCR <- qPCR %>%
  mutate(Quantity = case_when(
    str_detect(Target, "XCV") ~ qPCR$Quantity / 2,
    TRUE ~ Quantity
  ))
#+ echo=FALSE
datatable(qPCR)

#' **Replace NA to 0:**
qPCR_noNA <- qPCR %>% mutate_all(~ replace(., is.na(.), 0))

#' **Replace abbreviations for full names:**
qPCR_noNA <- qPCR_noNA %>% mutate(Target, Target = case_when(
  Target == "CPV" ~ "Culex orthophasmavirus 2 (CPV)",
  Target == "XCV" ~ "Xanthi chryso-like virus (XCV)",
  Target == "WMV6" ~ "Wuhan Mosquito Virus 6 (WMV6)",
  Target == "WMV4" ~ "Wuhan Mosquito Virus 4 (WMV4)",
  Target == "Daeseongdong" ~ "Daeseongdong virus 2 (DV2)",
  Target == "HMV4" ~ "Hubei mosquito virus 4 (HMV4)",
  TRUE ~ Target
))

#' **Merge metadata with qPCR data:**
#' Abrreviate locations to fit the names on the plots later on.
metadata <- meta %>% mutate(Municipality = case_when(
  Municipality == "Frameries" ~ "Fr",
  Municipality == "Dilsen-Stokkem" ~ "DS",
  Municipality == "Kallo" ~ "K",
  TRUE ~ Municipality
))
metaqPCR <- merge.data.frame(qPCR_noNA, metadata, by = "Sample")
#+ echo=FALSE
datatable(metaqPCR)

#' **Count number of tested samples per mosquito species:**
metaqPCR %>%
  select(Sample, SKA_Subspecies) %>%
  distinct() %>%
  group_by(SKA_Subspecies) %>%
  summarise(n = n())

#' **Create infection rate table**

inf_rate_loc <- metaqPCR %>%
  select(Sample, Target, Municipality, SKA_Subspecies, Quantity) %>%
  group_by(Target, Municipality, SKA_Subspecies) %>%
  summarise(
    positive = sum(Quantity > 0, na.rm = TRUE), total = n(),
    perc = positive / total * 100
  ) %>%
  mutate(Municipality = case_when(
    Municipality == "Fr" ~ "Frameries",
    Municipality == "DS" ~ "Dilsen-Stokkem",
    Municipality == "K" ~ "Kallo",
    TRUE ~ Municipality
  )) %>%
  filter(perc > 0) %>%
  rename(
    Virus = Target, Species = SKA_Subspecies,
    "Positive mosquitoes" = positive,
    "Total mosquitoes" = total,
    "Infection rate" = perc
  )
xlsx::write.xlsx(
  x = data.frame(inf_rate_loc), file = "tables/supplementary_table1.xlsx",
  row.names = F, showNA = F
)


inf_rate <- metaqPCR %>%
  select(Sample, Target, SKA_Subspecies, Quantity) %>%
  group_by(Target, SKA_Subspecies) %>%
  summarise(
    positive = sum(Quantity > 0, na.rm = TRUE), total = n(),
    perc = positive / total * 100
  ) %>%
  rename(
    Virus = Target, Species = SKA_Subspecies,
    "Positive mosquitoes" = positive,
    "Total mosquitoes" = total,
    "Infection rate" = perc
  )

xlsx::write.xlsx(
  x = data.frame(inf_rate), file = "tables/supplementary_table2.xlsx",
  row.names = F, showNA = F
)


#' ## Plot data
legendlabels <- c(
  expression(paste(italic("Aedes japonicus"), " (n=8)")),
  expression(paste(italic("Culex pipiens molestus"), " (n=47)")),
  expression(paste(italic("Culex pipiens pipiens"), " (n=127)")),
  expression(paste(italic("Culex torrentium"), " (n=16)"))
)

qPCRviolin <- ggplot(metaqPCR, aes(x = SKA_Subspecies, y = Quantity)) +
  geom_violin(aes(fill = SKA_Subspecies, color = SKA_Subspecies), position = "identity", scale = "width") +
  geom_jitter(aes(color = SKA_Subspecies), size = 1, width = 0.1) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.text.align = 0
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0), add = c(0, 0.1)),
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^seq(2, 8, 2)),
    labels = scientific_10(c(0, 10^seq(2, 8, 2)))
  ) +
  xlab(label = "") +
  ylab(label = "Viral copies per mosquito") +
  scale_fill_viridis_d(
    begin = 0, end = 0.9, name = "", alpha = 0.5,
    labels = legendlabels
  ) +
  scale_color_viridis_d(
    begin = 0, end = 0.9, name = "",
    labels = legendlabels
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~Target, nrow = 2, ncol = 3)
ggsave("figures/qPCR-violin.pdf", dpi = 300, width = 11)
ggsave("~/Documenten/Doctoraat/Doctoral school/Research seminar/violin.pdf", dpi = 300, width = 11)
#+ echo=FALSE
qPCRviolin


metaqPCR %>%
  select(Sample, SKA_Subspecies) %>%
  distinct() %>%
  left_join(data.frame(SKA_Subspecies = names(myColors), color = myColors))

color_name <- metaqPCR %>%
  select(Sample, SKA_Subspecies) %>%
  distinct() %>%
  left_join(data.frame(SKA_Subspecies = names(myColors), color = myColors)) %>%
  pull(color, Sample)

#' **Plot qPCR overview per sample:**
qPCRoverview <- ggplot(
  metaqPCR,
  aes(x = Sample, y = Quantity)
) +
  facet_nested(Target ~ Municipality,
    space = "free",
    scales = "free_x",
    labeller = labeller(Target = label_wrap_gen(10), Municipality = label_wrap_gen(5))
  ) +
  geom_col(aes(fill = SKA_Subspecies), position = "identity", alpha = 1) +
  geom_tile(aes(y = -.6, fill = SKA_Subspecies), height = .5) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0), add = c(0, 0.1)),
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10^seq(2, 8, 2)),
    labels = scientific_10(c(0, 10^seq(2, 8, 2)))
  ) +
  xlab(label = "Sample") +
  ylab(label = "Viral copies per mosquito") +
  scale_fill_viridis_d(begin = 0, end = 0.9, name = "") +
  guides(fill = guide_legend(
    override.aes = list(alpha = 1),
    label.theme = element_text(size = 10, angle = 0, face = "italic")
  )) +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank()
  )

ggsave("figures/qPCR_overview.pdf", width = 18.7, height = 10.3, units = "in", dpi = 300)
#+ echo=FALSE, fig.width=18.7, fig.height=10.3, fig.cap="DS=Dilsen-Stokkem, Fr=Frameries, K=Kallo"
qPCRoverview

qPCR_vcount <- metaqPCR %>%
  mutate(Present = if_else(Quantity > 0, 1, 0, 0)) %>%
  select(Sample, Target, Present, SKA_Subspecies) %>%
  group_by(Sample, SKA_Subspecies) %>%
  summarise(qPCR_count = sum(Present))

vcount_df <- left_join(qPCR_vcount, rownames_to_column(vcount), by = c("Sample" = "rowname")) %>%
  replace(is.na(.), 0)

test <- vcount_df %>%
  group_by(SKA_Subspecies) %>%
  dplyr::count(qPCR_count, viral_count)

test$viral_count <- as.numeric(test$viral_count)
test$qPCR_count <- as.numeric(test$qPCR_count)

test1 <- test %>%
  pivot_wider(names_from = SKA_Subspecies, values_from = n) %>%
  replace(is.na(.), 0) %>%
  rowwise() %>%
  mutate(totalcount = sum(c_across(3:6)))

colnames(test1)[4] <- "Culex pipiens molestus"
colnames(test1)[5] <- "Culex pipiens pipiens"

pieCols <- myColors
names(pieCols) <- colnames(test1)[3:6]

# ggplot()+
#  scatterpie::geom_scatterpie(data=test1, aes(x=viral_count, y=qPCR_count, r=sqrt(totalcount)/15),
#                              cols=c("Culex pipiens molestus", "Culex pipiens pipiens", "Culex torrentium", "Aedes japonicus"))+
#  scatterpie::geom_scatterpie_legend(radius = sqrt(test1$totalcount)/15, x = 6.5, y=0, n=3, labeller = function(x) (15*x)^2)+
#  scale_fill_manual(values = pieCols)+
#  coord_equal()+
#  theme_minimal()

test2 <- test1 %>%
  pivot_longer(names_to = "Mosquito", values_to = "Amount", cols = 3:6)

points <- vcount_df %>%
  group_by(SKA_Subspecies)

p <- ggplot() +
  geom_rect(data = test2[test2$Mosquito == "Aedes japonicus", ], aes(xmin = viral_count - 0.5, xmax = viral_count, ymin = qPCR_count, ymax = qPCR_count + 0.5), fill = "#2A788EFF", color = NA) +
  geom_rect(data = test2[test2$Mosquito == "Culex pipiens molestus", ], aes(xmin = viral_count, xmax = viral_count + 0.5, ymin = qPCR_count, ymax = qPCR_count + 0.5), fill = "#1FA188FF", color = NA) +
  geom_rect(data = test2[test2$Mosquito == "Culex pipiens pipiens", ], aes(xmin = viral_count - 0.5, xmax = viral_count, ymin = qPCR_count - 0.5, ymax = qPCR_count), fill = "#54C568FF", color = NA) +
  geom_rect(data = test2[test2$Mosquito == "Culex torrentium", ], aes(xmin = viral_count, xmax = viral_count + 0.5, ymin = qPCR_count - 0.5, ymax = qPCR_count), fill = "#BBDF27FF", color = NA) +
  geom_text(data = test2[test2$Mosquito == "Aedes japonicus", ], aes(x = viral_count - 0.25, y = qPCR_count + 0.25, label = Amount), color = "white") +
  geom_text(data = test2[test2$Mosquito == "Culex pipiens molestus", ], aes(x = viral_count + 0.25, y = qPCR_count + 0.25, label = Amount), color = "white") +
  geom_text(data = test2[test2$Mosquito == "Culex pipiens pipiens", ], aes(x = viral_count - 0.25, y = qPCR_count - 0.25, label = Amount), color = "white") +
  geom_text(data = test2[test2$Mosquito == "Culex torrentium", ], aes(x = viral_count + 0.25, y = qPCR_count - 0.25, label = Amount), color = "white") +
  geom_label(data = test2[test2$Mosquito == "Culex torrentium", ], aes(x = viral_count, y = qPCR_count, label = totalcount), color = "black", fontface = "bold", size = 6) +
  geom_tile(data = test2, aes(x = viral_count, y = qPCR_count), fill = "transparent", col = "white", lwd = 1) +
  scale_x_continuous("# viral species NGS", breaks = c(0:8), expand = c(0, 0)) +
  scale_y_continuous("# viral species qPCR", expand = c(0, 0)) +
  # geom_point(data=points, aes(x=viral_count, y=qPCR_count), color="transparent")+
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    # aspect.ratio = 4/9,
    plot.margin = unit(c(.1, .1, .1, .1), "mm")
  )

ggsave("figures/viral_count.pdf", dpi = 300, width = 7.5, height = 5.25)

p1 <- points %>%
  ggplot(aes(x = viral_count)) +
  geom_histogram(binwidth = 1, fill = "grey", colour = "white") +
  stat_bin(binwidth = 1, geom = "text", aes(label = after_stat(count)), vjust = -.5) +
  # ylim(c(0,140))+
  coord_cartesian(ylim = c(0, 140)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))

p2 <- points %>%
  ggplot(aes(x = qPCR_count)) +
  geom_histogram(binwidth = 1, fill = "grey", colour = "white") +
  stat_bin(binwidth = 1, geom = "text", aes(label = after_stat(count)), angle = -90, vjust = -.1) +
  # ylim(c(0,110))+
  coord_cartesian(ylim = c(0, 110)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_flip() +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))

ggarrange(p1, NULL, NULL,
  NULL, NULL, NULL,
  p, NULL, p2,
  widths = c(.9, -.035, .1), heights = c(.2, -.05, .8), align = "hv"
)

ggsave("figures/viral_count_barplot.pdf", dpi = 300, width = 8, height = 4)

#' ***
#' # Phylogenetics
#+ results="hide"
source("BMV_phylogenetics.R")

#' ***
#' # Phageome & Wolbachia bacteria
#' ## Coverage plot

phageblast <- read.table("data/phage_blast.tsv", header = F, sep = "\t")
phageblast <- phageblast[phageblast$V4 > 2000, ]

rRNA <- data.frame(
  name = c("16S", "23S", "5S"),
  start = c(1136001, 1235938, 1238768),
  end = c(1137446, 1238682, 1238878)
)

prophage <- data.frame(
  name = c("p1", "p2", "p3", "p4", "p5"),
  start = c(250652, 320006, 346157, 444716, 1376656),
  end = c(266846, 337904, 360657, 483289, 1413160)
)

chunk_size <- 2500

MEMO_wolb <- read.table("data/depth_files/MEMO050.bowtie.depth", header = T)
# MEMO_wolb <- read.table("data/depth_files/MEMO043.bowtie.depth", header = T)

averages <- MEMO_wolb %>%
  group_by(chunk = ceiling(row_number() / chunk_size)) %>%
  summarize(across(starts_with("MEMO"), ~ mean(.[!is.na(.)]), .names = "Average")) %>%
  mutate(pos = chunk * chunk_size)

covplot <- ggplot(averages, aes(x = pos, y = Average)) +
  geom_area(color = "black", fill = "black", alpha = 0.5) +
  geom_rect(data = phageblast, inherit.aes = FALSE, aes(
    xmin = V9, xmax = V10, ymin = -max(averages$Average) * 0.066, ymax = -max(averages$Average) * 0.043,
    fill = "#3977AF", color = "#3977AF"
  ), alpha = 1) +
  geom_rect(data = rRNA, inherit.aes = FALSE, aes(
    xmin = start, xmax = end, ymin = -max(averages$Average) * 0.036, ymax = -max(averages$Average) * 0.012,
    fill = "red", color = "red"
  ), alpha = 1) +
  geom_rect(data = prophage, inherit.aes = FALSE, aes(
    xmin = start, xmax = end, ymin = -max(averages$Average) * 0.036, ymax = -max(averages$Average) * 0.012,
    fill = "limegreen", color = "limegreen"
  ), alpha = 1) +
  labs(x = "Nucleotide position", y = "Read depth") +
  # ggtitle("Wolbachia coverage plot of sample MEMO050")+
  scale_fill_identity(
    name = "",
    breaks = c("#3977AF", "limegreen", "red"),
    labels = c("phageWO prophage region (BLASTn)", "prophage regions (Genbank)", "rRNA genes"),
    guide = "legend"
  ) +
  scale_color_identity(
    name = "",
    # breaks = c("#3977AF", "red"),
    # labels = c("phageWO prophage region", "rRNA genes"),
    guide = "none", # "legend"
  ) +
  theme_bw() +
  theme(
    legend.position = c(.2, .97),
    legend.background = element_rect(fill = NA, color = NA),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.key.size = unit(3, "mm"),
    plot.title = element_text(size = 9, face = "bold"),
    # panel.grid.major = element_line(size = .25),
    panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

covplot
# ggsave("figures/wolbachia_coverageplot_MEMO043.pdf", dpi=300, width = 7.5, height=5)

# Define the function
process_file <- function(file_path, sample) {
  MEMO_wolb <- read.table(file_path, header = TRUE)

  averages <- MEMO_wolb %>%
    group_by(chunk = ceiling(row_number() / chunk_size)) %>%
    summarize(across(starts_with(c("MEMO", "NEMO")), ~ mean(.[!is.na(.)]), .names = "Average")) %>%
    mutate(pos = chunk * chunk_size)

  max_y <- max(averages$Average)
  if (max_y > 200) {
    max_y <- 200
  }

  if (sample == "MEMO100" | sample == "MEMO139") {
    max_y <- 15
  }

  covplot <- ggplot(averages, aes(x = pos, y = Average)) +
    geom_area(color = "black", fill = "black", alpha = 0.5) +
    geom_rect(data = phageblast, inherit.aes = FALSE, aes(
      xmin = V9, xmax = V10,
      ymin = -max_y * 0.066, ymax = -max_y * 0.043,
      fill = "#3977AF", color = "#3977AF"
    ), alpha = 1) +
    geom_rect(data = rRNA, inherit.aes = FALSE, aes(
      xmin = start, xmax = end,
      ymin = -max_y * 0.036, ymax = -max_y * 0.012,
      fill = "red", color = "red"
    ), alpha = 1) +
    geom_rect(data = prophage, inherit.aes = FALSE, aes(
      xmin = start, xmax = end,
      ymin = -max_y * 0.036, ymax = -max_y * 0.012,
      fill = "limegreen", color = "limegreen"
    ), alpha = 1) +
    labs(x = "Nucleotide position", y = "Read depth", title = sample) +
    scale_fill_identity(
      name = "",
      breaks = c("#3977AF", "limegreen", "red"),
      labels = c("phageWO prophage region (BLASTn)", "prophage regions (Genbank)", "rRNA genes"),
      guide = "legend"
    ) +
    scale_color_identity(
      name = "",
      guide = "none"
    ) +
    theme_bw() +
    theme(
      legend.position = c(.15, .97),
      legend.background = element_rect(fill = NA, color = NA),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.key.size = unit(3, "mm"),
      plot.title = element_text(size = 9, face = "bold"),
      panel.grid = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))

  covplot <- covplot + coord_cartesian(ylim = c(NA, max_y))

  return(covplot)
}


samples <- c(
  "MEMO005", "MEMO078", "NEMO46", "MEMO141", "NEMO42", "MEMO085", "MEMO073", "NEMO35", "MEMO108",
  "MEMO139", "MEMO057", "MEMO056", "MEMO105", "MEMO129",
  "MEMO100", "MEMO043"
)

# Apply the function to each file path using map
covplots <- purrr::map2(paste0("data/depth_files/", samples, ".bowtie.depth"), samples, process_file)

p <- ggarrange(plotlist = covplots, ncol = 4, nrow = 4, common.legend = T)
p

ggsave(file = "figures/wolbachia_coverageplots.pdf", plot = p, dpi = 300, height = 10, width = 15)

# covplot <- ggplot(MEMO_wolb, aes_string(x=pos, y=names(MEMO_wolb[3])))+
#  geom_line()+
#  geom_rect(data = phageblast, inherit.aes=FALSE, aes(xmin=V9, xmax=V10, ymin=-max(MEMO_wolb[3])*0.066, ymax=-max(MEMO_wolb[3])*0.043,
#                                                      fill="#3977AF", color="#3977AF"), alpha=1)+
#  geom_rect(data = rRNA, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=-max(MEMO_wolb[3])*0.036, ymax=-max(MEMO_wolb[3])*0.012,
#                                                      fill="red", color="red"), alpha=1)+
#  geom_rect(data = prophage, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=-max(MEMO_wolb[3])*0.036, ymax=-max(MEMO_wolb[3])*0.012,
#                                                      fill="limegreen", color="limegreen"), alpha=1)+
#  labs(x="Nucleotide position", y="Read depth")+
#  #ggtitle("Wolbachia coverage plot of sample MEMO050")+
#  scale_fill_identity(name = "",
#                      breaks = c("#3977AF", "limegreen", "red"),
#                      labels = c("phageWO prophage region (BLASTn)", "prophage regions (Genbank)", "rRNA genes"),
#                      guide = "legend")+
#  scale_color_identity(name = "",
#                      #breaks = c("#3977AF", "red"),
#                      #labels = c("phageWO prophage region", "rRNA genes"),
#                      guide = "none",#"legend"
#                      )+
#  theme_bw()+
#  theme(legend.position = c(.3,.93), #"top"
#        legend.background=element_rect(fill=NA, color=NA),
#        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit='cm'),
#        legend.key.size = unit(3, 'mm'),
#        plot.title = element_text(size=9, face="bold"),
#        #panel.grid.major = element_line(size = .25),
#        panel.grid=element_blank())+
#  guides(fill = guide_legend(override.aes = list(alpha = 1)))
#
# covplot


#' ## Swarmplot
library(ggbeeswarm)
# phage <- read.table("data/wolbachia_phage_ratio.tsv", sep="\t", header = T)
# phage <- read.table("data/wolbachia_phage_bowtieratio_no_rRNA.tsv", sep="\t", header = T)
phage <- read.table("data/wolbachia_prophage_gb_bowtieratio_no_rRNA.tsv", sep = "\t", header = T)
phage <- merge(phage, meta, by = "Sample")

phage <- phage %>%
  mutate(ratio = case_when(
    ratio == NA ~ 0,
    nonphage < 1 ~ 0,
    TRUE ~ ratio
  ))

phage %>%
  filter(ratio > 0) %>%
  summarise(median(ratio))

phage <- phage %>%
  mutate(color = case_when(
    ratio == 0 ~ "NA",
    TRUE ~ color
  ))

swarmplot <- ggplot(phage, aes(x = SKA_Subspecies, y = ratio, color = color)) +
  geom_quasirandom(width = .4) +
  scale_color_manual(
    values = c("<3.5" = "#3977AF", ">3.5" = "darkorange", "NA" = "#7F7F7F"),
    name = "",
    labels = c("Prophage", "Real phage particle?", "Average depth < 1")
  ) +
  geom_text(aes(label = ifelse(ratio > 3.5 & nonphage > 1, Sample, "")),
    nudge_x = .5, size = 3,
    show.legend = F
  ) +
  ylab(paste0("<span style='font-size: 11pt'>Depth ratio</span><br><span style='font-size: 8pt'>*(average depth phage/average depth nonphage)*</span>")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "italic", size = 6),
    legend.position = c(.2, .93),
    legend.background = element_blank(),
    legend.key.size = unit(3, "mm"),
    panel.grid = element_blank(),
    axis.title.y = ggtext::element_markdown()
  )
swarmplot


#' ## Wolbachia
#' Wolbachia bwa-mem2 mapped
wMel <- read.table("data/wMel.txt", header = F, dec = ".", sep = " ")
wPip <- read.table("data/wPip.txt", header = F, dec = ".", sep = " ")
phageWO <- read.table("data/phageWO.txt", header = F, dec = ".", sep = " ")
colnames(wMel) <- c("Sample", "numreads", "covbases", "coverage", "meandepth")
colnames(wPip) <- c("Sample", "numreads", "covbases", "coverage", "meandepth")
colnames(phageWO) <- c("Sample", "supergroup", "numreads", "covbases", "coverage", "meandepth")
wPip$supergroup <- "wPip"
wMel$supergroup <- "wMel"

wolb <- rbind(wMel, wPip, phageWO)

wolb %<>%
  plyr::mutate(supergroup = str_replace(string = supergroup, pattern = "NODE_12_length_11674_cov_94.907217_MEMO129", replacement = "phageWO2")) %>%
  plyr::mutate(supergroup = str_replace(string = supergroup, pattern = "NODE_1_length_31159_cov_25.624799_MEMO050", replacement = "phageWO1"))
wolb$numreads[wolb$coverage < 5] <- 0
wolb$meandepth[wolb$coverage < 5] <- 0

left_join(wolb, meta) %>%
  filter(supergroup == "wPip") %>%
  group_by(SKA_Subspecies) %>%
  summarise(n = n(), wol_pos = sum(numreads > 0)) %>%
  mutate(wol_freq = round(wol_pos / n * 100, 3))


wviolin <- pivot_wider(wolb, id_cols = "Sample", names_from = "supergroup", values_from = "numreads") %>%
  ggplot(aes(x = meta[meta$Sample != "NEMO39", ]$SKA_Subspecies, y = wPip)) +
  geom_boxplot(aes(color = meta[meta$Sample != "NEMO39", ]$SKA_Subspecies), width = 0.1, outlier.shape = NA) +
  geom_violin(aes(fill = meta[meta$Sample != "NEMO39", ]$SKA_Subspecies), color = NA, show.legend = F, trim = T) +
  geom_jitter(aes(color = meta[meta$Sample != "NEMO39", ]$SKA_Subspecies),
    size = 1, alpha = 0.6,
    position = position_jitter(w = 0.1, h = 0)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab(label = "") +
  ylab(label = "Wolbachia read count") +
  scale_fill_manual(values = alpha(myColors, 0.1), name = "") +
  scale_color_viridis_d(begin = 0.4, end = 0.9, name = "") +
  guides(fill = guide_legend( # override.aes = list(alpha = 1),
    label.theme = element_text(size = 8, angle = 0, face = "italic")
  )) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 10^seq(2, 8, 2)))
wviolin
ggsave("figures/wolbachia_readcount.pdf", dpi = 300)

ggarrange(
  ggarrange(covplot + theme(legend.key.size = unit(3, "mm")), swarmplot,
    ncol = 1, labels = "AUTO", align = "v"
  ),
  wviolin +
    theme(
      legend.position = c(.15, .92),
      legend.background = element_blank()
    ),
  labels = c("", "C"), widths = c(1, 1.25), align = "v"
)
ggsave("figures/combined_wolbachia.pdf", dpi = 300, width = 12)

p1 <- ggarrange(covplot + theme(legend.key.size = unit(3, "mm")), swarmplot,
  ncol = 1, labels = "AUTO"
)
p2 <- wviolin +
  theme(
    legend.position = c(.15, .93),
    legend.background = element_blank(),
    panel.background = element_blank()
  )
p2

library(patchwork)
patch <- p2 + plot_spacer() + covplot + swarmplot + plot_layout(
  widths = c(1.5, .01, 1),
  design = "
                                                123
                                                124"
)
patch + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/combined_wolbachia.pdf", dpi = 300, width = 12)

#' Second layout
patch2 <- covplot + theme(legend.position = c(0.3, 0.97)) + swarmplot + plot_spacer() + p2 + plot_layout(
  widths = c(1, .01, 1.5),
  design = "
                                                134
                                                234"
)
patch2 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/combined_wolbachia_2.pdf", dpi = 300, width = 12)

#' barplot
wb <- pivot_wider(wolb, id_cols = "Sample", names_from = "supergroup", values_from = "numreads") %>%
  ggplot(aes(x = Sample, y = wPip)) +
  geom_col(aes(fill = meta[meta$Sample != "NEMO39", ]$SKA_Subspecies)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01), add = c(0, 0.1))) +
  xlab(label = "Sample") +
  ylab(label = "Wolbachia read count") +
  labs(shape = "Virus", color = "Virus") +
  scale_fill_viridis_d(begin = 0.4, end = 0.9, name = "") +
  guides(fill = guide_legend(
    override.aes = list(alpha = 1),
    label.theme = element_text(size = 10, angle = 0, face = "italic")
  )) +
  facet_grid(~ meta[meta$Sample != "NEMO39", ]$Municipality, scales = "free", space = "free")
wb


#' ## Phageome
#'
phageotu <- read.table("data/BEmosq_classification-85-taxfile-1000nt.tsv", sep = "\t", header = T)
phage <- phageotu[phageotu$Phylum == "Phage", ]
rownames(phage) <- phage$Contig
phage <- subset(phage, select = -Contig)

newtax <- tax[!rownames(tax) %in% rownames(phage), ]

newtax <- rbind(newtax, phage)

tax_table(BMV) <- as.matrix(newtax) # Maybe better solution

# newtax.UF <- tax_table(as.matrix(newtax))

# BMV <- phyloseq(OTU.UF, newtax.UF, meta.UF)

BMV.V <- subset_taxa(BMV, Kingdom == "Viruses")

phageome <- subset_taxa(BMV.V, Phylum == "Phage" | Order == "Caudovirales" | Family == "Inoviridae")
phageome

otu_table(phageome) <- as.data.frame(otu_table(phageome)) %>%
  mutate(otu = rownames(.)) %>%
  pivot_longer(cols = 1:197, names_to = "Sample", values_to = "abundance") %>%
  mutate(abundance = if_else(otu %in% c(
    "NODE_1_length_31159_cov_25.624799_MEMO050|full",
    "NODE_12_length_11674_cov_94.907217_MEMO129|full"
  ) &
    Sample != "NEMO039", 0, as.double(abundance))) %>%
  pivot_wider(names_from = Sample, values_from = abundance) %>%
  column_to_rownames("otu") %>%
  otu_table(taxa_are_rows = T)

# sample_data(phageome) <- data.frame(wv.meta, row.names = 1)
# sample_data(phageome)

phage_pruned <- prune_samples(sample_sums(phageome) > 0, phageome)
phage_pruned <- prune_taxa(taxa_sums(phage_pruned) > 0, phage_pruned)
phage_pruned

pp <- as.data.frame(otu_table(phage_pruned))

p.merge <- merge_taxa(phage_pruned, 1:2, 2)
p.merge <- merge_taxa(p.merge, 2:24, 2)
p.merge <- merge_taxa(p.merge, 3:4, 2)
tax_table(p.merge)
otu_table(p.merge)
plot_heatmap(p.merge, taxa.label = "Species", low = "#FFF7B8", high = "#800026", na.value = "#FFFFCC") +
  theme_bw() +
  theme(axis.title = element_blank())


phage_tax <- as.data.frame(tax_table(p.merge))
phage_tax <- replace_na(phage_tax, list(
  Phylum = "Uroviricota", Class = "Caudoviricetes", Order = "Caudovirales",
  Family = "unclassified", Genus = "unclassified", Species = "unclassified"
))
tax_table(p.merge) <- as.matrix(phage_tax)

#' ### Wolbachia vs virus
viral_reads <- as.data.frame(colSums(otu_table(BMV.V2)))
viral_reads <- rownames_to_column(viral_reads, "Sample")
colnames(viral_reads) <- c("Sample", "viral_reads")

wol_virus_df <- left_join(wPip, viral_reads)

fit <- lm(log10(viral_reads + 1) ~ coverage, data = wol_virus_df)
summary(fit)

wol_virus_df %>%
  left_join(meta[meta$Sample != "NEMO39", ]) %>%
  # mutate(viral_reads=viral_reads+1) %>%
  ggplot(aes(x = coverage, y = viral_reads)) +
  geom_point(aes(col = SKA_Subspecies)) +
  geom_smooth(method = lm) +
  geom_vline(
    xintercept = 5, linetype = "dashed",
    color = "red"
  ) +
  # scale_x_continuous(trans='log10')+
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 10^seq(2, 8, 2))) +
  # scale_y_continuous(trans='log10')+
  scale_color_manual(values = myColors) +
  ylab("log10(viral reads+1)") +
  xlab("Wolbachia wPip genome coverage") +
  labs(
    color = "Mosquito species",
    title = paste(
      "Adj R2 = ", signif(summary(fit)$adj.r.squared, 5),
      "Intercept =", signif(fit$coef[[1]], 5),
      " Slope =", signif(fit$coef[[2]], 5),
      " P =", signif(summary(fit)$coef[2, 4], 5)
    )
  ) +
  theme_bw() +
  guides(col = guide_legend(
    override.aes = list(shape = 15, size = 5),
    label.theme = element_text(size = 10, face = "italic")
  ))

ggsave("figures/wolbachia_virus_correlation.pdf", dpi = 300)

left_join(wPip, viral_reads) %>%
  ggplot(aes(x = coverage > 5, y = viral_reads)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 10^seq(2, 8, 2))) +
  geom_boxplot() +
  geom_point()

left_join(wPip, rownames_to_column(vcount, "Sample")) %>%
  ggplot(aes(x = coverage, y = viral_count)) +
  # geom_boxplot()+
  geom_point() +
  geom_smooth(method = lm)

#' ### pWCP coverage
pWCP_cov <- read.table("data/pWCP_stats.tsv", header = T)
# names(pWCP_cov) <- c("Sample", "pWCP_cov")
left_join(wPip, pWCP_cov, by = "Sample") %>%
  left_join(meta[meta$Sample != "NEMO39", ]) %>%
  mutate(
    meandepth.x = case_when(
      coverage.x < 5 ~ 0,
      TRUE ~ meandepth.x
    ),
    meandepth.y = case_when(
      coverage.x < 5 ~ 0,
      TRUE ~ meandepth.y
    ),
    completeness = cut(pWCP_cov$coverage, breaks = c(0, 25, 50, 75, 100), include.lowest = T)
  ) %>%
  ggplot(aes(x = meandepth.x, y = meandepth.y)) +
  geom_smooth(method = "lm", se = F, col = "lightgrey") +
  geom_point(aes(col = SKA_Subspecies, shape = completeness)) +
  scale_color_manual(values = myColors) +
  scale_shape_manual(values = c(3, 17, 15, 19)) +
  scale_y_continuous(trans = pseudo_log_trans()) +
  scale_x_continuous(trans = pseudo_log_trans()) +
  labs(
    color = "Mosquito species",
    x = "Wolbachia sequencing depth",
    y = "pWCP plasmid sequencing depth"
  ) +
  theme_bw() +
  guides(col = guide_legend(
    override.aes = list(shape = 15, size = 5),
    label.theme = element_text(size = 10, face = "italic")
  ))

ggsave("figures/pWCP_coverage.pdf", dpi = 300)

left_join(wPip, pWCP_cov, by = "Sample") %>%
  left_join(meta[meta$Sample != "NEMO39", ]) %>%
  mutate(
    meandepth.x = case_when(
      coverage.x < 5 ~ 0,
      TRUE ~ meandepth.x
    ),
    meandepth.y = case_when(
      coverage.x < 5 ~ 0,
      TRUE ~ meandepth.y
    ),
    completeness = cut(pWCP_cov$coverage, breaks = c(0, 25, 50, 75, 100), include.lowest = T)
  ) %>%
  ggplot(aes(x = meandepth.x, y = meandepth.y)) +
  geom_smooth(method = "lm", se = F, col = "lightgrey") +
  geom_point(aes(col = coverage.y, shape = SKA_Subspecies)) +
  scale_color_gradientn(colours = viridis(10), guide = guide_colourbar(direction = "horizontal", title.position = "top")) +
  scale_y_continuous(trans = pseudo_log_trans()) +
  scale_x_continuous(trans = pseudo_log_trans()) +
  labs(
    color = "pWCP completeness (%)",
    shape = "Mosquito species",
    x = "Wolbachia sequencing depth",
    y = "pWCP plasmid sequencing depth"
  ) +
  theme_bw()

# ggsave("figures/pWCP_coverage.pdf", dpi=300)

#+ include=FALSE
spin("BelgianMosquitoVirome.R", knit = F, format = "Rmd")
#' <div class="tocify-extend-page" data-unique="tocify-extend-page" style="height: 0;"></div>
