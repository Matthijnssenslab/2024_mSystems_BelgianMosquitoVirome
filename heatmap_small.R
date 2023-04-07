#' Run main R script first
left_ra=rowAnnotation('Family'=gcanno$Family,
                      "Blastx"= rowanno$blastx,
                      simple_anno_size = unit(2, "mm"),
                      col=list('Family'=myFamCol,
                               "Blastx"=col_fun),
                      show_annotation_name=T,
                      annotation_name_side = "top",
                      annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),
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
                                                  pch=1:length(colnames(p.fungi)),
                                                  size=unit(1, "mm")),
                              simple_anno_size = unit(2, "mm"),
                              height = unit(1, "cm"),
                              show_annotation_name = T,
                              annotation_label = gt_render(c("Location", "Species", "log<sub>10</sub>(Fungi)")),
                              annotation_name_side = "right",
                              annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),
                              annotation_name_rot = 0,
                              annotation_legend_param=list('Mosquito species'= list(title ='Mosquito species', 
                                                                                    labels_gp = gpar(fontface="italic", 
                                                                                                     size=2)),
                                                           'Location' = list(title='Location',
                                                                             labels_gp=gpar(size=2))),
                              col = list('Mosquito species'=c(myColors), Location=c(locColors)))

#n depends on the amount of taxa in 'BMV_metaseq_species'
hm <- plot_abundance(BMV_metaseq_species, n = n_species, log = T, norm = F, colclust = "bray",
                     col = heatmapCols, name = "Log2 Read Counts", 
                     top_annotation=column_ha,
                     row_names_gp = gpar(fontsize = 4.5),
                     show_column_names = FALSE,
                     heatmap_legend_param = list(direction = "horizontal",
                                                 legend_gp = gpar(size = 2)),
                     cluster_rows = FALSE,
                     cluster_row_slices = FALSE,
                     row_split=factor(gcanno$GC, 
                                      levels = c("dsDNA","ssDNA","dsRNA","ssRNA(+)","ssRNA(-)","unknown")),
                     row_order=tax_gc[with(tax_gc, order(tax_gc$Genome.Composition, Family)),"Species"],
                     left_annotation=left_ra,
                     row_title_gp = gpar(fontsize=7),
                     row_title_rot=0,
                     border=T)

pdf("figures/BEmosquitoes_hm_baltimore_test.pdf", width = 6.875, height=4.125)
draw(hm, show_heatmap_legend=F, show_annotation_legend=F)
dev.off()
