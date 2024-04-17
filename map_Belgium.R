library(scatterpie)
library(ggpubr)
library(maps)
library(tidyverse)
library(raster)
library(sf)
library(readxl)
library(ggrepel)
library(rnaturalearth)
library(ggforce)
library(here)
here::i_am("map_Belgium.R")

geom_scatterpie_legend2 <- function (radius, x, y, n = 5, labeller, linetype=1) 
{
  if (length(radius) > n) {
    radius <- unique((seq(min(radius), max(radius), 
                                length.out = n)))
  }
  label <- FALSE
  if (!missing(labeller)) {
    if (!inherits(labeller, "function")) {
      stop("labeller should be a function for converting radius")
    }
    label <- TRUE
  }
  dd <- data.frame(r = radius, start = 0, end = 2 * pi, x = x, 
                   y = y + radius - max(radius), maxr = max(radius))
  if (label) {
    dd$label <- round(labeller(dd$r))
  }
  else {
    dd$label <- dd$r
  }
  list(geom_circle(aes_(x0 = ~x, y0 = ~y, r0 = ~r, r = ~r, 
                         start = ~start, end = ~end), data = dd, inherit.aes = FALSE), 
       geom_segment(aes_(x = ~x, xend = ~x + maxr * 1.5, y = ~y + 
                           r, yend = ~y + r), data = dd, inherit.aes = FALSE, linetype=linetype), 
       geom_text(aes_(x = ~x + maxr * 1.6, y = ~y + r, label = ~label), 
                 data = dd, hjust = "left", inherit.aes = FALSE))
}

round_digit <- function (d) {
  if (d > 1) {
    round(d)
  } else {
    round(d, -as.integer(floor(log10(abs(d)))))
  }
}

environment(geom_scatterpie_legend2) <- asNamespace('scatterpie')

# Get map of Belgium
bel <- getData('GADM', country='BEL', level=2)
# Convert polygon map to ggplot friendly dataframe
belsf <- st_as_sf(bel)
p<-ggplot() + geom_sf(data = belsf)
p

beldf <-map_data(bel)
#world <- map_data('world')
# Read in dataframe
metadata <- read_csv("data/BEmosq_metadata.csv")

# Function for custom theme for maps
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      panel.background = element_rect(fill = "#f5f5f2", color = NA), 
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      ...
    )
}
metadata$lat <- as.numeric(metadata$lat)
metadata$long <- as.numeric(metadata$long)

metadata <- metadata %>% group_by(long, SKA_Subspecies) %>% mutate(count = n())

pie <- metadata %>% dplyr::select(SKA_Subspecies, Municipality, long, lat, count) %>% dplyr::distinct() %>% drop_na()
pie2 <- pie %>% group_by(Municipality, SKA_Subspecies) %>% dplyr::summarise(lat=mean(lat), long=mean(long), count=sum(count))

pie3 <- pie2 %>% dplyr::select(SKA_Subspecies, count) %>% pivot_wider(names_from = SKA_Subspecies, values_from = count)
finalpie <- left_join(pie2 %>% dplyr::select(Municipality, long, lat), pie3, by="Municipality") %>% 
  mutate_if(is.numeric, round, digits=3) %>% dplyr::distinct()

colnames(finalpie)[4]<- "Culex pipiens molestus"
colnames(finalpie)[5]<- "Culex pipiens pipiens"
finalpie

finalpie2 <- finalpie %>% 
  replace(is.na(.), 0) %>% 
  rowwise() %>% 
  dplyr::mutate(totalcount = sum(c_across("Culex pipiens molestus":"Aedes japonicus"))) %>% 
  dplyr::mutate(r = sqrt(totalcount)/50)

finalpie2 <- finalpie2[!duplicated(finalpie2$Municipality),]
colnames(finalpie)
colSums(finalpie2[4:5])
finalpie2
# finalpie2 <- finalpie %>% 
#   replace(is.na(.), 0) %>% 
#   rowwise() %>% 
#   dplyr::mutate(totalcount = sum(c_across(3:7))) %>% 
#   dplyr::mutate(r = totalcount/500)
# colnames(finalpie2)

# Plot of Belgium with mosquito species
p <- ggplot(beldf, aes(x=long, y=lat)) + 
  geom_map(map=beldf, aes(map_id=region), fill='lightgrey', color="white")+
  coord_quickmap()
#p+coord_quickmap()
p+theme_void()
#ggsave("Belgium-provinces.pdf", width=15, height=15, dpi = 300)

# p <- ggplot(beldf, aes(x=long, y=lat)) + 
#   geom_polygon(aes(group=group), fill='lightgrey', color='white')
# p

p2 <- p+ geom_point(data = finalpie2, aes(x=long, y=lat), size=.5)+
  #geom_text(data = finalpie2, aes(label=Municipality))+
  theme_void()#+
p2
#ggsave("Belgium.pdf")
#geom_scatterpie(data = finalpie2, aes(x=long, y=lat, group=Municipality, cols=count))

p3 <- p2+geom_scatterpie(data = finalpie2[!finalpie2$Municipality %in% c("Bertem", "Kallo"),], 
                    aes(x=long, y=lat, group=Municipality, r=r), 
                    cols = c("Culex pipiens molestus", "Culex pipiens pipiens", "Culex torrentium", "Aedes japonicus"),
                  color="transparent", size=0.15, legend_name = 'Subpecies',
                  sorted_by_radius = T)+
  geom_scatterpie(data = finalpie2[finalpie2$Municipality %in% c("Bertem", "Kallo"),], 
                   aes(x=long, y=lat, group=Municipality, r=r), 
                   cols = c("Culex pipiens molestus", "Culex pipiens pipiens", "Culex torrentium", "Aedes japonicus"),
                   color="black", size=0.15,
                   sorted_by_radius = T)+
  geom_scatterpie_legend2(finalpie2$r, x=3, y=49.5, n=3, labeller=function(x) (50*x)^2, linetype = "dotted")+
  scale_fill_viridis_d(begin=0.4, end =0.9, name="")+
  geom_text_repel(data = finalpie2, aes(x=long, y=lat, label=Municipality), size=3, box.padding = 1.5, 
                  min.segment.length = .7, segment.size=.2, seed=123)+
  theme_void()+
  theme(legend.text = element_text(face="italic"),
        legend.key = element_rect(color = NA),
        legend.background = element_rect(color = NA),
        legend.position = c(.15,.35))+
  coord_equal()
p3

ggplot()+
  geom_scatterpie(data = finalpie2, 
                          aes(x=long, y=lat, group=Municipality, r=r), 
                          cols = c("Culex pipiens molestus", "Culex pipiens pipiens", "Culex torrentium", "Aedes japonicus"),
                          color="transparent", size=0.15, legend_name = 'Subpecies',
                          sorted_by_radius = T)+
  #geom_scatterpie(data = finalpie2[!finalpie2$Municipality %in% c("Bertem", "Kallo"),], 
         #        aes(x=long, y=lat, group=Municipality, r=r), 
         #        cols = c("Culex pipiens molestus", "Culex pipiens pipiens", "Culex torrentium", "Aedes japonicus"),
         #        color="transparent", size=0.15, legend_name = 'Subpecies',
         #        sorted_by_radius = T)+
  #geom_scatterpie(data = finalpie2[finalpie2$Municipality %in% c("Bertem", "Kallo"),], 
  #                aes(x=long, y=lat, group=Municipality, r=r), 
  #                cols = c("Culex pipiens molestus", "Culex pipiens pipiens", "Culex torrentium", "Aedes japonicus"),
  #                color="black", size=0.15,
  #                sorted_by_radius = T)+
  geom_scatterpie_legend2(finalpie2$r, x=3, y=49.5, n=3, labeller=function(x) (50*x)^2, linetype = "dotted")+
  scale_fill_viridis_d(begin=0.4, end =0.9, name="")+
  geom_text_repel(data = finalpie2, aes(x=long, y=lat, label=Municipality), size=3, box.padding = 1.5, 
                  min.segment.length = .7, segment.size=.2, seed=123)+
  theme_void()+
  theme(legend.text = element_text(face="italic"),
        legend.key = element_rect(color = NA),
        legend.background = element_rect(color = NA),
        legend.position = c(.15,.35))+
  coord_equal()


finalpie2[!finalpie2$Municipality %in% c("Bertem", "Kallo"),]
finalpie2[finalpie2$Municipality %in% c("Bertem", "Kallo"),]

ggsave("figures/pies.pdf", bg = 'transparent')

?geom_scatterpie

p+geom_point(data = finalpie2, 
                  aes(x=long, y=lat, size=totalcount, color=Municipality), alpha=0.6)+ 
  #scatterpie::geom_scatterpie_legend(finalpie2$totalcount/500, x=3, y=49.75, n=3, labeller=function(x) 500*x)+
  theme_void()

?geom_scatterpie

install.packages("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf()+
  #geom_sf(data = belsf)+
  coord_sf(xlim=c(2.5, 6.5), ylim=c(49.5,51.5))+
  theme_void()

pie = data.frame(
  lon=c(-5.0,-3.5,-5.5,5.0), 
  lat=c(50.0,50.2,50.1,50.5), 
  A=c(0.25,0.75,0,0.25), 
  B=c(0.75,0.10,0,0.75), 
  C=c(0,0.15,1,0), 
  radius=0.05)
pie

finalpie3 <- finalpie2 %>% 
  rowwise() %>% 
  mutate(across(4:7, ~./totalcount))
finalpie3

finalpie3 %>% 
  tidyr::gather(type, value, -long, -lat, -r, -totalcount, -Municipality)

pie.list <- finalpie3[!finalpie3$Municipality %in% c("Bertem", "Kallo"),] %>% 
  tidyr::gather(type, value, -long, -lat, -r, -totalcount, -Municipality) %>%
  tidyr::nest(type, value) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = value, fill = type)) +
                                                        geom_col(color = "transparent",
                                                                 show.legend = FALSE, alpha=1) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_viridis_d(begin=.4, end = .9)+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(r = r*1.5) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = long - r, xmax = long + r,
                                          ymin = lat - r, ymax = lat + r)))

pie.list2 <- finalpie3[finalpie3$Municipality %in% c("Bertem", "Kallo"),] %>% 
  tidyr::gather(type, value, -long, -lat, -r, -totalcount, -Municipality) %>%
  tidyr::nest(type, value) %>%
  
  # make a pie chart from each row, & convert to grob
  mutate(pie.grob = purrr::map(data,
                               function(d) ggplotGrob(ggplot(d, 
                                                             aes(x = 1, y = value, fill = type)) +
                                                        geom_col(color = "transparent", size=.1,
                                                                 show.legend = FALSE, alpha=1) +
                                                        coord_polar(theta = "y") +
                                                        scale_fill_viridis_d(begin=.4, end = .9)+
                                                        theme_void()))) %>%
  
  # convert each grob to an annotation_custom layer. I've also adjusted the radius
  # value to a reasonable size (based on my screen resolutions).
  rowwise() %>%
  mutate(r = r*1.5) %>%
  mutate(subgrob = list(annotation_custom(grob = pie.grob,
                                          xmin = long - r, xmax = long + r,
                                          ymin = lat - r, ymax = lat + r)))

#pie.list

#pie.list$subgrob

p + 
 # Optional. this hides some tiles of the corresponding color scale BEHIND the
 # pie charts, in order to create a legend for them
 geom_tile(data = finalpie3 %>% tidyr::gather(type, value, -long, -lat, -r,-totalcount, -Municipality),
           aes(x = long, y = lat, fill=type), 
           color = "transparent", width = 0.01, height = 0.01, 
           inherit.aes = FALSE) +
  scale_fill_viridis_d(begin=.4, end = .9, name="")+
  pie.list$subgrob+
  pie.list2$subgrob+
  geom_text_repel(data=finalpie3[!finalpie3$Municipality %in% c("Vrasene", "Leuven","Villers-Le-Bouillet",
                                                               "Eupen", "Maasmechelen"),], 
                  aes(x=long, y=lat, label=Municipality), size=3, segment.size=.2,
                  seed=1, direction = "both", box.padding = 1)+
  geom_text_repel(data=finalpie3[finalpie3$Municipality %in% c("Vrasene", "Leuven","Villers-Le-Bouillet",
                                                                "Eupen", "Maasmechelen"),], 
                  aes(x=long, y=lat, label=Municipality), size=3, min.segment.length = .7, segment.size=.2,
                  seed=321, box.padding = 1, segment.colour="transparent", direction = "y")+
  geom_scatterpie_legend2(finalpie3$r, x=3, y=49.75, n=3, labeller=function(x) (50*x)^2, linetype = "dotted")+
  geom_point(data = finalpie2, aes(x=long, y=lat), size=.5)+
  theme_void()+
  theme(legend.position = c(.2,.3),
        legend.text = element_text(face="italic"))

ggsave("figures/Belgium_species.pdf", dpi=300)
#ggsave("Belgium_subspecies_new.pdf", dpi=300)

ggplot()+
  geom_scatterpie_legend2(finalpie3$r, x=3, y=49.75, n=3, labeller=function(x) 500*x, linetype = "dotted")+
  #coord_quickmap()+
  theme_void()


