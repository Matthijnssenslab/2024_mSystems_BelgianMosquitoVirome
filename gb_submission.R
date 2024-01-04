library(tidyverse)
library(anytime)
library(lubridate)

gb_df <- read_tsv("data/gb_submission/gb_submission_fixed.tsv")
bioproject_df <- read_tsv("data/gb_submission/BMV_bioproject_attributes.tsv")
SRA_df <- read_tsv("data/gb_submission/BMV_SRA_accessions.tsv")
clusters_df <- read_tsv("data/gb_submission/gb_sub_clusters.tsv", col_names = c("representative", "contig"))

all_genomes <- clusters_df %>% 
  separate_rows(contig, sep=",") %>% 
  mutate(sample=str_split(contig, "_") %>%
           map_chr(~ ifelse(length(.) > 6, .[7], NA)) %>% 
           str_replace("\\|.*", ""),
         full_length=as.numeric(str_match(representative, "length_([0-9]+)_cov")[, 2]),
         contig_length=as.numeric(str_match(contig, "length_([0-9]+)_cov")[, 2]),
         contig_completeness=contig_length/full_length) %>% 
  filter(contig_completeness >= 0.9) %>% 
  group_by(sample, representative) %>%
  filter(contig_length == max(contig_length)) %>%
  ungroup() %>% 
  select(-full_length, -contig_length, -contig_completeness)

all_genomes_w_meta <- full_join(all_genomes,
          gb_df,
          by=join_by("representative" == "Contig"))

all_genomes_w_meta <- left_join(all_genomes_w_meta, 
          SRA_df, 
          by=join_by("sample" == "sample_name"))

all_genomes_w_meta <- left_join(all_genomes_w_meta, 
          bioproject_df, 
          by=join_by("sample" == "sample_name"))


source_table <- all_genomes_w_meta %>%
  arrange(sample, Name) %>%
  group_by(sample, Name) %>%
  mutate(Isolation_source = first(host),
         Isolate = paste(sample, stri_rand_strings(1, pattern = "[a-zA-Z0-9]", length = 5), sep = "_"),
         Organism = first(`MiUViG name`),
         `Collection date` = ifelse(nchar(collection_date) == 4,
                                    format(anytime(paste(collection_date, "01-01", sep = "-")), "%Y"),
                                    format(anytime(collection_date), "%d-%b-%Y"))) %>%
  ungroup() %>%
  select(contig, Name, sample, Organism, Isolate, Segment, Protein, Metagenomic, Metagenomic_source,  
         Country, `Collection date`, bioproject_accession.x, biosample_accession, accession.x) %>%
  rename(SRA_accession = accession.x,
         bioproject_accession = bioproject_accession.x)

write_tsv(source_table, "data/gb_submission/BMV_cleaned_gb_sub.tsv")

final_source_table <- source_table %>% 
  select(-sample, -Protein) %>% 
  rename(Sequence_ID=contig,
         Collection_date=`Collection date`,
         Bioproject=bioproject_accession,
         Biosample=biosample_accession,
         SRA=SRA_accession,
         Metagenome_source=Metagenomic_source)
write_tsv(final_source_table, "data/gb_submission/BMV_gb_submission.src")
