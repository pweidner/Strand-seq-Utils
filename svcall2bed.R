library(tidyverse)
library(grDevices)

# Load your sv call table of interest
df <- read_table("~/Downloads/ADCR11_stringent_filterTRUE.tsv")

# Extract per cell bed files with preserved RGB colors and write them to individual files as in write.table
for (i in unique(df$cell)) {
  message("Writing file ", paste0("data/proc/",i ,"_svcalls.bed", collapse = ""))
  dt <- df %>%
    filter(cell == i) %>% 
    transmute(chrom=chrom,
              chromStart=start,
              chromEnd=end,
              name=sv_call_name,
              score=0,
              strand = ".",
              thickStart = start,
              thickEnd = end,
              itemRgb = case_when(
                sv_call_name == "none" ~ "190,190,190",
                sv_call_name == "del_h1" ~ "119,170,221",
                sv_call_name == "del_h2" ~ "68,119,170",
                sv_call_name == "del_hom" ~ "17,68,119",
                sv_call_name == "dup_h1" ~ "204,153,187",
                sv_call_name == "dup_h2" ~ "170,68,136",
                sv_call_name == "dup_hom" ~ "119,17,85",
                sv_call_name == "inv_h1" ~ "221,221,119",
                sv_call_name == "inv_h2" ~ "170,170,68",
                sv_call_name == "inv_hom" ~ "119,119,17",
                sv_call_name == "idup_h1" ~ "221,170,119",
                sv_call_name == "idup_h2" ~ "170,119,68",
                sv_call_name == "complex" ~ "119,68,17")
              )
  # Make sure numbers are integers due to strict BED formatting
  dt$chromEnd <- as.integer(dt$chromEnd)
  dt$chromStart <- as.integer(dt$chromStart)
  dt$thickStart <- as.integer(dt$thickStart)
  dt$thickEnd <- as.integer(dt$thickEnd)
  # Change here to adapt where and how to save output BED files
  write.table(dt, paste0("data/proc/",i ,"_svcalls.bed", collapse = ""), col.names = F, sep='\t', quote=FALSE, row.names = F)
  assign(paste(i), dt)
}


# Add Track line to each file in the terminal before catentating for upload:

# for file in *; do
# name_start="track name='"
# name_end="' itemRgb="On"'"
# echo $name_start$file$name_end | cat - $file > temp && mv temp $file
# done

# Catenate files for upload as follows:
cat *svcall* > SAMPLE_svcalls.bed
gzip SAMPLE_svcalls.bed


#######################################################################
# Transform the table into sample level .bed file format (all cells in one track):
df <- df %>% 
  transmute(chrom=chrom,
            chromStart=start,
            chromEnd=end,
            name=sv_call_name,
            score=0,
            strand = ".",
            thickStart = start,
            thickEnd = end,
            itemRgb = case_when(
                sv_call_name == "none" ~ "190,190,190",
                sv_call_name == "del_h1" ~ "119,170,221",
                sv_call_name == "del_h2" ~ "68,119,170",
                sv_call_name == "del_hom" ~ "17,68,119",
                sv_call_name == "dup_h1" ~ "204,153,187",
                sv_call_name == "dup_h2" ~ "170,68,136",
                sv_call_name == "dup_hom" ~ "119,17,85",
                sv_call_name == "inv_h1" ~ "221,221,119",
                sv_call_name == "inv_h2" ~ "170,170,68",
                sv_call_name == "inv_hom" ~ "119,119,17",
                sv_call_name == "idup_h1" ~ "221,170,119",
                sv_call_name == "idup_h2" ~ "170,119,68",
                sv_call_name == "complex" ~ "119,68,17")
            )

# Make sure numbers are integers before exporting the file
df$chromEnd <- as.integer(df$chromEnd)
df$chromStart <- as.integer(df$chromStart)
df$thickStart <- as.integer(df$thickStart)
df$thickEnd <- as.integer(df$thickEnd)

# Export table into bed file which you can upload to the browser
write.table(df, "~/Downloads/test.bed", col.names = F, sep='\t', quote=FALSE, row.names = F)

# Add Track line to merged file in terminal later to enable RGB colors using echo 'track name='cell_sv_calls' itemRgb="On"' | cat - test.bed > temp && mv temp test.bed
# Or add "track name='cell_sv_calls' itemRgb="On"" manually in the UCSC broweser later for catenated files
