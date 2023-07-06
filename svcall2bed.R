library(tidyverse)
library(grDevices)

# Load your sv call table of interest
df <- read_table("~/Downloads/ADCR11_stringent_filterTRUE.tsv")

# Transform the table into .bed file format
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

# Add Track line to each file in terminal later to enable RGB colors using echo 'track name='cell_sv_calls' itemRgb="On"' | cat - test.bed > temp && mv temp test.bed
# Or add "track name='cell_sv_calls' itemRgb="On"" manually in the UCSC broweser later for catenated files
