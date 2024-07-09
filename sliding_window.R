sliding_window <- function(start, end, window, slide, chrom, bed = FALSE, file_name = "sliding_windows.bed") {
  
  # Create a sliding windows df
  windows <- data.frame(
    chrom = chrom,
    start = seq(start, end - window + 1, by = slide),
    end = seq(start + window - 1, end, by = slide)
  )
  
  # Ensure windows do not exceed the end position
  windows <- windows[windows$end <= end, ]
  
  # Save to BED file if required
  if (bed) {
    write.table(
      windows[1:3],
      file = file_name,
      quote = FALSE,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE
    )
  }
  
  return(windows)
}
