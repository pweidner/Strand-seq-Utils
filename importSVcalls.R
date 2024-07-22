import.svcalls <- function(data_location, filter = c("lenient", "stringent")) {
  # Grab all lenient df from sample folders
  dirs <- list.dirs(data_location, recursive = TRUE)
  dirs <- grep("mosaiclassifier/sv_calls", dirs, value = TRUE)
  if (filter == "lenient"){
    file.names <- list.files(dirs, pattern="lenient_filterFALSE.tsv", full.names = TRUE, recursive = TRUE)
  } else {
    file.names <- list.files(dirs, pattern="stringent_filterFALSE.tsv", full.names = TRUE, recursive = TRUE)
  }
  df.list <- lapply(file.names, function(file.name) {
    df <- readr::read_tsv(file.name, col_names = TRUE)
    return(df)
  })
  
  df <- dplyr::bind_rows(df.list)
  return(df)
}
