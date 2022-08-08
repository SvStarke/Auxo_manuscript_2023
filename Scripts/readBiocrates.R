#' A tidy parsing of a Biocrates XLSX-file
#'
#' @description Reads and parses biocrates data excel tables. Does not only
#' extract numerical values but also he measurment quality status.
#' 
#' @return A list with 5 elements: `BCdat$sample_meas` - a data.table with the
#' samples' meta information. `mets` a data.table with all the metabolites and 
#' metabolic indicator (if any) columns from the input excel. `dat` - The actual
#' numeric results from the metabolite quantification. `dat_format` - Id of the 
#' measurements' quality status. `format_legends` - Legend for the quality status.
#'
#' @author Silvio Waschina (s.waschina@nutrinf.uni-kiel.de)
readBiocrates <- function(file, sheet, incl_val_status = c("Valid","< LLOQ","> ULOQ")) {
  format_colors <- c(Valid    = "FF00CD66", 
                     `< LLOQ` = "FF87CEEB", 
                     `> ULOQ` = "FF6666FF",
                     `< LOD`  = "FF00FFCC",
                     `Warning! ISTD out of Range` = "FFFFFF33",
                     `no intercept: calibration curve asymptotic. Value cannot be specified.` = "FF808080")
  acc_colors <- format_colors[incl_val_status]
  
  # ex 
  #file <- "data/biocrates/2020_034_2020-11-30_UKSH_Feces_Oktober2020.xlsx"
  #sheet <- "Data_Export_2"
  
  #file <- "data/biocrates/2020-034_2020-11-10_UKSH_Serum_Oktober_Ergebnisse_ITEM.xlsx"
  #sheet <- "recoded"
  
  sheet_keep <- sheet
  
  require(tidyxl)
  require(data.table)
  
  formats <- xlsx_formats(file)
  formats <- formats$local$fill$patternFill$fgColor$rgb
  accepted_formats <- which(formats %in% acc_colors)
  
  format_color_ids <- lapply(names(format_colors), function(x) {
    format_id <- which(formats == format_colors[x])
    if(length(format_id) > 0) {
      format_id <- data.table(id = format_id,
                              status = x)
    } else {
      format_id <- data.table(id = integer(0),
                              status = character(0))
    }
    
    return(format_id)
  })
  format_color_ids <- rbindlist(format_color_ids)
  
  #acc_row_tmp <- cells[character %in% incl_val_status, row]
  #accepted_formats <- cells[row %in% acc_row_tmp & col == 1, local_format_id]
  

    
  cells <- xlsx_cells(file)
  cells <- data.table(cells)
  cells <- cells[sheet == sheet_keep]
  
  # Extracting concentration data part
  dat_start_row <- cells[character == "Measurement Time", row] + 1
  dat_start_col <- cells[character == "Measurement Time", col] + 1
  #dat_end_col <- cells[grepl("^Choline$", character), col]
  dat_end_col <- max(cells$col)
  
  dat <- copy(cells[row >= dat_start_row & col >= dat_start_col & col <= dat_end_col])
  dat$conc <- NA_real_
  dat[local_format_id %in% accepted_formats, conc := numeric]
  dat <- dcast(dat, row ~ col, value.var = "conc")
  dat <- as.matrix(dat[,2:ncol(dat)])
  
  met_names <- cells[row == dat_start_row -1 & col >= dat_start_col & col <= dat_end_col,
                     character]
  colnames(dat) <- met_names
  
  # Extracting data format (valid, < LOD, < LLOQ, ...)
  dat_format <- copy(cells[row >= dat_start_row & col >= dat_start_col])
  dat_format <- dcast(dat_format, row ~ col, value.var = "local_format_id")
  dat_format <- as.matrix(dat_format[,2:ncol(dat_format)])
  colnames(dat_format) <- met_names
  
  format_color_ids <- format_color_ids[id %in% unique(as.vector(dat_format))]
  if(!all(unique(as.vector(dat_format)) %in% format_color_ids$id)) {
    warning("Not all local formats (i.e. measurment status) were recognised.")
  }
  
  # Extract col data (metabolites)
  met_start_row <- cells[character == "Class" & row < dat_start_row & col == dat_start_col-1, row]
  mets <- data.table(met_name  = met_names,
                     met_class = cells[row == met_start_row & col >= dat_start_col, character])
  
  # Extract Sample names and their measurement date
  sample_start_col <- cells[row == dat_start_row-1 & !is.na(character), min(col)]
  
  sample_meas <- data.table(row = cells[col == sample_start_col & row >= dat_start_row, row])
  for(i in sample_start_col:(dat_start_col-1)) {
    tmp_type <- cells[col == i & row == dat_start_row, data_type]
    tmp_name <- cells[col == i & row == dat_start_row-1, character]
    sample_meas[[tmp_name]] <- cells[col == i & row >= dat_start_row][[tmp_type]]
  }
  
  #rownames(dat) <- rownames(dat_format) <- sample_meas$sample
  
  return(list(sample_meas = sample_meas,
              mets        = mets,
              dat         = dat,
              dat_format  = dat_format,
              format_legend = format_color_ids))
}



