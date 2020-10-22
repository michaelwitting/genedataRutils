#' @title Writing template for Metabolights sample file
#'
#' @description
#'
#' `createMtblsSampleFile` Exports data read from a .gda file to a Metabolights 
#'    sample file. A study ID is required to generate consistent files. To be 
#'    used this function requires the following column annotations to be present:
#'    "Source [C]", "Sample [C]", "Assay [C]", "File [C]", "Sample Type [C]"
#'
#' @param x `list` List with data read from .gda file
#' @param study_id `character` Study ID from Metabolights or place holder 
#'     (default is MTBLSXXXX)
#' @param study_factors `character` Potentially present study factors
#' @param metadata `list` Additional metadata, such a organism etc.
#'#'
#' @author Michael Witting
#' 
#' @export
#'
#' @examples
#'
#'
createMtblsSampleFile <- function(x, study_id = "MTBLSXXXX", study_factors = NA,
                                  metadata = list("organism_characteristics" = NA,
                                                  "organism_source_ref" = NA,
                                                  "organism_term_access" = NA,
                                                  "organism_part_characteristics" = NA,
                                                  "organism_part_source_ref" = NA,
                                                  "organism_part_term_access" = NA,
                                                  "variant_characteristics" = NA,
                                                  "variant_source_ref" = NA,
                                                  "variant_term_access" = NA,
                                                  "protocol_ref" = NA)) {
  
  # sanity checks
  if(!length(x) == 3) {
    
    stop("Input is not of length 3. Sure it contains data from a .gda file?")
    
  }
  
  # check if column annotation contains correct values
  col_anno <- getColAnno(x)
  
  if(!all(c("Source [C]", "Sample [C]", "Assay [C]", "File [C]", "Sample Type [C]") %in% colnames(col_anno))) {
    
    stop("Missing one or more required column annotations")
    
  }
  
  if(!is.na(study_factors) && !all(study_factors %in% colnames(col_anno))) {
    
    stop("Study factors are not present in column annotation")
    
  }
  
  # add fixed part
  sample_df <- data.frame("source_name" = col_anno$`Source [C]`,
                          "organism_characteristics" = rep(metadata$organism_characteristics, nrow(col_anno)),
                          "organism_source_ref" = rep(metadata$organism_source_ref, nrow(col_anno)),
                          "organism_term_access" = rep(metadata$organism_term_access, nrow(col_anno)),
                          "organism_part_characteristics" = rep(metadata$organism_part_characteristics, nrow(col_anno)),
                          "organism_part_source_ref" = rep(metadata$organism_part_source_ref, nrow(col_anno)),
                          "organism_part_term_access" = rep(metadata$organism_part_term_access, nrow(col_anno)),
                          "variant_characteristics" = rep(metadata$variant_characteristics, nrow(col_anno)),
                          "variant_source_ref" = rep(metadata$variant_source_ref, nrow(col_anno)),
                          "variant_term_access" = rep(metadata$variant_term_access, nrow(col_anno)),
                          "sample_type" = col_anno$`Sample Type [C]`,
                          "sample_type_ref" = NA_character_,
                          "sample_type_term_access" = NA_character_,
                          "protocol_ref" = rep(metadata$protocol_ref, nrow(col_anno)),
                          "sample_name" = col_anno$`Sample [C]`)
  
  # add variable part
  if(!is.na(study_factors)) {
    
    for(study_factor in study_factors) {
      
      # remove type description from Genedata
      study_factor_new <- gsub("\\[C\\]|\\[N\\]", "", study_factor)
      
      # add to df
      sample_df[study_factor_new] <- col_anno[study_factor]
      
    }
  }
  
  # create header for file (fixed section)
  header_fixed_fields <- c("Source Name",
                           "Characteristics[Organism]",
                           "Term Source REF",
                           "Term Accession Number",
                           "Characteristics[Organism part]",
                           "Term Source REF",
                           "Term Accession Number",
                           "Characteristics[Variant]",
                           "Term Source REF",
                           "Term Accession Number",
                           "Characteristics[Sample type]",
                           "Term Source REF",
                           "Term Accession Number",
                           "Protocol REF",
                           "Sample Name")
  
  # create header for file (variable section)
  header_variable_fields <- c()
  
  if(!is.na(study_factors)) {
    
    for(study_factor in study_factors) {
      
      header_variable_fields <- c(header_variable_fields,
                                  paste0("Factor Value[",
                                         gsub("\\[C\\]|\\[N\\]", "", study_factor),
                                         "]"),
                                  "Term Source REF",
                                  "Term Accession Number")
      
    }
  }
  
  # combined headers
  header <- c(header_fixed_fields, header_variable_fields)
  
  file <- paste0("s_", study_id, ".txt")
  
  # write file
  .write_sample_file(sample_df, header, file)
  
}

#' helper function
#' @noRd
.write_sample_file <- function(sample_df, header, file) {
  
  cat(paste(header, collapse = "\t"), '\n', file = file)
  write.table(sample_df, file, append = T, sep = "\t", col.names = FALSE, row.names = FALSE)
  
}

#' #' @title Writing template for Metabolights sample file
#' #'
#' #' @description
#' #'
#' #' `createMtblsSampleFile` Exports data read from a .gda file to a Metabolights 
#' #'    sample file. A study ID is required to generate consistent files. To be 
#' #'    used this function requires the following column annotations to be present:
#' #'    "Source [C]", "Sample [C]", "Assay [C]", "File [C]", "Sample Type [C]"
#' #'
#' #' @param x `list` List with data read from .gda file
#' #' @param study_id `character` Study ID from Metabolights or place holder 
#' #'     (default is MTBLSXXXX)
#' #' @param study_factors `character` Potentially present study factors
#' #' @param metadata `list` Additional metadata, such a organism etc.
#' #'#'
#' #' @author Michael Witting
#' #' 
#' #' @export
#' #'
#' #' @examples
#' #'
#' #'
#' create_assay_file <- function(x, study_id) {
#'   
#' }

#' @title Writing template for Metabolights metabolite file
#'
#' @description
#'
#' `createMtblsMetaboliteFile` Exports data read from a .gda file to a Metabolights 
#'    MAF file. A study ID is required to generate consistent files. To be 
#'    used this function requires the following column annotations to be present:
#'    "Source [C]", "Sample [C]", "Assay [C]", "File [C]", "Sample Type [C]"
#'
#' @param x `list` List with data read from .gda file
#' @param study_id `character` Study ID from Metabolights or place holder 
#'     (default is MTBLSXXXX)
#'
#' @author Michael Witting
#' 
#' @export
#'
#' @examples
#'
#'
createMtblsMetaboliteFile <- function(x, study_id) {
  
  # sanity checks
  if(!length(x) == 3) {
    
    stop("Input is not of length 3. Sure it contains data from a .gda file?")
    
  }
  
  # check if column annotation contains correct values
  row_anno <- getRowAnno(x)
  ms_data <- getMsData(x)
  
  # order both dataframes to have same order
  row_anno <- row_anno[order(row.names(row_anno)),]
  ms_data <- ms_data[order(row.names(ms_data)),]
  

  if(!all(row.names(row_anno) == row.names(ms_data))) {
    
    stop("Row annotations and MS data ist not mathcing")
    
  }

  # add fixed part
  metabolite_df <- data.frame("database_identifier" = if(is.null(row_anno$DBID)) {NA_character_} else {row_anno$DBID},
                              "chemical_formula" = if(is.null(row_anno$Formula)) {NA_character_} else {row_anno$Formula},
                              "smiles" = if(is.null(row_anno$SMILES)) {NA_character_} else {row_anno$SMILES},
                              "inchi" = if(is.null(row_anno$InChI)) {NA_character_} else {row_anno$InChI},
                              "metabolite_identification" = if(is.null(row_anno$Name)) {NA_character_} else {row_anno$Name},
                              "mass_to_charge" = row_anno$`m/z`,
                              "fragmentation" = NA_character_,
                              "modifications" = NA_character_,
                              "charge" = NA_character_,
                              "retention_time" = row_anno$RT,
                              "taxid" = NA_character_,
                              "species" = NA_character_,
                              "database" = NA_character_,
                              "database_version" = NA_character_,
                              "reliability" = if(is.null(row_anno$MSI)) {NA_character_} else {row_anno$MSI},
                              "uri" = NA_character_,
                              "search_engine" = NA_character_,
                              "search_engine_score" = NA_character_,
                              "smallmolecule_abundance_sub" = NA_character_,
                              "smallmolecule_abundance_stdev_sub" = NA_character_,
                              "smallmolecule_abundance_std_error_sub" = NA_character_,
                              "id" = NA_character_)
  
  metabolite_df <- cbind.data.frame(metabolite_df, ms_data)
  
  write.table(metabolite_df, file = paste0("m_", study_id), row.names = FALSE)
  
}
