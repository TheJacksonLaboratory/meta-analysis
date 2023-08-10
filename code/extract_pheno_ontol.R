#extract_pheno_ontol
#' extract phenotype measures associated with ontology term
#'
#' @param url url of MPD. Default: "https://phenome.jax.org/"
#' @param ontology the name ontology term. for example "VT:0010488"
#'
#' @return one data.frame for all phenotype measures associated with ontology term
#'

extract_pheno_ontol <- function(url = "https://phenome.jax.org/",
                                ontology) {
  #if(missing(ontology))
  if (missing(ontology)) {
    stop(paste("The ontology cannot be missing. Please enter the ontology term!"))
  }
  #api
  api <- paste0(url, "api/pheno/measures_by_ontology/", ontology, "?csv=yes")
  #read
  message(paste0(" - Phenotype measures associated with ontology term ", ontology))
  #get phenotype measures associated with ontology term from api
  df = tryCatch(readr::read_csv(api, col_types = cols()),
                HTTPError = function(e) {
                  return(paste("HTTP error: ", trimws(e$message)))
                })

  #get the individual data for each measure in the df.
  dat <- map(df$measnum, safely(function(x){
    #get the data from api
    readr::read_csv(paste0(url, "api/pheno/animalvals/", x, "?csv=yes"), col_types = cols())
  })) %>%
    purrr::map(~.x$result) %>%
    setNames(df$measnum) %>%
    compact() #remove NULL from list
  #rbind
  dat_all <- do.call(rbind, dat) %>%
    dplyr::mutate(strain = gsub(" ", "", strain)) %>% #remove space from strain name
    mutate(strain = case_when(
      str_starts(strain, "CC")    ~ gsub('.{1}$', '', strain), #remove last character "J" from CC name
      TRUE                        ~ as.character(strain)
    ))
  return(dat_all)
}
