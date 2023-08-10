#extract_pheno_measnum
#' extract phenotype measures by measnum
#'
#' @param url url of MPD. Default: "https://phenome.jax.org/"
#' @param measnum the name measure_id. for example "15403"
#'
#' @return one data.frame for measnum
#'
extract_pheno_measnum <- function(url = "https://phenome.jax.org/",
                                  measnum) {
  #if(missing(measnum))
  if (missing(measnum)) {
    stop(paste("The measnum cannot be missing!"))
  }
  #api
  api <- paste0(url, "api/pheno/animalvals/", measnum, "?csv=yes")
  #read
  message(paste0(" - Phenotype measure ", measnum))
  #get phenotype measures associated with ontology term from api
  df = tryCatch(readr::read_csv(api, col_types = cols()),
                HTTPError = function(e) {
                  return(paste("HTTP error: ", trimws(e$message)))
                })
  return(df)
}
