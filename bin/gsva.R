library(GSVA)
library(limma)

gsva_limma <- function (eset, design, db, obs_id = NULL, correlation = NULL) {
  length_receptors <- sapply(db, length)
  gsvapar <- gsvaParam(eset, db, minSize = min(length_receptors), maxDiff = TRUE)
  gsva_eset <- gsva(gsvapar)
  if (!is.null(obs_id)) {
    fit <- eBayes(lmFit(gsva_eset, design, block = obs_id, correlation = correlation))
    message("fitting model with paired samples.")
  }
  else {
    fit <- eBayes(lmFit(gsva_eset, design))
    message("fitting model without paired sample consideration.")
  }
  top <- topTable(fit, coef = 2, number = nrow(fit))
  return(top)
}