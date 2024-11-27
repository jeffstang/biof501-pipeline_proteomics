library(limma)

create_design <- function (treatment, obs_id = NULL, eset = NULL) {
    if (is.null(obs_id)) {
    design <- model.matrix(~treatment)
  }
  else {
    design <- model.matrix(~treatment)
    dupcor <- duplicateCorrelation(eset, design = model.matrix(~treatment), 
      block = obs_id)
    message("Within-block correlation: ", dupcor$consensus.correlation)
  }
  return(list(design = design, dupcor = if (!is.null(eset) && 
    !is.null(obs_id)) dupcor else NULL))
}