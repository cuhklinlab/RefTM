#Optimizing beta
opt.beta <- function(beta.ss, kappa, settings,refBeta = NULL) {
  #if its standard lda just row normalize
  if(is.null(kappa)) return(list(beta=list(beta.ss[[1]]/rowSums(beta.ss[[1]]))))

  #If its a SAGE model use the distributed poissons
  if(settings$tau$mode=="L1") {
    out <- mnreg(beta.ss,settings,refBeta)
  } else {
    out <- stm:::jeffreysKappa(beta.ss, kappa, settings)
  }
  return(out)
}

