RefTM_postprocess = function(object,k1,erase.BF = TRUE){
  if(class(object) == "list"){
    theta = object$theta
    theta[,1:k1] =theta[,1:k1]*(sum(colSums(theta)[-(1:k1)])/sum(colSums(theta)[1:k1]))
  }else if(class(object) == "STM" & erase.BF == TRUE){
    eta = cbind(object$eta - t(object$mu$mu),rep(0,dim(object$eta)[1]))
    theta = exp(eta)/rowSums(exp(eta))
  }else{
    theta = object$theta
  }
  return(theta)
}
