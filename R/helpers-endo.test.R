#' Summary of endotest
#'
#' @description Summary function for endo.test
#' @keywords internal
#' @return No return value, called for summary.
#' @export
summary.endotest<- function(object,...){
  endotest <- object
  if (typeof(endotest$VHat)=="list") {
    result<-matrix(NA, ncol=3, nrow=length(endotest$VHat))
    result <- data.frame(result)
    colnames(result)<-c("P-value","Test","Valid IVs")
    rownames(result)<-paste0("MaxClique",1:length(endotest$VHat))
    result[,1] <- round(unlist(endotest$p.value),4)
    result[,2] <- ifelse(unlist(endotest$check),"H0 rejected","H0 not rejected")
    for (i in 1:length(endotest$VHat)) {
      result[i,3] <- paste(endotest$VHat[[i]], collapse = " ")
    }
    VHat.union <- Reduce(union,endotest$VHat)
    invalidIV <- setdiff(endotest$SHat,VHat.union)
    print(result,right=F)
    cat(rep("_", 30), "\n")
    if (length(invalidIV)==0) {
      cat("No invalid IV is detected","\n")
    } else{
      cat("Detected invalid IVs:",paste(invalidIV,collapse = " "),"\n")
    }
    # cat("Test result with significance level",endotest$alpha,"\n")
    # cat(rep("_", 30), "\n")


  } else {
    result<-matrix(NA, ncol=3, nrow=1)
    result <- data.frame(result)
    colnames(result)<-c("P-value","Test","Valid IVs")
    rownames(result) <- ""
    result[,1] <- round(endotest$p.value,4)
    result[,2] <- ifelse(endotest$check,"H0 rejected","H0 not rejected")
    result[,3] <- paste(endotest$VHat, collapse = " ")
    # cat("Test result with significance level",endotest$alpha,"\n")
    # cat(rep("_", 30), "\n")
    invalidIV <- setdiff(endotest$SHat,endotest$VHat)

    print(result,right=F)
    cat(rep("_", 30), "\n")
    if (length(invalidIV)==0) {
      cat("No invalid IV is detected","\n")
    } else{
      cat("Detected invalid IVs:",paste(invalidIV,collapse = " "),"\n")
    }
    # cat("\nValid IVs:", endotest$VHat, "\n");
    # cat(rep("_", 30), "\n")
    # cat("P-value = ",endotest$p.value,"\n")
    # if (endotest$check) {
    #   cat("'H0 : Sigma12 = 0' is rejected at the significance level",endotest$alpha,".\n")
    # } else {
    #   cat("'H0 : Sigma12 = 0' is not rejected at the significance level",endotest$alpha,".\n")
    # }
  }

}
endo.SHat <- function(n, ITT_D, V.gamma, method='OLS', tuning.1st=NULL){
  pz = nrow(V.gamma)
  if(method=="OLS"){
    Tn1 = Tn2 = sqrt(log(n))
  }else{
    Tn1 = Tn2 = max(sqrt(2.01*log(pz)), sqrt(log(n)))
  }
  if(!is.null(tuning.1st)) Tn1 = tuning.1st
  # if(!is.null(tuning.2nd)) Tn2 = tuning.2nd
  ## First Stage
  SHat = (1:pz)[abs(ITT_D) > (Tn1 * sqrt(diag(V.gamma)/n))]

  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak.
            TSHT with these IVs will give misleading CIs, SEs, and p-values.
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  return(SHat)
}
