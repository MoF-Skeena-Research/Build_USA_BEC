# 
# dat = X1
# alpha = .5
# numIDvars = 2
# curr = "BGdh_OR"
# removeOutlier <- function(dat, alpha,numIDvars){
#   out <- foreach(curr = unique(as.character(dat$BGC)), .combine = rbind) %do% {
#     temp <- as.data.frame(dat[dat$BGC == curr,])
#     md <- tryCatch(mahalanobis(temp[,-c(1:numIDvars)],
#                                center = colMeans(temp[,-c(1:numIDvars)]),
#                                cov = cov(temp[,-c(1:numIDvars)])), error = function(e) e)
#     if(!inherits(md,"error")){
#       ctf <- qchisq(1-alpha, df = ncol(temp)-1)
#       outl <- which(md > ctf)
#       cat("Removing", length(outl), "outliers from",curr, "; ")
#       if(length(outl) > 0){
#         temp <- temp[-outl,]
#       }
#     }
#     temp
#     
#   }
#   return(out)
# }

removeOutlier <- function(dat, alpha, numIDvars){
  
  library(dplyr)
  
  id_cols  <- names(dat)[1:numIDvars]
  var_cols <- names(dat)[-(1:numIDvars)]
  
  groups <- dat %>% group_split(BGC)
  
  out_list <- lapply(groups, function(g){
    
    curr <- unique(g$BGC)
    temp <- as.data.frame(g)
    
    # Too few rows → skip
    if(nrow(temp) < 5){
      message("Skipping ", curr, " (too few rows)")
      return(temp)
    }
    
    X <- temp[, var_cols, drop = FALSE]
    
    # Compute covariance
    covmat <- cov(X)
    
    # Determine rank
    qr_obj <- qr(covmat)
    r <- qr_obj$rank
    
    # If full rank, use all variables
    if(r == ncol(X)){
      keep <- seq_len(ncol(X))
    } else {
      # Drop collinear variables
      keep <- qr_obj$pivot[1:r]
      message("Dropping ", ncol(X) - r, " collinear variables for ", curr)
    }
    
    Xr <- X[, keep, drop = FALSE]
    
    # Compute Mahalanobis distances in reduced space
    center <- colMeans(Xr)
    covr   <- cov(Xr)
    
    md <- mahalanobis(Xr, center, covr)
    
    # Threshold
    ctf <- qchisq(1 - alpha, df = ncol(Xr))
    
    outl <- which(md > ctf)
    message("Removing ", length(outl), " outliers from ", curr)
    
    if(length(outl) > 0){
      temp <- temp[-outl, ]
    }
    
    temp
  })
  
  bind_rows(out_list)
}
