#' Print percentage
#' @param ind a vector of for loop iterator
#' @param tot a vector of for loop length
#' @return on screen output of percentage
#' @description From github.com/wangtengyao/putils
printPercentage <- function (ind, tot){
  ind <- as.vector(ind); tot <- as.vector(tot)
  if ((length(tot) > 1) & (length(ind) == 1)) {ind <- match(ind, tot); tot <- length(tot)}
  len <- length(ind)
  contrib <- rep(1,len)
  if (len > 1) {
    for (i in (len-1):1) contrib[i] <- contrib[i+1] * tot[i+1]
  }
  grand_tot <- contrib[1] * tot[1]
  count <- (sum(contrib * (ind - 1)) + 1)
  out <- ""
  if (sum(ind-1)>0) out <- paste0(rep("\b", nchar(round((count-1)/grand_tot * 100))+1), collapse = "")
  out <- paste0(out, round(count/grand_tot*100), "%")
  if (identical(ind, tot)) out <- paste0(out, '\n')
  cat(out)
  return(NULL)
}

