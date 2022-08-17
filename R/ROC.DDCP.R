##' ROC curve areas and variances estimated as per Delong et al, 1988.
##'
##' This is a utility function that users will probably not want to
##' invoke directly.
##' @title ROC.DDCP
##' @param response logical vector
##' @param indep matrix or \code{data.frame} of
##'     \code{length(response)} rows containing numeric values.
##' @seealso DeLong, Elizabeth R., David M. DeLong, and Daniel
##'     L. Clarke-Pearson. \dQuote{Comparing the areas under two or
##'     more correlated receiver operating characteristic curves: a
##'     nonparametric approach.}  Biometrics (1988): 837-845.
##' @return \code{list} with element \code{theta} giving the areas
##'     under the ROC curve for the variables in \code{indep} and
##'     element \code{var} giving the variance-covariance matrix of
##'     the estimates.
##' @author Charles Berry
ROC.DDCP<-
  function(response, indep)
{
    if(!is.logical(response))
        stop("response must be a boolean variable")
    if(any(is.na(response)))
        stop("NA's not allowed in response")
    indep <- as.matrix(indep)
    x <- indep[!response,  , drop = FALSE ]
    y <- indep[response,  , drop = FALSE ]
    m <- sum(!response)
    n <- sum(response)
    x.ranks <- apply(x,2,rank)
    y.ranks <- apply(y,2,rank)
    is.x <- rep( c(TRUE,FALSE), c( m, n ))
    xy.ranks <- lapply(1:ncol(indep),
                       function(z) split(rank(c(x[,z],y[,z])), is.x ))
    V.10 <- 1 - sapply( 1:ncol(indep),
                       function(z) 2*(xy.ranks[[z]][[2]] - x.ranks[, z ])/ n )/2
    V.01 <- sapply( 1:ncol(indep),
                   function(z) 2*(xy.ranks[[z]][[1]]-y.ranks[, z])/ m )/2
    theta.hat <- colMeans(V.10)
    S.10 <- stats::var(V.10)
    S.01 <- stats::var(V.01)
    ## ensure S is posdef
    S <- S.10/nrow(x) + S.01/nrow(y) + 1e-10*diag(ncol(S.10))
    list(theta = theta.hat, var = S)
}
