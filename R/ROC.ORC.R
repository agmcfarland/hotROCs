
##' ROC curve areas, variances, and p-values for datasets in which the
##' cases have ordinary random controls (i.e. not matched on any
##' characteristic).  The ROC areas and corresponding variances are
##' estimated by the method of Delong et al.  The data may derive from
##' multiple datasets. Comparisons of different variables in a dataset
##' compare the ROC curve areas relative to the controls.  Comparisons
##' of different datasets uses only the responses from each dataset
##' and do not involve the random controls.  It is the responsiibility
##' of the caller to remove \code{NA} values from the arguments.
##'
##' A matrix of ROC curve areas and variances.
##' @title ROC.ORC - ROC area matrix
##' @param response logical vector identifying cases or character
##'     vector or factor vector with \code{"insertion"} marking the
##'     cases.
##' @param variables a \code{matrix} or \code{data.frame} with
##'     \code{length(response)} rows and two or more columns.
##' @param origin \code{NULL} or a vector of \code{length(response)}
##'     elements to identify different data sources.
##' @param origin.levels optional character vector of origin levels
##' @return \code{list} with elements \sQuote{ROC} giving a matrix of
##'     ROC curve areas, \sQuote{var} giving a list of variance
##'     matrices, and \sQuote{pvalues} giving a list of matrices
##'     containing pvalues for various contrasts.
##' @export
##' @examples
##' case <- rep(rep(c(TRUE,FALSE),c(1,3)),1000)
##' var1 <- case + rnorm(4000)
##' var2 <- case + rexp(4000)
##' var3 <- case + runif(4000,-1,1)
##' origin <- rep(factor(letters[1:4]),each=1000)
##' ROC.ORC(case,cbind(var1,var2,var3),origin)
##' @author Charles Berry
ROC.ORC <-
    function(response,variables,origin=NULL,origin.levels=NULL)
{
    if (any(is.na(variables))){
        res <- colSums(is.na(variables))>0
        stop("NA values not allowed. \nFound in:",
             paste(names(res)[res],collapse="\n\t"))
    }
    stopifnot(all(!is.na(response)))
    stopifnot(all(!is.na(origin)))
    if (!is.logical(response)) response <- response == "insertion"
    if (is.null(origin))
        origin <- rep(1,length(response))
    if (is.null(origin.levels))
        origin.levels <- as.character(unique(origin))
    ## res provides all the ROC areas, variances for comparison to
    ## null, and for variable vs variable
    res <-
        lapply(origin.levels,
	       function(lev){
                   resp <- response[origin==lev]
                   indep <- variables[origin==lev,]
                   ROC.DDCP(resp,indep)
	       })
    ## res2 provides variances for origin to origin comparisons
    pairs <- utils::combn(origin.levels,2)
    res2 <-
        apply(pairs,2,
	      function(x) {
                  ok.rows <- response & (origin %in% x)
                  ROC.DDCP(origin[ok.rows]==x[1],
                           variables[ok.rows,])})
    rocz <-
        do.call(cbind, lapply(res,"[[","theta"))
    colnames(rocz) <- origin.levels
    rownames(rocz) <- colnames(variables)

    nullVars <- sapply(res,function(x) diag(x$var))
    nullStats <- (rocz-0.50)^2/nullVars
    nullPvals <- stats::pchisq(nullStats,df=1L,lower.tail=FALSE)
    ncv <- nrow(rocz)
    variableDVars <-
        sapply(res,
	       function(resElt) utils::combn(ncv,2,
				      function(x) sum(resElt$var[x,x]*c(1,-1,-1,1))))
    variableDiffs <-
        do.call(rbind,
                utils::combn(ncv,2,
		      function(x) rocz[x[1],]-rocz[x[2],],simplify=FALSE))
    variableDStats <- variableDiffs^2/variableDVars
    variablePvals <- stats::pchisq(variableDStats,df=1,lower.tail=FALSE)
    attr(variablePvals,"whichRow") <- utils::combn(ncv,2)

    originDiffs <- sapply(res2,"[[","theta")-0.50
    originVars <- lapply(res2,function(x) diag(x$var))
    originStats <- originDiffs^2 / do.call(cbind,originVars)
    originPvals <- stats::pchisq(originStats,df=1,lower.tail=FALSE)
    rownames(originPvals) <- rownames(rocz)
    attr(originPvals,"whichCol") <- utils::combn(length(origin.levels),2)
    list(ROC=rocz,
         var=
             list(within.origin=lapply(res,"[[","var"),
                  between.origin=lapply(res2,"[[","var")),
         pvalues=list(op=originPvals,vp=variablePvals,
		      np=nullPvals))
}
