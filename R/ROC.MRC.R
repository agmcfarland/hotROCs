##' ROC curve areas, variances, and p-values for datasets in which
##' each case has a collection of matching controls.
##'
##' When data are collected under a scheme in which matched controls
##' are used to allow for biased selection and every case has \code{k}
##' controls, the ROC area is estimated as the average fraction of
##' controls whose variable value is less than that of its reference
##' case. Variances are computed from the usual sample variance
##' divided by the number of cases.
##' @title ROC.MRC - ROC areas for matched control data
##' @param response logical vector identifying cases or character
##'     vector or factor vector with \dQuote{insertion} marking the
##'     cases. Alternatively, a factor or character vector with
##'     elements matching \code{"insertion"} denoting the cases.
##' @param stratum a vector of \code{length(response)} elements whose
##'     unique values correspond to cases and their matched
##'     controls. There must be exactly one case and \code{k} controls
##'     in every stratum.
##' @param variables numeric a \code{matrix} or \code{data.frame} with
##'     \code{length(response)} rows and two or more columns.
##' @param origin \code{NULL} or a vector of \code{length(response)}
##'     elements to identify different data sources.
##' @param origin.levels optional character vector of origin levels
##' @param ragged.OK logical - \code{TRUE}, if differing numbers of
##'     MRCs of any one origin are acceptable.
##' @return \code{list} with elements \sQuote{ROC} giving a matrix of
##'     ROC curve areas, \sQuote{var} giving a list of variance
##'     matrices, and \sQuote{pvalues} giving a list of matrices
##'     containing pvalues for various contrasts.
##' @export
##' @examples
##' case <- rep(rep(c(TRUE,FALSE),c(1,3)),1000)
##' group <- rep(1:1000,each=4)
##' var1 <- case + rnorm(4000)
##' var2 <- case + rexp(4000)
##' var3 <- case + runif(4000,-1,1)
##' origin <- rep(factor(letters[1:4]),each=1000)
##' ROC.MRC(case,group,cbind(var1,var2,var3),origin)
##' @author Charles Berry
ROC.MRC <-
    function(response,stratum,variables,origin=NULL,origin.levels=NULL,
             ragged.OK=TRUE)
{

    if (!is.logical(response)) response <- response == "insertion"
    if (is.null(origin)) origin <- rep(1,length(response))
    if (is.null(origin.levels))
        origin.levels <- as.character(unique(origin))
    stopifnot(is.logical(response))
    nvars <- ncol(variables)
    origin.levels <-
        if (is.factor(origin))
            levels(origin)
        else
            unique(as.character(origin))
    phi.fun <-
        function(x)
    {
        ok.rows <- origin==x
        stratum <- stratum[ok.rows]
        response <- response[ok.rows]
        variables <- variables[ok.rows,]
        rstrata <- stratum[response]
        mstrata <- stratum[!response]
        stopifnot(all(!duplicated(rstrata)))
        nsites <- length(rstrata)
        cmtab <- table(stratum,response)
        nMRCs <- cmtab[1,1]
        if (any(cmtab[,1]!=nMRCs)){
            if (ragged.OK){
                morder <- order(mstrata)
                mindex <- split(which(!response)[morder],
                                mstrata[morder])
                lenMRCs <- lengths(mindex)
                nMRCs <- max(lenMRCs)
                mindex[lenMRCs<nMRCs] <-
                    lapply(mindex[lenMRCs<nMRCs],
                           function(x) c(x,rep(NA,nMRCs-length(x))))
                mindex <- unlist(mindex,use.names=FALSE)
            } else {
                stop("Differing Numbers of MRCs not allowed.")
            }
        } else {
            mindex <- which(!response)[order(mstrata)]
        }
        
        if (any(cmtab[,2]!=1)) stop("MRCs with no matching Integration Site.")
        phi <- 
            sweep(
                array(variables[mindex,],
                      dim=c(nMRCs, nsites, nvars),
                      dimnames=list(NULL,NULL,colnames(variables))),2:3,
                variables[which(response)[order(rstrata)],],
                function(x,y) (sign(y-x)+1)/2)
        colMeans(phi,na.rm=TRUE)
    }

    phi.list <-
        sapply(origin.levels,phi.fun,simplify=FALSE)
    rocz <- sapply(phi.list,colMeans)
    ## inflate the variance by 1e-10 to avoid diff/0.0 in Stats 
    roczVar <- lapply(phi.list,function(x) 1e-10*diag(ncol(x)) + var(x)/nrow(x))
    nullStats <- (rocz-0.5)^2/sapply(roczVar,diag)
    nullPvals <- pchisq(nullStats,df=1,lower.tail=FALSE)
    variableDiffs <-
        do.call(rbind,combn(nvars,2,function(x) rocz[x[1],]-rocz[x[2],],
                simplify=FALSE))
    variableDVars <-
        sapply(roczVar,
	       function(x) combn(nvars,2,
                                 function(y) sum(x[y,y]*c(1,-1,-1,1))))
    variableDStats <- variableDiffs^2/variableDVars
    variablePvals <- pchisq(variableDStats,df=1,lower.tail=FALSE)
    attr(variablePvals,"whichRow") <- combn(nvars,2)
    if (length(origin.levels)>1){
        originDVars <- combn(roczVar,2,function(x) diag(x[[1]])+diag(x[[2]]))
        originDiffs <- combn(origin.levels,2,function(x) rocz[,x[1]]-rocz[,x[2]])
        originStats <- originDiffs^2/originDVars
        originPvals <- pchisq(originStats,df=1,lower.tail=FALSE)
        attr(originPvals,"whichCol") <- combn(length(origin.levels),2)
        rownames(originPvals) <- rownames(rocz)
    } else {
        originPvals <- matrix(NA,nrow=nvars,ncol=1L)
    }
    pvals <- list(op=originPvals,vp=variablePvals,np=nullPvals)
    list(ROC=rocz,var=roczVar,pvalues=pvals)
}
