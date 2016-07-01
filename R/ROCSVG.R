##' Create an interactive ROC curve heatmap
##'
##'  A file \code{file.path(svg.file.base, "main.svg" ) } is created
##'  which can be viewed with and SVG viewer like the FireFox
##'  browser. There are many other files linked to that file that
##'  allow the user to inspect significance test results.
##' @title ROC curve heatmaps
##' @param roc.res.list An object such as returned by
##'     \code{\link{ROC.ORC}} or \code{\link{ROC.MRC}}
##' @param svg.file.base The name of the directory for in which to
##'     store the resulting files
##' @return \code{NULL} -- the function is run for its side effects.
##' @importFrom grDevices dev.off
##' @importFrom graphics axis box image layout mtext par text title
##' @importFrom stats pchisq var
##' @importFrom utils combn
##' @importFrom RSVGTipsDevice devSVGTips setSVGShapeURL setSVGShapeToolTip
##' @importFrom colorspace diverge_hcl
##' @export
##' @author Charles Berry
ROCSVG <-
    function(roc.res.list,svg.file.base='roc')
{
    ## Purpose: Produce ROC heatmap with dynamic  p-value display
    ## ----------------------------------------------------------------------
    ## Arguments: roc.res - object of ROC.strata()
    ##            file - where to save results
    ## ----------------------------------------------------------------------
    ## Author: CCB, Date:  7 Oct 2008, 13:14

    roc.res <- roc.res.list$ROC
    
    ## colormap:

    
    dcol <- diverge_hcl(21, c = c(100, 0), l = c(50, 90), power = 1)


    ## main svg file:

    ## mainfile <- paste(svg.file.base,"main.svg",sep='-')

    mainfile <- file.path(svg.file.base,"main.svg")

    ## open up left margin

    bigmarg <- c(3,10,4,2)+.1

      ### using a subdir for most svg files and putting the unadorned
    ### version above it causes some headaches with keeping button
    ### relative URLs pointing to the right place. hence the use of
    ### strip.dirname and mk.image(use.base=...)

    
    strip.dirname <- function(buttons) lapply(buttons,function(x) {if(dirname(x['URL'])==svg.file.base) x['URL'] <- basename(x['URL']);x})
    
    ## button for the right margin

    null.URL <- file.path(svg.file.base,"H50.svg")
    relative.null.URL <- basename(null.URL)
    button.list <- list(
                        c(text="<Show Plain Heatmap>",URL=mainfile,tiptitle="Click to:",tipdesc="Clear Annotations"),
                        c(text="<Compare to Area == 0.50>",URL=relative.null.URL,tiptitle="Test Each Area",tipdesc="vs Chance Discrimination")
                        )
    
    ## transform to make image 'look right'

    troc.res <- t(roc.res)
    roc.rows <- nrow(roc.res)
    roc.cols <- ncol(roc.res)

    ### generic call to image

    mk.image <- function(title.=NULL,use.base=TRUE){ # TRUE for subdir images
        image(1:roc.cols,1:roc.rows,troc.res,zlim=c(0,1),axes=FALSE,col=dcol,
              xlab="",ylab="")
            if (!is.null(title.)) title(title.,line=4)
        box()
        sapply(1:roc.rows,
               function(x) {
                    setSVGShapeToolTip(title="Compare rows to:",desc=rownames(roc.res)[x])
                   ru <- if (use.base) basename(rowURL[x]) else rowURL[x]
                   setSVGShapeURL( ru )
                   mtext(rownames(roc.res)[x],side=2,adj=1,at=x,las=1,line=1)
               })

        sapply(1:roc.cols,
               function(x) {
                    setSVGShapeToolTip(title="Compare columns to:",desc=colnames(roc.res)[x])
                   cu <- if (use.base) basename(colURL[x]) else colURL[x]
                   setSVGShapeURL( cu )
                   mtext(colnames(roc.res)[x],side=3,adj=0.5,at=x,las=1,line=1)
               })
    }

    mk.buttons <- function(blist,descend){

      ## button list is a list with each element a vector with
      ## elements named "text", "URL", "tiptitle","tipdesc"

      n.button <- length(blist)
      at.pos <- if (roc.rows<2*n.button) {
        seq(0,roc.rows, length=n.button+2)[-c(1,n.button+2)]
      } else {
        seq(roc.rows-0.5,by=-2,length=n.button)
      }
      for (i in 1:n.button){
        button <- blist[[i]]
        if (!is.null(button['tiptitle'])&nchar(button['tiptitle'])!=0){
          setSVGShapeToolTip(title=button['tiptitle'],desc=button['tipdesc'])
        }
        setSVGShapeURL( button['URL'] ) 
        mtext(button['text'],side=4,adj=0,at=at.pos[i],las=1,line=2,
              col='blue',cex=0.8)
      }
    }

    
    mk.stripe <- function(stripe.ht=1){
      new.mai <- par()$mai
      new.mai[1] <- stripe.ht * 0.65
      new.mai[3] <- 0.25
      par(mai=new.mai)
      image(seq(0,1,length=length(dcol)),1,matrix(seq(0,1,length=length(dcol))),
            zlim=c(0,1),axes=FALSE,col=dcol,
            xlab="",ylab="",main="Color Key",cex.main=0.8)
      axis(1)

    }
    
    rowURL <- file.path(svg.file.base,paste('row',1:roc.rows,'svg',sep='.'))
    colURL <- file.path(svg.file.base,paste('col',1:roc.cols,'svg',sep='.'))

    ### functions to add text overlays:
    
    do.overlay <- function(ov) do.call(text,mk.overlay(ov)) 
    mk.overlay <- function(ovmat) # ovmat is a character matrix
      c(
        subset( data.frame(x=as.vector(col(ovmat)),y=as.vector(row(ovmat)),labels=as.vector(ovmat)),
               labels!=""),
        adj=0.5)

    mk.stars <- function(pvmat){
      x <- array("",dim(pvmat))
      x[] <- as.character( cut( pvmat, c(0, 0.001, 0.01, 0.05, 1), c("***","**","*",""), include.lowest=TRUE))
      x
    }

### here set up the overlays:

    opvals <- roc.res.list$pvalues$op
    vpvals <- roc.res.list$pvalues$vp
    nullpvals <- roc.res.list$pvalues$np
    
    opstars <- mk.stars( opvals )
    vpstars <- mk.stars( vpvals )
    nullpstars <- mk.stars( nullpvals )

    matchCol <- function(x,indx) x[2:1,][x==indx] 
    isCol <- function(x,indx) colSums(x==indx)==1

    omasks <- lapply(1:roc.cols,function(x){
      res <- array("",dim(roc.res))
      res[,x] <- "--"
      wc <- attr(opvals,'whichCol')
      res[,matchCol( wc, x )] <- opstars[ , isCol( wc,x)  ]
      res
    })
 
     vmasks <- lapply(1:roc.rows,function(x){
      res <- array("",dim(roc.res))
      res[x,] <- "|"
      wr <- attr(vpvals,'whichRow')
      res[matchCol(wr,x),] <- vpstars[  isCol(wr,x), ]
      res
    })


    ## draw the images:
    
    ## base image

    ## allow enough space for row labels to fit
    max.rowchar <- max(nchar(rownames(roc.res)))+10
    char.inch <- c(0.10,0.15) ### this is a guess - cin seems not to work for SVGDevice
    right.padding <- 1.0 + char.inch[1]*max(sapply(button.list,function(x) nchar(x['text'])))
    par.args <- list(xpd=NA,mai=c(4,max.rowchar,8,2)*rep(rev(char.inch),2)+c(0,0,0,right.padding))
    
    ## be sure columns are wide enough

    max.colchar <- max(nchar(colnames(roc.res)))+1
    min.width <- max(6,max.colchar) * char.inch[1] * roc.cols

    ## be sure to make enough room for text or it will not display
    ## properly

    vpad <- 0.025*roc.rows
    dev.width <- ceiling(min.width + sum(par.args$mai[c(2,4)]))
    min.height <- ceiling(  char.inch[2]*roc.rows )
    dev.height <- vpad+min.height + sum(par.args$mai[c(1,3)])

    layout.args <- list( matrix(1:2,ncol=1), heights=c(dev.height,1) )

    dir.create( svg.file.base, showWarnings = FALSE )

    
    devSVGTips(file = mainfile, width=dev.width,height=dev.height+1)
    do.call(layout, layout.args)
    do.call(par, par.args )
    mk.image("ROC Curve Areas")
    mk.buttons( button.list[ -1 ] )
    mk.stripe()
    dev.off()

    ## populate subdir with overlays:
    
    for (i in seq(along=vmasks)){
      devSVGTips(file = rowURL[i],width=dev.width,height=dev.height+1)
      do.call(layout, layout.args)
      do.call( par, par.args )
      mk.image("Rows Compared",use.base=TRUE)
      mk.buttons( strip.dirname(button.list) )
      do.overlay( vmasks[[i]])
      mk.stripe()
      dev.off()
    }

    for (i in seq(along=omasks)){
      devSVGTips(file = colURL[i],width=dev.width,height=dev.height+1)
      do.call(layout, layout.args)
      do.call( par, par.args )
      mk.image("Columns Compared",use.base=TRUE)
      mk.buttons( strip.dirname(button.list) )
      do.overlay( omasks[[i]])
      mk.stripe()
      dev.off()
    }


    devSVGTips(file = null.URL,width=dev.width,height=dev.height+1)
    do.call(layout, layout.args)
    do.call( par, par.args )
    mk.image("Compare to Chance Discrimination",use.base=TRUE)
    mk.buttons( strip.dirname( button.list[ -grep(relative.null.URL,sapply(button.list,'[','URL')) ] ))
    do.overlay( nullpstars )
    mk.stripe()
    dev.off()
  
    
    invisible()

}
