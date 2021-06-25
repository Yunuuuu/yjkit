# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
#
#       Options:
#       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii)  QN = Quantile normalization of input mixture (default = TRUE)
#       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
#               - note that cell subsets will be scaled by their absolute levels and will not be
#                 represented as fractions (to derive the default output, normalize absolute
#                 levels such that they sum to 1 for each mixture sample)
#               - the sum of all cell subsets in each mixture sample will be added to the ouput
#                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
#       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#               - sig.score = for each mixture sample, define S as the median expression
#                 level of all genes in the signature matrix divided by the median expression
#                 level of all genes in the mixture. Multiple cell subset fractions by S.
#               - no.sumto1 = remove sum to 1 constraint
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


#dependencies
# require(e1071)
# require(parallel)
# require(preprocessCore)

#' a modified CIBERSORT function
#' @description  Robust enumeration of cell subsets from tissue expression
#'   profiles
#' @details By default, CIBERSORT estimates the relative fraction of each cell
#'   type in the signature matrix, such that the sum of all fractions is equal
#'   to 1 for a given mixture sample.
#'
#'   CIBERSORT can also be used to produce a score that quantitatively measures
#'   the overall abundance of each cell type (as described in
#'   \href{https://www.nature.com/articles/nmeth.3337}{Analysis of deconvolution
#'   consistency}). Briefly, the absolute immune fraction score is estimated by
#'   the median expression level of all genes in the signature matrix divided by
#'   the median expression level of all genes in the mixture. more details see
#'   \href{https://cibersort.stanford.edu/manual.php#run}{CIBERSORT Manual}
#'
#' @param mixture_data a \code{data.fram} with gene names consistent with
#'   \code{sig_data} gene names in column 1; Data should be in non-log space.
#'   Note: if maximum expression value is <50; CIBERSORT will assume that data
#'   are in log space, and will anti-log all expression values by 2x. If gene
#'   symbols are redundant, CIBERSORT will choose the one with highest mean
#'   expression across the mixtures. CIBERSORT performs a feature selection and
#'   therefore typically does not use all genes in the signature matrix. It is
#'   generally ok if some genes are missing from the user’s mixture file. If
#'   <50% of signature matrix genes overlap, CIBERSORT will issue a warning.
#'   Normal quantification of RNA-seq like FPKM, and TPM can be used
#' @param sig_data CIBERSORT requires an input data.frame of reference gene
#'   expression signatures, or signature data.frame with gene names in column 1,
#'   for routine analysis. This is stored in a Signature Genes File and consists
#'   of a table with groups of "barcode" genes whose expression values
#'   collectively define a unique gene expression signature for each component
#'   pure cell population that will be used to deconvolute the mixture. if
#'   \code{NULL} uses default \code{LM22} with HUGO gene symbols as gene names
#' @param perm No. permutations; set to >=100 to calculate p-values default:
#'   \code{200}
#' @param quantile_norm Quantile normalization of input mixture default:
#'   \code{TRUE}
#' @param absolute Run CIBERSORT in absolute mode default: \code{FALSE}. note
#'   that cell subsets will be scaled by their absolute levels and will not be
#'   represented as fractions (to derive the default output, normalize absolute
#'   levels such that they sum to 1 for each mixture sample); the sum of all
#'   cell subsets in each mixture sample will be added to the ouput ('Absolute
#'   score'). If LM22 is used, this score will capture total immune content.
#' @param abs_method if absolute is set to \code{TRUE}, abs_method choose
#'   method: 'no.sumto1' or 'sig.score' - sig.score = for each mixture sample,
#'   define S as the median expression level of all genes in the signature
#'   matrix divided by the median expression level of all genes in the mixture.
#'   Multiple cell subset fractions by S. - no.sumto1 = remove sum to 1
#'   constraint
#' @return a data.frame
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @references \itemize{ \item Newman, A., Liu, C., Green, M. et al. Robust
#'   enumeration of cell subsets from tissue expression profiles. Nat Methods
#'   12, 453–457 (2015). \url{https://doi.org/10.1038/nmeth.3337} \item Chen B.,
#'   Khodadoust M.S., Liu C.L., Newman A.M., Alizadeh A.A. (2018) Profiling
#'   Tumor Infiltrating Immune Cells with CIBERSORT. In: von Stechow L. (eds)
#'   Cancer Systems Biology. Methods in Molecular Biology, vol 1711. Humana
#'   Press, New York, NY. \url{https://doi.org/10.1007/978-1-4939-7493-1_12}}
#' @export
run_cibersort <- function(mixture_data, sig_data = NULL,
                         perm = 200, quantile_norm = TRUE, absolute = FALSE,
                         abs_method = 'sig.score'){

    if (!requireNamespace("e1071", quietly = TRUE)){
        stop("e1071 needed for this function to work. Please install it",
             call. = FALSE)
    }

    if (!requireNamespace("preprocessCore", quietly = TRUE)){
        stop("preprocessCore needed for this function to work. Please install it",
             call. = FALSE)
    }

    if(absolute) abs_method <- match.arg( 'no.sumto1', 'sig.score' )

    # read in data

    if (is.null(sig_data)) X <- run_cibersort_lm22 else {
        X <- sig_data
    }
    X <- as.data.frame(X, make.names = FALSE)
    rownames(X) <- X[[1]]
    X <- X[, -1, drop = FALSE]


    Y <- as.data.frame(mixture_data, make.names = FALSE)

    #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
    dups <- dim(Y)[1] - length(unique(Y[[1]]))
    if(dups > 0) {
        warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
        rownames(Y) <- make.unique(Y[[1]], sep = ".")
    } else {
        rownames(Y) <- Y[[1]]
    }

    Y <- Y[, -1, drop = FALSE]

    X <- data.matrix(X)
    Y <- data.matrix(Y)

    #order
    X <- X[order(rownames(X)),]
    Y <- Y[order(rownames(Y)),]

    P <- perm #number of permutations

    #anti-log if max < 50 in mixture file
    if(max(Y) < 50) {Y <- 2^Y}

    #quantile normalization of mixture file
    if(quantile_norm == TRUE){
        tmpc <- colnames(Y)
        tmpr <- rownames(Y)
        Y <- preprocessCore::normalize.quantiles(Y)
        colnames(Y) <- tmpc
        rownames(Y) <- tmpr
    }

    #store original mixtures
    Yorig <- Y
    Ymedian <- max(stats::median(Yorig),1)

    #intersect genes
    Xgns <- row.names(X)
    Ygns <- row.names(Y)
    YintX <- Ygns %in% Xgns
    Y <- Y[YintX,]
    XintY <- Xgns %in% row.names(Y)
    X <- X[XintY,]

    #standardize sig matrix
    X <- (X - mean(X)) / stats::sd(as.vector(X))

    #empirical null distribution of correlation coefficients
    if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}

    header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
    if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))

    output <- matrix()
    itor <- 1
    mixtures <- dim(Y)[2]
    pval <- 9999

    #iterate through mixtures
    while(itor <= mixtures){

        y <- Y[,itor]

        #standardize mixture
        y <- (y - mean(y)) / stats::sd(y)

        #run SVR core algorithm
        result <- CoreAlg(X, y, absolute, abs_method)

        #get results
        w <- result$w
        mix_r <- result$mix_r
        mix_rmse <- result$mix_rmse

        if(absolute && abs_method == 'sig.score') {
            w <- w * stats::median(Y[,itor]) / Ymedian
        }

        #calculate p-value
        if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

        #print output
        out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
        if(absolute) out <- c(out, sum(w))
        if(itor == 1) {output <- out}
        else {output <- rbind(output, out)}

        itor <- itor + 1

    }

    #return matrix object containing all results
    obj <- rbind(header,output)
    obj <- obj[, -1, drop = FALSE]
    obj <- obj[-1, , drop = FALSE]
    obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
    rownames(obj) <- colnames(Y)
    if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
    else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
    obj
}

#Core algorithm
CoreAlg <- function(X, y, absolute, abs_method){

    #try different values of nu
    svn_itor <- 3

    res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model <- e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
    }

    if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)

    nusvm <- rep(0,svn_itor)
    corrv <- rep(0,svn_itor)

    #do cibersort
    t <- 1
    while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- stats::cor(k, y)
        t <- t + 1
    }

    #pick best model
    rmses <- nusvm
    mn <- which.min(rmses)
    model <- out[[mn]]

    #get and normalize coefficients
    q <- t(model$coefs) %*% model$SV
    q[which(q<0)]<-0
    if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
    if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)

    mix_rmse <- rmses[mn]
    mix_r <- corrv[mn]

    newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#do permutations
doPerm <- function(perm, X, Y, absolute, abs_method){
    itor <- 1
    Ylist <- as.list(data.matrix(Y))
    dist <- matrix()

    while(itor <= perm){
        #print(itor)

        #random mixture
        yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

        #standardize mixture
        yr <- (yr - mean(yr)) / stats::sd(yr)

        #run CIBERSORT core algorithm
        result <- CoreAlg(X, yr, absolute, abs_method)

        mix_r <- result$mix_r

        #store correlation
        if(itor == 1) {dist <- mix_r}
        else {dist <- rbind(dist, mix_r)}

        itor <- itor + 1
    }
    newList <- list("dist" = dist)
}
