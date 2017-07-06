
#' Incindence from testing history
#'
#' @param report_hiv_pos A logical vector indicating whether each subject reported a positive hiv status
#' @param biomarker_art A logical vector indicating whether ART antibodies are present. NA if test not done.
#' @param low_viral A logical vector indicating whether viral load is <=1000
#' @param hiv A logical vector indicating hiv status
#' @param last_hiv_test A numeric vector indicating the lower bound of the time since last hiv test bin in months. See details.
#' @param ever_hiv_test A logical vector indicating whether the subject had eer had an hiv test.
#' @param weights Survey weights
#' @param distribution Either "weibull" or "gamma." This controls the family of distribution used to
#' model time since last test.
#' @param ptruth The proportion of the diagnosed hiv positive population that would report being hiv positive.
#' If NULL, this is estimated using biomarker_art and low_viral.
#' @param ptreated The proportion of hiv positive individuals with postive biomarker_art or low_viral.
#' If NULL, this is estimated from the data.
#' @details
#' last_hiv_test assumes that last hiv test information is captured in bins,
#' for example 0-6 month ago, 6-12 months ago, 12-24 month ago, and 24 or greater months ago.
#' Subjects in the first bin would have a last_hiv_test value of 0, those in the second would have 6,
#' the third would have 12 and the last would have 24. Those with missing data or who never had an hiv test
#' should be assigned NA.
#' @return A list with elements:
#' 'inc': the estimated incidence.
#' 'pundiag': The estimated proportion of positive cases that are undiagnosed.
#' 'psay_undiag': The proportion of positive cases that report being undiagnosed.
#' 'pmiss_class': the proportion whose diagnosis status is incorrectly reported by the individual .
#' 'phiv': The proportion with a positive diagnosis.
#' 'ptester': The proportion who have ever been tested.
#' 'mean_time_since_last_test': the mean tie since last test in years.
#' 'tid': mean time between infection and diagnosis.
#' 'ptruth' the proportion of positive individuals that correctly report their status.
#' 'ptreated': the proportion of positive indivudals who are identified as treated by viral load or biomarker.
#' @export
testing_incidence <- function(report_hiv_pos, biomarker_art, low_viral, hiv,
                              last_hiv_test, ever_hiv_test,
                              weights=rep(1, length(report_hiv_pos)) / length(report_hiv_pos),
                              distribution="weibull", ptruth=NULL, ptreated=NULL){

  treated <- low_viral | biomarker_art
  treated[is.na(treated)] <- FALSE


  tab <- wtd.table(!report_hiv_pos, hiv, weights=weights)

  psay_undiag <- tab[2,2] / sum(tab[,2])

  phiv <- as.vector(prop.table(wtd.table(hiv, weights=weights))[2])
  ln_sub <- !hiv & !is.na(hiv) & !is.na(last_hiv_test) & !is.na(weights)
  last_neg <- last_hiv_test[ln_sub]
  ln_wts <- weights[ln_sub]
  rate <- 1/14.03188
  tln <- wtd.table(last_neg, weights = ln_wts)
  t <- as.numeric(names(tln))

  #The data likelihood for time since last negative test
  lik2 <- function(scale,k){
    if(distribution == "weibull")
      p <- function(x) suppressWarnings(pweibull(x, scale=scale,shape=k))
    else
      p <- function(x) suppressWarnings(pgamma(x, scale=scale,shape=k))
    result <- 0
    n <- length(t)
    for(i in 1:(n - 1)){
      result <- result - tln[i] * log(p(t[i+1]) - p(t[i]))
    }
    result <- result - tln[n] * log(1 - p(t[n]))
    result
  }
  opt <- optim(function(x)lik2(x[1],x[2]),par = c(1/rate,1))
  opt

  prop_gt <- rev(cumsum(rev(tln))) / sum(tln)
  #plot(t,prop_gt,ylim=c(0,1))
  #lines(t,pweibull(t,scale=opt$par[1],shape=opt$par[2],lower.tail = FALSE),col="red")

  #mean time between infection and test
  if(distribution == "weibull"){
    m2 <- opt$par[1] * gamma(1 + 1/ opt$par[2])
  }else{
    m2 <- opt$par[1] * opt$par[2]
  }

  # For those who test, time between infection and test is m2, for those who don't test,
  # we assume diagnosis at AIDS, and use the natural history distribution
  ptester <- as.vector(prop.table(wtd.table(ever_hiv_test[!hiv],weights=weights[!hiv]))[2])
  tid <- m2 * ptester +  12 * (1 / 0.086) * gamma(1 + 1 / 2.516) * (1 - ptester)

  # Adjust for missreporting of HIV status
  tlie <- wtd.table(treated, !report_hiv_pos, weights = weights)
  if(is.null(ptruth))
    ptruth <- tlie[2,1] / sum(tlie[2,])
  if(is.null(ptreated))
    ptreated <- tlie[2,2] / sum(tlie[,2])
  pmiss_class <- 1 - (ptruth + ptreated * ( 1 - ptruth))

  pundiag <- (psay_undiag - pmiss_class) / (1 - pmiss_class)

  # Calculate incidence
  tid <- tid / 12
  inc <- pundiag * phiv / (tid * (1-phiv))

  result <- list(inc=inc,
                 pundiag=pundiag,
                 psay_undiag=psay_undiag,
                 pmiss_class=pmiss_class,
                 phiv=phiv,
                 ptester=ptester,
                 mean_time_since_last_test=m2 / 12,
                 tid=tid,
                 ptruth=ptruth,
                 ptreated=ptreated)
  class(result) <- "test-inc"
  result
}


#' Perorm a survey bootstrap
#' @param design an object of class svydesign
#' @param fun a function taking a dataframe as a parameter
#' @param ... additional parameters passed to as.srvrepdesign
#' @return
#' a list containing
#' 'value': the value of the function applied to dat.
#' 'var': the boostrap variance.
#' 'replicates': The bootstrap values of fun.
#' 'nrep': the number of replicates.
#' @examples
#' library(survey)
#' data(api,package="survey")
#' ## one-stage cluster sample
#' dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
#' ## convert to JK1 jackknife
#' rclus1<-as.svrepdesign(dclus1)
#' ## convert to bootstrap
#' set.seed(1)
#' bclus1<-as.svrepdesign(dclus1,type="bootstrap", replicates=100)
#' attr(svymean(~api00, bclus1),"var")
#' set.seed(1)
#' survey_bootstrap(dclus1, fun=function(dat, wts) {sum(wts*dat$api00) / sum(wts)},
#'   type="bootstrap", replicates=100)$var
#' @export
survey_bootstrap <- function(design, fun, ...){
  dat <- design$variables
  weights <- weights(design)

  coef <- fun(dat, weights)

  des1 <- as.svrepdesign(design, compress=FALSE, ...)
  scale <- des1$scale
  rscales <- des1$rscales
  mse <- des1$mse

  rep_weights <- weights(des1)
  nrep <- ncol(rep_weights)
  estimates <- rep(NA, nrep)
  for(i in 1:nrep){
    cat(".")
    estimates[i] <- fun(dat, weights * rep_weights[,i])
  }
  cat("\n")
  list(value=coef,
       var = svrVar(estimates, scale, rscales, mse=mse,coef=coef),
       replicates=estimates,
       nrep=nrep)
}


#' Performs a bootstrap assuming a simple random sample
#' @param dat a data.frame
#' @param fun a function taking a dataframe as a parameter
#' @param replicates the number of replicates
#' @return
#' a list containing
#' 'value': the value of the function applied to dat.
#' 'var': the boostrap variance.
#' 'replicates': The bootstrap values of fun.
#' @examples
#' df <- data.frame(a=rnorm(100))
#' srs_bootstrap(df, colMeans,replicates=100)$var
#' @export
srs_bootstrap <- function(dat, fun, replicates=1000){
  value <- fun(dat)
  n <- nrow(dat)
  estimates <- rep(NA, replicates)
  for(i in 1:replicates){
    samp <- sample.int(n, n, replace=TRUE)
    boot <- dat[samp, , drop=FALSE]
    estimates[i] <- fun(boot)
  }
  list(value=value,
       var=var(estimates),
       replicates=estimates)
}
