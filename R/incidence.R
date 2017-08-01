
#' Incindence from testing history
#'
#' @param report_hiv_pos A logical vector indicating whether each subject reported a positive hiv status
#' @param biomarker_art A logical vector indicating whether ART antibodies are present. NA if test not done.
#' @param low_viral A logical vector indicating whether viral load is <=1000
#' @param hiv A logical vector indicating hiv status
#' @param last_hiv_lower A numeric vector indicating the lower bound of the months since last hiv test bin in months.
#' @param last_hiv_upper A numeric vector indicating the upper bound of the months since last hiv test bin in months.
#' @param ever_hiv_test A logical vector indicating whether the subject had eer had an hiv test.
#' @param weights Survey weights
#' @param distribution Either "weibull" or "gamma." This controls the family of distribution used to
#' model time since last test.
#' @param ptruth The proportion of the diagnosed hiv positive population that would report being hiv positive.
#' If NULL, this is estimated using biomarker_art and low_viral.
#' @param ptreated The proportion of hiv positive individuals with postive biomarker_art or low_viral.
#' If NULL, this is estimated from the data.
#' @param non_tester_tid mean nummber of months between infection and diagnosis for non_testers.
#' @details
#' last_hiv_upper should always be greater than last_hiv_lower, and may be infinite (e.g. test was  > 24 months ago). Those with missing data or who never had an hiv test
#' should be assigned NA for both their lower and upper bound.
#' @return A data.frame of class test_inc with elements:
#' 'inc': the estimated incidence.
#' 'pundiag': The estimated proportion of positive cases that are undiagnosed.
#' 'psay_undiag': The proportion of positive cases that report being undiagnosed.
#' 'pmiss_class': the proportion whose diagnosis status is incorrectly reported by the individual .
#' 'phiv': The proportion with a positive diagnosis.
#' 'ptester': The proportion who have ever been tested.
#' 'mean_time_since_last_test': the mean tie since last test in years.
#' 'tid': mean time between infection and diagnosis in years.
#' 'ptruth' the proportion of positive individuals that correctly report their status.
#' 'ptreated': the proportion of positive indivudals who are identified as treated by viral load or biomarker.
#' @export
testing_incidence <- function(report_hiv_pos, biomarker_art, low_viral, hiv,
                              last_hiv_lower, last_hiv_upper, ever_hiv_test,
                              weights=rep(1, length(report_hiv_pos)) / length(report_hiv_pos),
                              distribution="weibull", ptruth=NULL, ptreated=NULL,
                              non_tester_tid= 123.824){

  treated <- low_viral | biomarker_art
  treated[is.na(treated)] <- FALSE


  tab <- wtd.table(!report_hiv_pos, hiv, weights=weights)

  psay_undiag <- tab[2,2] / sum(tab[,2])

  phiv <- as.vector(prop.table(wtd.table(hiv, weights=weights))[2])
  ln_sub <- !hiv & !is.na(hiv) & !is.na(last_hiv_upper) & !is.na(last_hiv_lower) & !is.na(weights)

  #get counts for testing windows
  test_window <- paste(last_hiv_lower[ln_sub],last_hiv_upper[ln_sub],sep="_")
  ln_count <- wtd.table(test_window, weights = weights[ln_sub])
  ln_windows <- strsplit(names(ln_count),"_")
  ln_lower <- as.numeric(sapply(ln_windows, function(x) x[1]))
  ln_upper <- as.numeric(sapply(ln_windows, function(x) x[2]))
  ln_count <- as.vector(ln_count)

  rate <- 1/14.03188

  #The data likelihood for time since last negative test
  lik2 <- function(scale,k){
    if(distribution == "weibull"){
      p <- function(x) suppressWarnings(pweibull(x, scale=scale,shape=k))
      log_d <- function(x) suppressWarnings(dweibull(x, scale=scale,shape=k, log = TRUE))
    }else{
      p <- function(x) suppressWarnings(pgamma(x, scale=scale,shape=k))
      log_d <- function(x) suppressWarnings(dgamma(x, scale=scale,shape=k, log = TRUE))
    }
    result <- 0

    n <- length(ln_count)
    exact_time <- ln_upper == ln_lower
    for(i in 1:n){
      if(!exact_time[i]){
        result <- result - ln_count[i] * log(p(ln_upper[i]) - p(ln_lower[i]))
      }
    }
    result <- result - sum(ln_count[exact_time] * log_d(ln_upper[exact_time]) )
    result
  }
  opt <- optim(function(x)lik2(x[1],x[2]),par = c(1/rate,1))

  #mean time between infection and test
  if(distribution == "weibull"){
    m2 <- opt$par[1] * gamma(1 + 1/ opt$par[2])
  }else{
    m2 <- opt$par[1] * opt$par[2]
  }

  # For those who test, time between infection and test is m2, for those who don't test,
  # we assume diagnosis at AIDS, and use the natural history distribution
  ptester <- as.vector(prop.table(wtd.table(ever_hiv_test[!hiv],weights=weights[!hiv]))[2])
  tid <- m2 * ptester +  non_tester_tid * (1 - ptester)

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

  result <- data.frame(list(incidence=inc,
                 pundiag=pundiag,
                 psay_undiag=psay_undiag,
                 pmiss_class=pmiss_class,
                 phiv=phiv,
                 ptester=ptester,
                 mean_time_since_last_test=m2 / 12,
                 tid=tid,
                 ptruth=ptruth,
                 ptreated=ptreated))
  row.names(result) <- "estimate"
  class(result) <- c("test_inc","data.frame")
  result
}


#' Perorm a survey bootstrap
#' @param design an object of class svydesign
#' @param fun a function taking a dataframe as a parameter
#' @param show_progress print progress
#' @param ... additional parameters passed to as.srvrepdesign
#' @return
#' a list with class boot_est containing
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
survey_bootstrap <- function(design, fun, show_progress=TRUE, ...){
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
    if(show_progress)
      cat(".")
    estimates[i] <- fun(dat, weights * rep_weights[,i])
  }
  if(show_progress)
    cat("\n")
  bb <- list(value=coef,
       var = svrVar(estimates, scale, rscales, mse=mse,coef=coef),
       replicates=estimates,
       nrep=nrep)
  class(bb) <- c("boot_est","list")
  bb
}

#' Bootstrap variance for incidence
#' @param frame a dataframe conatining the variables needed for id_formula, strata_formula and wieghts_formula
#' @param strata_formula a survey formula specifying the strata structure
#' @param weights_formula a survey formula specifying the weights
#' @param id_formula a survey formula specifying the unit of sampling
#' @param report_hiv_pos A logical vector indicating whether each subject reported a positive hiv status
#' @param biomarker_art A logical vector indicating whether ART antibodies are present. NA if test not done.
#' @param low_viral A logical vector indicating whether viral load is <=1000
#' @param hiv A logical vector indicating hiv status
#' @param last_hiv_lower A numeric vector indicating the lower bound of the months since last hiv test bin in months.
#' @param last_hiv_upper A numeric vector indicating the upper bound of the months since last hiv test bin in months.
#' @param ever_hiv_test A logical vector indicating whether the subject had eer had an hiv test.
#' @param replicates The number of bootstrap replicated
#' @param type The type of survey bootstrap
#' @param show_progress print progress
#' @param ... additional parameters for testing_incidence
#' @export
bootstrap_incidence <- function(report_hiv_pos, biomarker_art, low_viral, hiv,
                                last_hiv_lower, last_hiv_upper, ever_hiv_test, replicates=5000,
                                frame=NULL, id_formula=NULL, strata_formula=NULL, weights_formula=NULL, type="mrbbootstrap", show_progress=TRUE, ...){
  if(is.null(frame) && is.null(id_formula) && is.null(strata_formula) && is.null(weights_formula)){
    dat <- data.frame(report_hiv_pos,
      biomarker_art,
      low_viral,
      hiv,
      last_hiv_lower,
      last_hiv_upper,
      ever_hiv_test,
      weights = rep(1,length(report_hiv_pos))
      )
    bb <- srs_bootstrap(dat,
                           fun=function(dat) testing_incidence(dat$report_hiv_pos, dat$biomarker_art, dat$low_viral, dat$hiv,
                                                                    dat$last_hiv_lower, dat$last_hiv_upper, dat$ever_hiv_test, dat$weights, ...)$incidence,
                           replicates=replicates,
                        show_progress=show_progress)
  }else{
    include <- rowSums(is.na(frame)) == 0
    frame$report_hiv_pos <- report_hiv_pos
    frame$biomarker_art <- biomarker_art
    frame$low_viral <- low_viral
    frame$hiv <- hiv
    frame$last_hiv_lower <- last_hiv_lower
    frame$last_hiv_upper <- last_hiv_upper
    frame$ever_hiv_test <- ever_hiv_test

    frame <- frame[include,]
    des <- svydesign(ids = id_formula, strata=strata_formula, weights =weights_formula, data=frame)
    bb <- survey_bootstrap(des,
                           fun=function(dat, wts) testing_incidence(dat$report_hiv_pos, dat$biomarker_art, dat$low_viral, dat$hiv,
                                                                    dat$last_hiv_lower, dat$last_hiv_upper, dat$ever_hiv_test, wts, ...)$incidence,
                           show_progress=show_progress,
                           type=type,
                           replicates=replicates)
  }
  bb
}

#' Performs a bootstrap assuming a simple random sample
#' @param dat a data.frame
#' @param fun a function taking a dataframe as a parameter
#' @param replicates the number of replicates
#' @param show_progress print progress
#' @return
#' a list with class boot_est containing
#' 'value': the value of the function applied to dat.
#' 'var': the boostrap variance.
#' 'replicates': The bootstrap values of fun.
#' 'nrep': the number of replicates.
#' @examples
#' df <- data.frame(a=rnorm(100))
#' srs_bootstrap(df, colMeans,replicates=100)$var
#' @export
srs_bootstrap <- function(dat, fun, replicates=1000, show_progress=TRUE){
  value <- fun(dat)
  n <- nrow(dat)
  estimates <- rep(NA, replicates)
  for(i in 1:replicates){
    if(show_progress)
      cat(".")
    samp <- sample.int(n, n, replace=TRUE)
    boot <- dat[samp, , drop=FALSE]
    estimates[i] <- fun(boot)
  }
  if(show_progress)
    cat("\n")
  bb <- list(value=value,
       var=var(estimates),
       replicates=estimates,
       nrep=replicates)
  class(bb) <- c("boot_est","list")
  bb
}

#' print bootstrap
#' @param x a boot_est object
#' @param probs A list of quantiles for the bootstrap
#' @param ... additional parameters for print.data.frame
#' @export
print.boot_est <- function(x,probs=c(.025,.975), ...){
  res <- data.frame(value=x$value,
                    se=sqrt(x$var))
  res <- cbind(res,t(quantile(x$replicates, probs=probs)))
  print(res, ...)
}

