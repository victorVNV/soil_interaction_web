#' Evaluate continuos heavy tailed distributions in a vector
#'
#' @param br vector with distribution   
#' @param returnEWS if TRUE returns the early warnings, FALSE returns the patch distribution   
#' @param returnOBJ if TRUE returns the early warnings, FALSE returns the patch distribution   
#' @param est_xmin  if TRUE 
#' @return a data frame with results
#' @export
#'
#' @examples
evaluate_distr <- function(distr,returnOBJ=FALSE,est_xmin=FALSE){
  
  source("R/powerlaw/discpowerexp.R")
  source("R/powerlaw/discexp.R")
  source("R/powerlaw/zeta.R")
  source("R/powerlaw/powerexp.R")
  source("R/powerlaw/exp.R")
  source("R/powerlaw/pareto.R")
  
  require(ggplot2)
  
  ## Convert to TRUE/FALSE matrix
  #

  require(poweRlaw)
  m_pl = conpl$new(distr)
  if(est_xmin) {
    est = estimate_xmin(m_pl)
    m_pl$setXmin(est)
  } else {
    est = estimate_pars(m_pl)
    m_pl$setPars(est)
    est$xmin <- m_pl$xmin
  }
  dd <- plot(m_pl) 
  ll <- lines(m_pl, col = 2)
  m_ln <- conlnorm$new(distr)
  m_ln$setXmin(est$xmin)
  pars <- estimate_pars(m_ln)
  m_ln$setPars(pars)
  pp <- lines(m_ln)
  
  m_exp <- conexp$new(distr)
  m_exp$setXmin(est$xmin)
  pars <- estimate_pars(m_exp)
  m_exp$setPars(pars)
  ee <- lines(m_exp)
  
  m_pexp <- powerexp.fit(distr,est$xmin)
  pe <- data.frame(x=dd$x,y=ppowerexp(dd$x,est$xmin,m_pexp$exponent,m_pexp$rate,lower.tail=F))
  if( est_xmin ){
    pe$y <- pe$y* dd$y[est$xmin]
  }
  gp <- ggplot( dd, aes(x,y) ) + geom_point() + theme_bw() + scale_x_log10() + scale_y_log10() + 
    geom_line(data=ll,aes(x,y, color="Pl.")) +  
    geom_line(data=ee,aes(x,y, color="Exp.")) + ylab("Frequency (P>=x)") + xlab( "Patch size") +
    geom_line(data=pe,aes(x,y, color="P.Exp.")) + 
    geom_line(data=pp,aes(x,y, color="Ln.")) + coord_cartesian(ylim=range(dd$y)) + scale_color_viridis_d(name="")
  
  if(returnOBJ)
    patch_df <- gp
  else {

    df1 <-tibble(type="pl", expo=m_pl$pars[1],rate=m_pl$pars[2],xmin=m_pl$xmin,AICc = calc_AICc(m_pl), range=max(distr)-xmin )
    df2 <-tibble(type="exp", expo=m_exp$pars[1],rate=m_exp$pars[2],xmin=m_exp$xmin,AICc = calc_AICc(m_exp)) 
    df3 <- tibble(type="ln", expo=m_ln$pars[1],rate=m_ln$pars[2],xmin=m_ln$xmin,AICc = calc_AICc(m_ln))
    df4 <- tibble(type="pexp", expo=m_pexp$exponent,rate=m_pexp$rate,xmin=est$xmin,AICc = calc_AICc(m_pexp))
    patch_df <- bind_rows(df1,df2,df3,df4)
  }
  
  return(patch_df)
}

# Calculate the AICc from a poweRlaw distribution object
#
calc_AICc <- function( d_obj){
  if( class(d_obj)=="list") {
    LL <- d_obj$loglike
    n <- d_obj$samples.over.threshold
    k <- 2
  } else {
    LL <- dist_ll(d_obj)
    n <- get_ntail(d_obj)
    k <- d_obj$no_pars
  }
  AICc <- (2*k-2*LL)+2*k*(k+1)/(n-k-1)
}