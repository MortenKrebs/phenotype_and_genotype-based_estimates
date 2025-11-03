library( data.table )
library(PAFGRS)

setwd("~/Desktop/AJHG Revisions/plots251103/")

r2_quant = function(N,M,h2_l){
  R2gg <- 1/(1+M/(N*h2_l))
  R2gg*h2_l}

# Functions

r2l_to_cor <- function( r2l_pgs, r2l_fgrs, h2_l ) {
  
  # Expected correlation between instruments given their squared correlations with liability
  
  r2l_pgs <- r2l_pgs				# Liability variance explained by PGS
  r2l_fgrs <- r2l_fgrs			# Liability variance explained by FGRS
  h2_l <- h2_l					# Heritability of liability
  
  sqrt( r2l_pgs*r2l_fgrs ) / h2_l
  
}

r2l_to_semipartial_cor_l <- function( r2l_1, r2l_2, h2_l ) {
  
  # Expected semi-partial correlation between instrument 1 and liability after adjustment by instrument 2 
  
  r2l_1 <- r2l_1					# Liability variance explained by instrument 1 (to be adjusted)
  r2l_2 <- r2l_2					# Liability variance explained by instrument 2 (to be adjusted by)
  h2_l <- h2_l					# Heritability of liability
  
  acc_1_l <- sqrt( r2l_1 )
  acc_2_l <- sqrt( r2l_2 )
  cor_1_2 <- r2l_to_cor( r2l_1, r2l_2, h2_l )
  
  ( acc_1_l - cor_1_2*acc_2_l ) / ( 1 - (cor_1_2)^2 )
  
}

r2l_marginal_to_joint <- function( r2l_1, r2l_2, h2_l ) {
  
  # Expected joint variance explained, given marginal
  
  r2l_1 <- r2l_1					# Liability variance explained by instrument 1
  r2l_2 <- r2l_2					# Liability variance explained by instrument 2
  h2_l <- h2_l					# Heritability of liability
  
  spcor <- r2l_to_semipartial_cor_l( r2l_2, r2l_1, h2_l )
  
  r2l_1 + spcor^2
  
}

FGRS_by_nSib <- function( n_sib, prev, h2_l ) {
  
  n_rel <- n_sib						# number of relatives
  r_ir <- 0.5							# relatedness of relative to index
  r_rr <- 0.5							# relatedness of relatives to each other
  k <- prev							# population lifetime prevalence
  h2 <- h2_l							# liability scale h2
  t <- qnorm( 1-k )					# case threshold
  phi_t <- dnorm( qnorm( 1-k ) )		# average case liability
  
  fgrs_rel <- (n_rel*r_ir^2) / ( (k*(1-k))/(h2*phi_t^2) + (n_rel-1)*r_rr*(1 + (r_rr*h2*t^2)/2) ) 
  fgrs_acc <- sqrt( fgrs_rel )
  fgrs_r2l <- fgrs_rel*h2
  
  c( acc=fgrs_acc, rel=fgrs_rel, r2l=fgrs_r2l )
  
}

gg_color_hue <- function( n, l=65, c=100 ) {
	hues = seq(15, 375, length = n + 1)
	hcl( h=hues, l=l, c=c )[1:n]
}


h2s <- c( 0.7, 0.7, 0.4, 0.5 ,0.5 )
prevs <- c( 0.01, 0.05, 0.125, 0.5 ,0 )

for ( m in c(10000, 100000, 60000, 1000000 ) ) {
for ( h2.p in c( 1,4/3,2,3,4 ) ) {
  
  pdf( paste0( "Jointr2l_", round(h2.p,1), "_M_", m, ".pdf" ), height=6, width=5 )
  par(mfrow=c(3,2))
  for ( trait in 1:5 ) {
    

t <- trait
h2 <- h2s[ trait ]
h2_snp <- h2/h2.p 
prev <- prevs[ trait ]
M <- m

N <- 10^seq( 3, 7, by=0.001 )
N_case_pop <- N*prev
N_control_pop <- N*(1-prev)
N_case_cc <- N_control_cc <- N*0.5

if(!trait==5){
  fgrs.1 <- fgrs_accuracy( nrel=1, prev=prev , h2=h2, 
                           rel_rr=0.5, rel_ir=0.5, method="pa", estimate="theory" )
  fgrs.3 <- fgrs_accuracy( nrel=3, prev=prev , h2=h2, 
                           rel_rr=0.5, rel_ir=0.5, method="pa", estimate="theory" )
  fgrs.5 <- fgrs_accuracy( nrel=5, prev=prev , h2=h2, 
                           rel_rr=0.5, rel_ir=0.5, method="pa", estimate="theory" )
  pgs.pop <- PAFGRS:::r2_cc_sample( N_case_pop, N_control_pop, M, prev, h2_snp)
  pgs.cc <- PAFGRS:::r2_cc_sample( N_case_cc, N_control_cc, M, prev, h2_snp)
}else {
  fgrs.1 <- fgrs_accuracy( nrel=1, prev=prev , h2=h2, 
                           rel_rr=0.5, rel_ir=0.5, method="quantitative", estimate="theory" )
  fgrs.3 <- fgrs_accuracy( nrel=3, prev=prev , h2=h2, 
                           rel_rr=0.5, rel_ir=0.5, method="quantitative", estimate="theory" )
  fgrs.5 <- fgrs_accuracy( nrel=5, prev=prev , h2=h2, 
                           rel_rr=0.5, rel_ir=0.5, method="quantitative", estimate="theory" )
  
  
  pgs.pop <- r2_quant(N, M, h2_snp)
  pgs.cc <- r2_quant(N, M, h2_snp)
  
}

er2l_pop_1 <- r2l_marginal_to_joint( pgs.pop, fgrs.1[[1]]^2*h2, h2 )
er2l_pop_1.1 <- r2l_marginal_to_joint( fgrs.1[[1]]^2*h2, pgs.pop, h2 )
er2l_cc_1 <- r2l_marginal_to_joint( pgs.cc, fgrs.1[[1]]^2*h2, h2 )
er2l_cc_1.1 <- r2l_marginal_to_joint( fgrs.1[[1]]^2*h2, pgs.cc, h2 )
er2l_pop_3 <- r2l_marginal_to_joint( pgs.pop, fgrs.3[[1]]^2*h2, h2 )
er2l_pop_3.1 <- r2l_marginal_to_joint( fgrs.3[[1]]^2*h2, pgs.pop, h2 )
er2l_cc_3 <- r2l_marginal_to_joint( pgs.cc, fgrs.3[[1]]^2*h2, h2 )
er2l_cc_3.1 <- r2l_marginal_to_joint( fgrs.3[[1]]^2*h2, pgs.cc, h2 )
er2l_pop_5 <- r2l_marginal_to_joint( pgs.pop, fgrs.5[[1]]^2*h2, h2 )
er2l_pop_5.1 <- r2l_marginal_to_joint( fgrs.5[[1]]^2*h2, pgs.pop, h2 )
er2l_cc_5 <- r2l_marginal_to_joint( pgs.cc, fgrs.5[[1]]^2*h2, h2 )
er2l_cc_5.1 <- r2l_marginal_to_joint( fgrs.5[[1]]^2*h2, pgs.cc, h2 )

	par( mar=c(3.5,4.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )
	if(!trait==5)
	  ylab=expression( 'E(' ~ R[l]^2 ~ ')' ) else
	    ylab=expression( 'E(' ~ R^2 ~ ')' )
	plot( log10( N ), pgs.pop, type='l', 
			lty=1, 
			col=c(gg_color_hue( 4, l=85 ),"black")[t],
			xlab=expression(log[10] ~ (Sample ~ Size) ),
			ylab="",
			yaxt='n',
			ylim=c(0,0.4), 
			lwd=1.5
			)
	mtext(ylab, side = 2, line = 2.5)
	mtext(LETTERS[trait],side=2, line=3, at=par('usr')[4], las=2,cex=1.7)
	
	axis(2, at=seq( 0,0.5,by=0.05 ), las=2 )

	lines( log10( N ), er2l_pop_1, lty="13", col=c(gg_color_hue( 4, l=85 ),"black")[t], lwd=1.5 )
	lines( log10( N ), er2l_pop_3, lty="34", col=c(gg_color_hue( 4, l=85 ),"black")[t], lwd=1.5 )
	lines( log10( N ), er2l_pop_5, lty="69", col=c(gg_color_hue( 4, l=85 ),"black")[t], lwd=1.5 )

	if( prev != 0.5 ) {

		lines( log10( N ), pgs.cc, lty=1, col=c(gg_color_hue( 4, l=45 ),"black")[t], lwd=1.5 )
		lines( log10( N ), er2l_cc_1, lty="13", col=c(gg_color_hue( 4, l=45 ),"black")[t], lwd=1.5 )
		lines( log10( N ), er2l_cc_3, lty="34", col=c(gg_color_hue( 4, l=45 ),"black")[t], lwd=1.5 )
		lines( log10( N ), er2l_cc_5, lty="69", col=c(gg_color_hue( 4, l=45 ),"black")[t], lwd=1.5 )

	}
	if(!trait==5){
	text( 6.8, er2l_cc_5[3801], 'c.c.', pos=3, offset=0.5, col=c(gg_color_hue( 4, l=45 ),"black")[t] )
	text( 6.8, pgs.pop[3801], 'pop.', pos=1, offset=0.5, col=c(gg_color_hue( 4, l=85 ),"black")[t] )
	}
	
	l1 <- as.expression( substitute( paste( h^2, " = ", h2t, sep='' ), list( h2t=h2 ) ) )
	l2 <- as.expression( substitute( paste( h[SNP]^2, " = ", h2t, sep='' ), list( h2t=round( h2_snp,2 ) ) ) )
	if(!trait==5){l3 <- paste( "prev = ", prev, sep='' )}else l3 <- NULL
	if(m<1e6) l4 <- paste( "m = ", m/1e3,",","000", sep='' ) else l4 <- paste( "m = ", m/1e6,",000,000", sep='' )
	legend(
	        legend=c( "",l1, l2, l3 ,l4 ),
	        bty='n',
	        pch=NULL,x=2.5,y=.47,
	        text.col=c(gg_color_hue( 4, l=45 ),"black")[t])		
	

  }					
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend( 'left',
          legend=c( expression(5*" "*N[sibe] (+PGS)), 
                    expression(3*" "*N[sibe] (+PGS)),
                    expression(1*" "*N[sibe] (+PGS)), 'PGS' ),
          bty='n',
          lty=c( "69", "34", "13","solid"),
          col=c( 'black', 'black', 'black','black' ),
          cex=1.5,
          lwd = 1.5)
dev.off()


}
}



