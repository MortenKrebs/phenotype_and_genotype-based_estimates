setwd("~/Desktop/AJHG Revisions/plots251103/")

library(PAFGRS)
r2_quant = function(N,M,h2_l){
  R2gg <- 1/(1+M/(N*h2_l))
  R2gg*h2_l}

r2l_to_cor <- function( r2l_pgs, r2l_fgrs, h2_l ) {
  
  # Expected correlation between instruments given their squared correlations with liability
  
  r2l_pgs <- r2l_pgs				# Liability variance explained by PGS
  r2l_fgrs <- r2l_fgrs			# Liability variance explained by FGRS
  h2_l <- h2_l					# Heritability of liability
  
  sqrt( r2l_pgs*r2l_fgrs ) / h2_l
  
}


## Libiraries

library( data.table )
#source( '~/Downloads/Code/PGS_FGRS_Theory.R' )

# Functions

gg_color_hue <- function( n, l=65, c=100 ) {
	hues = seq(15, 375, length = n + 1)
	hcl( h=hues, l=l, c=c )[1:n]
}


h2s <- c( 0.7, 0.7, 0.4, 0.5,0.5,.8 )
prevs <- c( 0.01, 0.05, 0.125, 0.5,0,0 )

for ( m in c(10000, 100000, 60000,1000000 ) ) {
for ( h2.p in c( 1,4/3,2,3,4 ) ) {
  M <- m
  
pdf( paste0( "Correlations_h2p_", round(h2.p,1), "_M_", M, ".pdf" ), height=8, width=6 )
  par( mar=c(3.5,4.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25, mfrow=c(3,2) )
  
  for ( trait in 1:5 ) {

   # pdf( paste0( "Correlations_h2p_", round(h2.p,1), "_M_", M, ".pdf" ), height=5, width=5 )
    
# Figure 5A
# Expected Correlation as a functon of GWAS sample
# GWAS lines: population sample, case-control sample
# FGRS lines: 1, 3, 5 sibling equivalents

t <- trait
h2 <- h2s[ trait ]
h2_snp <- h2/h2.p 
prev <- prevs[ trait ]


N <- 10^seq( 3, 7, by=0.001 )
N_case_pop <- N*prev
N_control_pop <- N*(1-prev)
N_case_cc <- N_control_cc <- N*0.5

if(!trait%in%5:6){
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





er_pop_1 <- r2l_to_cor( pgs.pop, fgrs.1[[1]]^2*h2, h2 )
er_cc_1 <- r2l_to_cor( pgs.cc, fgrs.1[[1]]^2*h2, h2 )
er_pop_3 <- r2l_to_cor( pgs.pop, fgrs.3[[1]]^2*h2, h2 )
er_cc_3 <- r2l_to_cor( pgs.cc, fgrs.3[[1]]^2*h2, h2 )
er_pop_5 <- r2l_to_cor( pgs.pop, fgrs.5[[1]]^2*h2, h2 )
er_cc_5 <- r2l_to_cor( pgs.cc, fgrs.5[[1]]^2*h2, h2 )



	plot( log10( N ), er_pop_1, type='l', 
			lty=3, 
			col=c(gg_color_hue( 4, l=85 ),"black","black")[t],
			xlab=expression(log[10] ~ (Sample ~ Size) ),
			ylab="",
			yaxt='n',
			ylim=c(0,0.5), 
			lwd=1.5
			)
	mtext(expression( 'E(' ~ r[PGS ~ ',' ~ PA-FGRS] ~ ')' ), side = 2, line = 2.5)
	mtext(LETTERS[trait],side=2, line=3, at=par('usr')[4], las=2,cex=1.7)
	
	axis(2, at=seq( 0,0.5,by=0.05 ), las=2 )

	lines( log10( N ), er_pop_3, lty="34", col=c(gg_color_hue( 4, l=85 ),"black","black")[t], lwd=1.5 )
	lines( log10( N ), er_pop_5, lty="69", col=c(gg_color_hue( 4, l=85 ),"black","black")[t], lwd=1.5 )

	if( prev != 0.5 ) {

		lines( log10( N ), er_cc_1, lty=3, col=c(gg_color_hue( 4, l=45 ),"black","black")[t], lwd=1.5 )
		lines( log10( N ), er_cc_3, lty="34", col=c(gg_color_hue( 4, l=45 ),"black","black")[t], lwd=1.5 )
		lines( log10( N ), er_cc_5, lty="69", col=c(gg_color_hue( 4, l=45 ),"black","black")[t], lwd=1.5 )

	}
	if(!trait%in%5:6){
	text( 6.8, er_cc_5[3801], 'c.c.', pos=3, offset=0.5, col=c(gg_color_hue( 4, l=45 ),"black","black")[t] )
	text( 6.8, er_pop_1[3801], 'pop.', pos=1, offset=0.5, col=c(gg_color_hue( 4, l=85 ),"black","black")[t] )
	}
	l1 <- as.expression( substitute( paste( h^2, " = ", h2t, sep='' ), list( h2t=h2 ) ) )
	l2 <- as.expression( substitute( paste( h[SNP]^2, " = ", h2t, sep='' ), list( h2t=round( h2_snp,2 ) ) ) )
	if(!trait%in%5:6){l3 <- paste( "prev = ", prev, sep='' )}else 
	  l3 <- NULL
	if(m<1e6) l4 <- paste( "m = ", m/1e3,",","000", sep='' ) else 
	  l4 <- paste( "m = ", m/1e6,",000,000", sep='' )
	legend( x=2.5,y=.57, 
			legend=c( "",l1, l2, l3 ,l4 ),
			bty='n',
			pch=NULL,text.col = c(gg_color_hue( 4, l=65 ),"black","black")[t]  )		



  }
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend( 'left',
          legend=c( expression(5*" "*N[sibe]),
                    expression(3*" "*N[sibe]),
                    expression(1*" "*N[sibe])),
          bty='n',
          lty=c( "69", "34", "13" ),
          col=c( 'black', 'black', 'black' ),
          cex=1.5,
          lwd = 1.5)
  
dev.off()
  }
}  



