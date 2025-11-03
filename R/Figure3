library( data.table )

r2_quant = function(N,M,h2_l){
   R2gg <- 1/(1+M/(N*h2_l))
   R2gg*h2_l}
 


gg_color_hue <- function( n, l=65, c=100 ) {
	hues = seq(15, 375, length = n + 1)
	hcl( h=hues, l=l, c=c )[1:n]
}

h2s <- c( 0.7, 0.7, 0.4, 0.5 , .5 )
prevs <- c( 0.01, 0.05, 0.125, 0.5 ,0)

setwd("~/Desktop/AJHG Revisions/plots251103/")
for ( m in c(10000, 100000, 60000, 1000000 ) ) {
#for ( m in  ) {

for ( h2.p in c( 1,2,4 ) ) {
		
pdf( paste0( "FGRSvPGS_h2p_", 
						round(h2.p,1), "_M_", m, ".pdf" ), height=8, width=6 )

par( mar=c(3.5,4.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25, mfrow=c(3,2) )
	
for ( trait in 1:5 ) {
	
t <- trait
h2 <- h2s[ trait ]
h2_snp <- h2/h2.p 
prev <- prevs[ trait ]
M <- m
ymax <- max( 0.05 + h2s/h2.p ) 


if(trait<5){
fgrs.1s.M.new <- PAFGRS::fgrs_accuracy( nrel=1, prev=prev , h2=h2, 
                          rel_rr=0.5, rel_ir=0.5, method="pa", estimate="theory" )
fgrs.3s.M.new <- PAFGRS::fgrs_accuracy( nrel=3, prev=prev , h2=h2, 
                          rel_rr=0.5, rel_ir=0.5, method="pa", estimate="theory" )
fgrs.5s.M.new <- PAFGRS::fgrs_accuracy( nrel=5, prev=prev , h2=h2, 
                          rel_rr=0.5, rel_ir=0.5, method="pa", estimate="theory" )

N <- 10^seq( 3, 7, by=0.001 )
N_case_pop <- N*prev
N_control_pop <- N*(1-prev)
N_case_cc <- N_control_cc <- N*0.5

pgs.pop.M <- PAFGRS:::r2_cc_sample( N_case_pop, N_control_pop, M, prev, h2_snp)
pgs.cc.M <- PAFGRS:::r2_cc_sample( N_case_cc, N_control_cc, M, prev, h2_snp)
}else{
  pgs.pop.M <- pgs.cc.M <- r2_quant( N,M,h2_snp)
  fgrs.1s.M.new <- PAFGRS::fgrs_accuracy( nrel=1, prev=prev , h2=h2, 
                                    rel_rr=0.5, rel_ir=0.5, method="quantitative", estimate="theory" )
  fgrs.3s.M.new <- PAFGRS::fgrs_accuracy( nrel=3, prev=prev , h2=h2, 
                                  rel_rr=0.5, rel_ir=0.5, method="quantitative", estimate="theory" )
  fgrs.5s.M.new <- PAFGRS::fgrs_accuracy( nrel=5, prev=prev , h2=h2, 
                                  rel_rr=0.5, rel_ir=0.5, method="quantitative", estimate="theory" )
  
}

	plot( log10( N ), pgs.pop.M, type='l', 
			lty=1, 
			col=gg_color_hue( 4, l=85 )[t],
			xlab=expression(log[10] ~ (Sample ~ Size) ),
			yaxt='n',
			ylab="",
			ylim=c(0,ymax), 
			lwd=1.5
			)
	mtext(expression( 'E(' ~ R[Liab.]^2 ~ ')' ), side = 2, line = 2.5)
	mtext(LETTERS[trait],side=2, line=3, at=par('usr')[4], las=2,cex=1.7)
	
	axis(2, at=seq( 0,ymax,by=0.05 ), las=2 )
	lines( log10( N ), pgs.cc.M, lty=1, col=c(gg_color_hue(4,l=45),"black")[t] , lwd=1.5 )

	segments( 0, h2*fgrs.1s.M.new[ 1 ]^2-0.00075, 6.3, h2*fgrs.1s.M.new[ 1 ]^2-0.00075, lty=3, col=c(gg_color_hue(4),"black")[t] , lwd=1.5 )
	segments( 0, h2*fgrs.3s.M.new[ 1 ]^2, 6.3, h2*fgrs.3s.M.new[ 1 ]^2, lty=2, col=c(gg_color_hue(4),"black")[t] , lwd=1.5 )
	segments( 0, h2*fgrs.5s.M.new[ 1 ]^2, 6.3, h2*fgrs.5s.M.new[ 1 ]^2, lty=5, col=c(gg_color_hue(4),"black")[t] , lwd=1.5 )	
	
	text( 6.65, h2*fgrs.1s.M.new[ 1 ]^2, expression(1*" "*N[sibe]), pos=NULL, offset=0, col=c(gg_color_hue(4),"black")[t] )
	text( 6.65, h2*fgrs.3s.M.new[ 1 ]^2, expression(3*" "*N[sibe]), pos=NULL, offset=0, col=c(gg_color_hue(4),"black")[t]  )
	text( 6.65, h2*fgrs.5s.M.new[ 1 ]^2, expression(5*" "*N[sibe]), pos=NULL, offset=0, col=c(gg_color_hue(4),"black")[t]  )

	if(trait<5) {
	text( 6.8, pgs.cc.M[3801], 'c.c.', pos=3, offset=0.5, col=gg_color_hue( 4, l=45 )[t] )
	text( 6.8, pgs.pop.M[3801], 'pop.', pos=1, offset=0.5, col=gg_color_hue( 4, l=85 )[t] )
	}
	l1 <- as.expression( substitute( paste( h^2, " = ", h2t, sep='' ), list( h2t=h2 ) ) )
	l2 <- as.expression( substitute( paste( h[SNP]^2, " = ", h2t, sep='' ), list( h2t=round( h2_snp,2 ) ) ) )
	if(trait<5) l3 <- paste( "prev = ", prev, sep='' ) else(l3 <-NULL)
	l4 <- paste( "m = ", format(M, scientific=F), sep='' )
		
	legend( x=2.6,y=.42 ,
			legend=c( l1, l2, l3, l4 ),
			bty='n',
			text.col = c(gg_color_hue(4),"black")[t],
			pch=NULL )

	
	print(trait)
	print(pgs.pop.M[3801])
	
}

dev.off()

}
}



