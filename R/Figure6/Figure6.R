rm( list=ls() )
set.seed(123)

## Libraries

library( data.table )
 
gg_color_hue <- function( n, l=65, c=100, alpha=1 ) {
	hues = seq(15, 375, length = n + 1)
	hcl( h=hues, l=l, c=c, alpha=alpha )[1:n]
}


## Sham Equations

Sham_r2l_cc_sample = function( N_case, N_control, M=60000, prev, h2_l_train, h2_l_test=h2_l_train, rg_test_train=1 )  {
  
  #prev in training sample
  
  N <- ( N_case + N_control )
  w <- N_case / N
  c <- w*(1-w)*( dnorm( x=qnorm( 1-prev ) ) / ( prev*(1-prev) ) )^2
  
  Nquant <- N*c
  
  R2gg_train <- ( h2_l_train*Nquant ) / ( h2_l_train*Nquant + M )
  
  r2l_test <- R2gg_train * rg_test_train^2 * h2_l_test 
  
  r2l_test
  
} 

Sham_r2l_to_m = function( N_case_train, N_control_train, r2l_test, prev, h2_l_train, h2_l_test=h2_l_train, rg_test_train=1 )  {
  
  N <- ( N_case_train + N_control_train )
  w <- N_case_train / N
  c <- w*(1-w)*( dnorm( x=qnorm( 1-prev ) ) / ( prev*(1-prev) ) )^2
  
  Nquant <- N*c
  
  r2gg_train <- r2l_test / ( h2_l_test*rg_test_train^2 )
  
  M <- Nquant * h2_l_train * ( 1-r2gg_train ) / r2gg_train
  
  M
  
} 

Sham_r2l_to_n = function( r2l_test, M, prev, h2_l_train, h2_l_test=h2_l_train, rg_test_train=1, w_train=prev )  {
  
  c <- w_train*(1-w_train)*( dnorm( x=qnorm( 1-prev ) ) / ( prev*(1-prev) ) )^2		
  
  r2gg_train <- r2l_test / ( h2_l_test*rg_test_train^2 )
  
  Nq_train <- M * ( 1/h2_l_train ) * ( r2gg_train/( 1-r2gg_train ) )
  
  N_train = Nq_train / c
  
  N_train
  
} 


acc_to_se <- function( acc, prev=0.2 , h2=.5 ) {
  p1 <- 1+.5*h2*qnorm(1-prev)^2/2
  p2 <- (prev*(1-prev))/(h2*(dnorm(qnorm(1-prev)))^2)
  nrel=(p2*2-p1)/(.5/acc^2-p1)
  nrel
}



## PGS Stats
figPath <- "~/"
data <- fread("Stats.txt" )[ 1:20,]

data$train_rg_test[ data$train_rg_test > 1 ] <- 1

data$pgs_r2l_eSham_40k <- Sham_r2l_cc_sample( data$train_nCa_eff, data$train_nCo_eff, 4e4, data$prev, data$train_h2_snp, data$h2_snp, data$train_rg_test )
data$pgs_r2l_eSham_60k <- Sham_r2l_cc_sample( data$train_nCa_eff, data$train_nCo_eff, 6e4, data$prev, data$train_h2_snp, data$h2_snp, data$train_rg_test )
data$pgs_r2l_eSham_100k <- Sham_r2l_cc_sample( data$train_nCa_eff, data$train_nCo_eff, 1e5, data$prev, data$train_h2_snp, data$h2_snp, data$train_rg_test )

data$pgs_npope <- Sham_r2l_to_n( Sham_r2l_cc_sample( data$train_nCa_eff, data$train_nCo_eff, 6e4, data$prev, data$train_h2_snp ), 6e4, data$prev, data$train_h2_snp  )
data$pgs_m_e <- Sham_r2l_to_m( data$train_nCa_eff, data$train_nCo_eff, data$pgs_r2l_rand, data$prev, data$train_h2_snp, data$h2_snp, data$train_rg_test )
data$pgs_r2l_portabilityFactor <- ( data$h2_snp/data$train_h2_snp ) * data$train_rg_test^2
data$pgs_r2gg_portabilityFactor <- data$train_rg_test^2

data$fgrs_nsibe_h2_real <- acc_to_se( data$fgrs_eacc_all, data$prev, data$h2_real)
data$fgrs_nsibe_w_h2_real <- acc_to_se( data$fgrs_eacc_w_all, data$prev, data$h2_real )

data$ps_nsibe_h2_real <- acc_to_se( data$fgrs_eacc_ps, data$prev, data$h2_real )
data$ps_nsibe_w_h2_real <- acc_to_se( data$fgrs_eacc_w_ps, data$prev, data$h2_real )

data$fgrs_r2l_e_h2_real <- data$h2_real * data$fgrs_eacc_all^2
data$fgrs_r2l_e_w_h2_real <- data$h2_real * data$fgrs_eacc_w_all^2 

data$ps_r2l_e_h2_real <- data$h2_real * data$fgrs_eacc_ps^2
data$ps_r2l_e_w_h2_real <- data$h2_real * data$fgrs_eacc_w_ps^2 

data$pgs_fgrs_rand_cor_e_h2_real <- sqrt( data$pgs_r2l_rand*data$fgrs_r2l_rand ) / data$h2_real
data$pgs_fgrs_ipw_cor_e_h2_real <- sqrt( data$pgs_r2l_cc_wt*data$fgrs_r2l_cc_wt ) / data$h2_real
data$pgs_fgrs_rand_cor_e_h2_halfReal <- sqrt( data$pgs_r2l_rand*data$fgrs_r2l_rand ) / ( 0.5*data$h2_real )
data$pgs_fgrs_ipw_cor_e_h2_halfReal <- sqrt( data$pgs_r2l_cc_wt*data$fgrs_r2l_cc_wt ) / ( 0.5*data$h2_real )

data$h2_rand_e_pgs_fgrs_cor <- sqrt( data$pgs_r2l_rand*data$fgrs_r2l_rand ) / data$pgs_fgrs_cor
data$h2_ipw_e_pgs_fgrs_cor <- sqrt( data$pgs_r2l_cc_wt*data$fgrs_r2l_cc_wt ) / data$pgs_fgrs_cor
data$h2_rand_e_pgs_fgrs_cor_x2 <- 2*sqrt( data$pgs_r2l_rand*data$fgrs_r2l_rand ) / data$pgs_fgrs_cor
data$h2_ipw_e_pgs_fgrs_cor_x2 <- 2*sqrt( data$pgs_r2l_cc_wt*data$fgrs_r2l_cc_wt ) / data$pgs_fgrs_cor

data$pgs_rand_nsibe_h2_real <- acc_to_se( sqrt( data$pgs_r2l_rand/data$h2_real ), data$prev, data$h2_real )
data$pgs_ipw_nsibe_h2_real <- acc_to_se( sqrt( data$pgs_r2l_cc_wt/data$h2_real ), data$prev, data$h2_real )

data$fgrs_rand_npope_60k_in <- Sham_r2l_to_n( data$fgrs_r2l_rand, 6e4, data$prev, data$h2_snp  )
data$fgrs_ipw_npope_60k_in <- Sham_r2l_to_n( data$fgrs_r2l_cc_wt, 6e4, data$prev, data$h2_snp  )
data$fgrs_rand_npope_60k_out <- Sham_r2l_to_n( data$fgrs_r2l_rand, 6e4, data$prev, data$train_h2_snp, data$h2_snp, data$train_rg_test  )
data$fgrs_ipw_npope_60k_out <- Sham_r2l_to_n( data$fgrs_r2l_cc_wt, 6e4, data$prev, data$train_h2_snp, data$h2_snp, data$train_rg_test  )

fwrite( data, file=paste0( figPath, "Stats_230122_Augmented.txt" ), sep='\t' )

data$pch <- c( rep(16,10), rep(17,10) )
data$col <- c( gg_color_hue(5)[ c( 1,1,1,2,2,2,3,4,5,5 ) ], gg_color_hue(5)[ c( 1,1,1,2,2,2,3,4,5,5 ) ] )
data$good <- c( 0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,1,1,1,0,1 )
	
	
p.width <- 4
p.height <- 3
	
# PGS performance vs. Expected / All / Random

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_PGS_o_vs_e_r2l_allTraits_rand.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25
#			, mfrow=c(2,1)
			)

	plot( plot.data$pgs_r2l_eSham_60k, plot.data$pgs_r2l_rand, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.05 ),
		ylim=c( 0, 0.05 ),
#		xlab=expression( E( R[Liab.]^2 ~ '|' ~ 'N,' ~ h[SNP]^2 ~ ', ' ~ rho[G] ~ ', k' ~ ', m=60,000, p=1') ),
		xlab=expression( 'Expected' ~ R[Liab.]^2 ~ '(PGS)' ),
		ylab=expression( 'Observed' ~ R[Liab.]^2 ~ '(PGS)' ),
		yaxt='n'
		)
		abline( 0,1 )	

	axis( 2, at=seq( 0,0.05,by=0.01 ), las=2 )

dev.off()

# PGS performance vs. Expected / All / IPW

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_PGS_o_vs_e_r2l_allTraits_ipw.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25
#			, mfrow=c(2,1)
			)

	plot( plot.data$pgs_r2l_eSham_60k, plot.data$pgs_r2l_cc_wt, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.05 ),
		ylim=c( 0, 0.05 ),
#		xlab=expression( E( R[Liab.]^2 ~ '|' ~ 'N,' ~ h[SNP]^2 ~ ', ' ~ rho[G] ~ ', k' ~ ', m=60,000, p=1') ),
		xlab=expression( 'Expected' ~ R[Liab.]^2 ~ '(PGS)' ),
		ylab=expression( 'Observed' ~ R[Liab.]^2 ~ '(PGS)' ),
		yaxt='n'
		)
		abline( 0,1 )	

	axis( 2, at=seq( 0,0.05,by=0.01 ), las=2 )

dev.off()

# FGRS performance vs. Expected / All

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_FGRS_o_vs_e_r2l_allTraits_h2_real_rand.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$fgrs_r2l_e_w_h2_real, plot.data$fgrs_r2l_rand, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.05 ),
		ylim=c( 0, 0.05 ),
#		xlab=expression( E( R[Liab.]^2 ~ '|' ~ 'Weighted' ~ N[sibe] ~ ', ' ~ h^2 ~ ', k') ),
		xlab=expression( 'Expected' ~ R[Liab.]^2 ~ '(FGRS)' ),
		ylab=expression( 'Observed' ~ R[Liab.]^2 ~ '(FGRS)' ),
		yaxt='n'
		)
	abline( 0,1 )	

	axis( 2, at=seq( 0,0.05,by=0.01 ), las=2 )
	
	legend( 'topleft', 
			legend=c( "ADHD", "ASD", "BPD", "MDD", "SCZ", "iPSYCH 2012", "iPSYCH 2015i", "", "", "" ), 
			pch=c( rep(15,5),16,17,NA,NA,NA ), 
			col=c( gg_color_hue(5)[ 1:5 ],1,1,NA,NA,NA ),
			ncol=2,
			bty='n', cex=0.6 )

dev.off()

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_FGRS_o_vs_e_r2l_allTraits_h2_real_ipw.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$fgrs_r2l_e_w_h2_real, plot.data$fgrs_r2l_cc_wt, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.05 ),
		ylim=c( 0, 0.05 ),
#		xlab=expression( E( R[Liab.]^2 ~ '|' ~ 'Weighted' ~ N[sibe] ~ ', ' ~ h^2 ~ ', k') ),
		xlab=expression( 'Expected' ~ R[Liab.]^2 ~ '(FGRS)' ),
		ylab=expression( 'Observed' ~ R[Liab.]^2 ~ '(FGRS)' ),
		yaxt='n'
		)
	abline( 0,1 )	

	axis( 2, at=seq( 0,0.05,by=0.01 ), las=2 )
	
	legend( 'topleft', 
			legend=c( "ADHD", "ASD", "BPD", "MDD", "SCZ", "iPSYCH 2012", "iPSYCH 2015i", "", "", "" ), 
			pch=c( rep(15,5),16,17,NA,NA,NA ), 
			col=c( gg_color_hue(5)[ 1:5 ],1,1,NA,NA,NA ),
			ncol=2,
			bty='n', cex=0.6 )

dev.off()

# Expected cor vs Observed cor / random smple performance

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_Cor_FGRS_PGS_expected_rand.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$pgs_fgrs_rand_cor_e_h2_real, plot.data$pgs_fgrs_cor_rand, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.10 ),
		ylim=c( 0, 0.10 ),
		xlab=expression( Expected ~ rho[PGS ~ ',' ~ FGRS] ),
		ylab=expression( Observed ~ rho[PGS ~ ',' ~ FGRS] ),
		yaxt='n'
		)
	axis( 2, at=seq(0,0.15,0.05), labels=seq(0,0.15,0.05), las=2 )	
	abline( 0,1 )

#	legend( 'bottomright', 
#			legend=c( "ADHD", "ASD", "BPD", "MDD", "SCZ", "iPSYCH 2012", "iPSYCH 2015i", "", "", "" ), 
#			pch=c( rep(15,5),16,17,NA,NA,NA ), 
#			col=c( gg_color_hue(5)[ 1:5 ],1,1,NA,NA,NA ),
#			ncol=2,
#			bty='n', cex=0.8 )

dev.off()

# Expected cor vs Observed cor / IPW performance

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_Cor_FGRS_PGS_expected_ipw.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$pgs_fgrs_ipw_cor_e_h2_real, plot.data$pgs_fgrs_cor_rand, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.10 ),
		ylim=c( 0, 0.10 ),
		xlab=expression( Expected ~ rho[PGS ~ ',' ~ FGRS] ),
		ylab=expression( Observed ~ rho[PGS ~ ',' ~ FGRS] ),
		yaxt='n'
		)
	axis( 2, at=seq(0,0.15,0.05), labels=seq(0,0.15,0.05), las=2 )	
	abline( 0,1 )

dev.off()

# Expected cor vs Observed cor / random sample performance / E (half h2)

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_Cor_FGRS_PGS_expected_rand_halfh2.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$pgs_fgrs_rand_cor_e_h2_halfReal, plot.data$pgs_fgrs_cor_rand, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.125 ),
		ylim=c( 0, 0.125 ),
		xlab=expression( Expected ~ rho[PGS ~ ',' ~ FGRS] ),
		ylab=expression( Observed ~ rho[PGS ~ ',' ~ FGRS] ),
		yaxt='n'
		)
	axis( 2, at=seq(0,0.12,0.02), labels=c(0,'',0.04,'',0.8,'',0.12), las=2 )	
	abline( 0,1 )

dev.off()

# Expected cor vs Observed cor / IPW performance / E (half h2)

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_Cor_FGRS_PGS_expected_ipw_halfh2.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$pgs_fgrs_ipw_cor_e_h2_halfReal, plot.data$pgs_fgrs_cor_rand, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 0.125 ),
		ylim=c( 0, 0.125 ),
		xlab=expression( Expected ~ rho[PGS ~ ',' ~ FGRS] ),
		ylab=expression( Observed ~ rho[PGS ~ ',' ~ FGRS] ),
		yaxt='n'
		)
	axis( 2, at=seq(0,0.12,0.02), labels=c(0,'',0.04,'',0.8,'',0.12), las=2 )	
	abline( 0,1 )

dev.off()

## Sample size equivalents

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_PGS_npope_vs_nsibe_rand.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$pgs_npope, plot.data$pgs_rand_nsibe_h2_real, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 2000000 ),
		ylim=c( 0, 5 ),
		xlab=expression( Observed ~ N[ 'Pop. GWAS Equiv.' ] ~ '(PGS)' ),
		ylab=expression( Expected ~ N[ 'Full Sib. Equiv.' ] ~ '(FGRS)' ),
		yaxt='n', xaxt='n'
		)
	axis( 2, at=1:5, labels=1:5, las=2 )	
	axis( 1, at=seq(0,2250000,by=250000), 
			labels=c( '0','','500,000','','1,000,000','','1,500,000','','2,000,000',''), 
			cex.axis=0.7, las=1 )	
		
	legend( 'topleft', 
			legend=c( "", "", "", "iPSYCH 2012", "iPSYCH 2015i", "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 
			pch=c( NA,NA,NA,16,17,rep(15,5) ), 
			col=c( NA,NA,NA,1,1,gg_color_hue(5)[ 1:5 ] ),
			ncol=2,
			bty='n', cex=0.6 )

dev.off()

## Sample size equivalents

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_PGS_npope_vs_nsibe_ipw.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$pgs_npope, plot.data$pgs_ipw_nsibe_h2_real, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 2000000 ),
		ylim=c( 0, 5 ),
		xlab=expression( Observed ~ N[ 'Pop. GWAS Equiv.' ] ~ '(PGS)' ),
		ylab=expression( Expected ~ N[ 'Full Sib. Equiv.' ] ~ '(FGRS)' ),
		yaxt='n', xaxt='n'
		)
	axis( 2, at=1:5, labels=1:5, las=2 )	
	axis( 1, at=seq(0,2250000,by=250000), 
			labels=c( '0','','500,000','','1,000,000','','1,500,000','','2,000,000',''), 
			cex.axis=0.7, las=1 )	
		
	legend( 'topleft', 
			legend=c( "", "", "", "iPSYCH 2012", "iPSYCH 2015i", "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 
			pch=c( NA,NA,NA,16,17,rep(15,5) ), 
			col=c( NA,NA,NA,1,1,gg_color_hue(5)[ 1:5 ] ),
			ncol=2,
			bty='n', cex=0.6 )

dev.off()

## FGRS to PGS sample size

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_FGRS_nsibe_vs_npope_inPopulation_rand.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$fgrs_nsibe_w_h2_real, plot.data$fgrs_rand_npope_60k_in, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 4),
		ylim=c( 0, 2000000 ),
		xlab=expression( Observed ~ N[ 'Full Sib. Equiv.' ] ~ '(FGRS)' ),
		ylab=expression( Expected ~ N[ 'Pop. GWAS Equiv.' ] ~ '(PGS)' ),
		yaxt='n', xaxt='n'
		)
	axis( 1, at=0:10, labels=0:10, las=1 )	
	axis( 2, at=seq(0,2250000,by=250000), 
			labels=c( '0','','5e5','','1e6','','1.5e6','','2e6',''), 
			cex.axis=0.8, las=1 )	
	
	legend( 'topleft', 
			legend=c( "", "", "", "iPSYCH 2012", "iPSYCH 2015i", "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 
			pch=c( NA,NA,NA,16,17,rep(15,5) ), 
			col=c( NA,NA,NA,1,1,gg_color_hue(5)[ 1:5 ] ),
			ncol=2,
			bty='n', cex=0.6 )

dev.off()

plot.data <- data[ data$good == 1, ]
pdf( paste0( figPath, "Fig_Main_FGRS_nsibe_vs_npope_inPopulation_ipw.pdf" ), height=p.height, width=p.width )

	par( mar=c(3.5,3.5,1,1), mgp=c( 2,0.5,0 ), tcl=-0.25 )

	plot( plot.data$fgrs_nsibe_w_h2_real, plot.data$fgrs_ipw_npope_60k_in, 
		pch=plot.data$pch, 
		col=plot.data$col, 
		cex=1,
		xlim=c( 0, 4),
		ylim=c( 0, 2000000 ),
		xlab=expression( Observed ~ N[ 'Full Sib. Equiv.' ] ~ '(FGRS)' ),
		ylab=expression( Expected ~ N[ 'Pop. GWAS Equiv.' ] ~ '(PGS)' ),
		yaxt='n', xaxt='n'
		)
	axis( 1, at=0:10, labels=0:10, las=1 )	
	axis( 2, at=seq(0,2250000,by=250000), 
			labels=c( '0','','5e5','','1e6','','1.5e6','','2e6',''), 
			cex.axis=0.8, las=1 )	
	
	legend( 'topleft', 
			legend=c( "", "", "", "iPSYCH 2012", "iPSYCH 2015i", "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 
			pch=c( NA,NA,NA,16,17,rep(15,5) ), 
			col=c( NA,NA,NA,1,1,gg_color_hue(5)[ 1:5 ] ),
			ncol=2,
			bty='n', cex=0.6 )

dev.off()

