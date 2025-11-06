rm(list = ls())

## Packages

library( data.table )
library( ggplot2 )
library( knitr )
library( fmsb )
library( sandwich )

## Load data

load( '~/MDD_heterogeneity/backup/data/psychData_AJS_230207.RData' )

## Fixed Parameters

N_cohort=1472762 # Pedersen2017

#c( "SCZ", "BP", "MDD", "ASD", "ADHD" )
life_time_risk <- NULL
life_time_risk$F9200 <- 0.03235
life_time_risk$F8100 <- 0.01965
life_time_risk$F3100 <- 0.0158
life_time_risk$F3300 <- 0.12285
life_time_risk$F2100 <- 0.01745

pNames <- c( "ADHD", "ASD", "BPD", "MDD", "SCZ" )

results <- NULL
ors <- NULL
ses <- NULL
ses.sw <- NULL

for ( c in c( 'i2012', 'i2015' ) ) {
  for ( p in c( 'F9200', 'F8100', 'F3100', 'F3300', 'F2100' ) ) {

	print( paste0( c, '_', p ) )
   
	i <- which( names( dat ) == paste0( c, "_", p ) )
	data <- dat[[ i ]]
  data <- dat[[ i ]][ !is.na(PC1) & !is.na(FGRS) & !is.na(PGS) ] 

	n <- length( data$dx )
	nCase <- data[ ,sum( dx==1 ) ]
	nControl <- data[ ,sum( dx==0 ) ]

	# scale using weights 
	P <- data[ ,mean( dx, na.rm=T ) ]
	P <- data[ dx==1, sum(1/prob) ] / data[ ,sum(1/prob) ]
	K <- as.numeric( life_time_risk[ p ] )
	wv <- (1-P)*K/( P*(1-K) ) 		#weighting factor
	wt <- data$dx + 1
	wt[ wt==2 ] <- wv 			#weighting array 
	wt.o <- wt
	wt <- wt / data$prob
	wt <- wt / max(wt)
 
	# Null model
	glmv0 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
			weights=wt, family=binomial )
	glmv1 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
			PGS, weights=wt, family=binomial ) 
	glmv2 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
			FGRS, weights=wt, family=binomial ) 
	glmv3 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
			FH, weights=wt, family=binomial ) 
	glmv4 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
			FGRS + PGS, weights=wt, family=binomial )
	glmv5 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
			FH + PGS, weights=wt, family=binomial )
	glmv6 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
			FGRS + FH, weights=wt, family=binomial )
	glmv7 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
			FGRS + FH + PGS, weights=wt, family=binomial )

	coef <- function( x ) { summary( x )$coefficients }
	t.ors <- c( c, p, 
			coef(glmv1)[14,1], coef(glmv4)[15,1], coef(glmv5)[15,1], coef(glmv7)[16,1], 
			coef(glmv2)[14,1], coef(glmv4)[14,1], coef(glmv6)[14,1], coef(glmv7)[14,1],
			coef(glmv3)[14,1], coef(glmv5)[14,1], coef(glmv6)[15,1], coef(glmv7)[15,1] ) 
	t.ses <- c( c, p, 
			coef(glmv1)[14,2], coef(glmv4)[15,2], coef(glmv5)[15,2], coef(glmv7)[16,2], 
			coef(glmv2)[14,2], coef(glmv4)[14,2], coef(glmv6)[14,2], coef(glmv7)[14,2],
			coef(glmv3)[14,2], coef(glmv5)[14,2], coef(glmv6)[15,2], coef(glmv7)[15,2] ) 
	t.ses.sw <- c( c, p, 
			sqrt(diag(sandwich(glmv1)))[14], sqrt(diag(sandwich(glmv4)))[15], 
			sqrt(diag(sandwich(glmv5)))[15], sqrt(diag(sandwich(glmv7)))[16], 
			sqrt(diag(sandwich(glmv2)))[14], sqrt(diag(sandwich(glmv4)))[14], 
			sqrt(diag(sandwich(glmv6)))[14], sqrt(diag(sandwich(glmv7)))[14],
			sqrt(diag(sandwich(glmv3)))[14], sqrt(diag(sandwich(glmv5)))[14], 
			sqrt(diag(sandwich(glmv6)))[15], sqrt(diag(sandwich(glmv7)))[15] ) 

	R2_glmv0 <- NagelkerkeR2( glmv0 )$R2
	R2_glmv1 <- NagelkerkeR2( glmv1 )$R2
	R2_glmv2 <- NagelkerkeR2( glmv2 )$R2
	R2_glmv3 <- NagelkerkeR2( glmv3 )$R2
	R2_glmv4 <- NagelkerkeR2( glmv4 )$R2
	R2_glmv5 <- NagelkerkeR2( glmv5 )$R2
	R2_glmv6 <- NagelkerkeR2( glmv6 )$R2
	R2_glmv7 <- NagelkerkeR2( glmv7 )$R2

	stat.test <-  data.frame( rbind(
				c( "I", data.frame( anova( glmv0, glmv1, test="LRT" ) )[1,] ),	
				c( "PGS", data.frame( anova( glmv0, glmv1, test="LRT" ) )[2,] ),
        c( "FGRS", data.frame( anova( glmv0, glmv2, test="LRT" ) )[2,] ),
        c( "FH", data.frame( anova( glmv0, glmv3, test="LRT" ) )[2,] ),
        c( "PGS_FGRS", data.frame( anova( glmv1, glmv4, test="LRT" ) )[2,] ),
        c( "PGS_FH", data.frame( anova( glmv1, glmv5, test="LRT" ) )[2,] ),
        c( "FGRS_PGS", data.frame( anova( glmv2, glmv4, test="LRT" ) )[2,] ),
        c( "FGRS_FH", data.frame( anova( glmv2, glmv6, test="LRT" ) )[2,] ),
        c( "FH_PGS", data.frame( anova( glmv3, glmv5, test="LRT" ) )[2,] ),
        c( "FH_FGRS", data.frame( anova( glmv3, glmv6, test="LRT" ) )[2,] ),
        c( "PGS_FGRS_FH", data.frame( anova( glmv4, glmv7, test="LRT" ) )[2,] ),
        c( "PGS_FH_FGRS", data.frame( anova( glmv5, glmv7, test="LRT" ) )[2,] ),
        c( "FGRS_FH_PGS", data.frame( anova( glmv6, glmv7, test="LRT" ) )[2,] )
				) )

	stat.test$p <- pchisq( as.numeric(stat.test$Deviance), as.numeric(stat.test$Df), lower.tail = F )

	stat.test[,3] <- round( as.numeric( stat.test[,3] ),2 )
	stat.test[,5] <- round( as.numeric( stat.test[,5] ),2 )
	stat.test[,6] <- signif( as.numeric( stat.test[,6] ),3 )
	stat.test[,7] <- signif( as.numeric( stat.test[,7] ),3 )

	n <- length( data$dx )
	nCase <- sum( data$dx == 1 )
	nControl <- sum( data$dx == 0 )

	results <- rbind( results, c( c, p, n, nCase, nControl, 
					R2_glmv0, R2_glmv1, R2_glmv2, R2_glmv3, R2_glmv4, R2_glmv5, R2_glmv6, R2_glmv7,
					stat.test$p[ 2:13 ] ) )
	ors <- rbind( ors, t.ors )
	ses <- rbind( ses, t.ses )
	ses.sw <- rbind( ses.sw, t.ses.sw )

	fwrite( stat.test, 
			file=paste0( "/home/anscho/MDD_heterogeneity/backup/data/", c, "_", 
					p, "_modelStats_logistic_withCovs_weighted_230217.txt"  ), 
			row.names=F, col.names=T, sep='\t' )

}
}

results <- data.table( results )
ors <- data.table( ors )
ses <- data.table( ses )
ses.sw <- data.table( ses.sw )

for ( i in 3:25 ) {
	results[[i]] <- as.numeric( results[[i]] )
}
for ( i in 3:14 ) {
	ors[[i]] <- as.numeric( ors[[i]] )
	ses[[i]] <- as.numeric( ses[[i]] )
	ses.sw[[i]] <- as.numeric( ses.sw[[i]] )
}

names( results ) <- c( 'cohort', 'trait', 'n', 'n_case', 'n_control', 
			'r2N_covs', 'r2N_PGS', 'r2N_FGRS', 'r2N_FH', 
			'r2N_PGS_FGRS', 'r2N_PGS_FH', 'r2N_FGRS_FH', 
			'r2N_PGS_FGRS_FH',
			paste0( "p_", stat.test$V1[2:13] ) )

names( ors ) <- names( ses ) <- names( ses.sw ) <- c( 'cohort', 'trait', 
			'pgs', 'pgs_fgrs', 'pgs_fh', 'pgs_both', 
			'fgrs', 'fgrs_pgs', 'fgrs_fh', 'fgrs_both', 
			'fh', 'fh_pgs', 'fh_fgrs', 'fh_both' )

write.table( results, "/home/anscho/MDD_heterogeneity/backup/data/R2N_logistic_withCovs_230217_unrel_weights.txt", 
		row.names=F, col.names=T, quote=F )
write.table( ors, "/home/anscho/MDD_heterogeneity/backup/data/OR_logistic_withCovs_230217_unrel_weights.txt", 
		row.names=F, col.names=T, quote=F )
write.table( ses, "/home/anscho/MDD_heterogeneity/backup/data/SE_logistic_withCovs_230217_unrel_weights.txt", 
		row.names=F, col.names=T, quote=F )
write.table( ses.sw, "/home/anscho/MDD_heterogeneity/backup/data/SE_sandwich_logistic_withCovs_230217_unrel_weights.txt", 
		row.names=F, col.names=T, quote=F )


## Make some plots

results <- fread( "/home/anscho/MDD_heterogeneity/backup/data/R2N_logistic_withCovs_230217_unrel_weights.txt" )
ors <- fread( "/home/anscho/MDD_heterogeneity/backup/data/OR_logistic_withCovs_230217_unrel_weights.txt" )
ses <- fread( "/home/anscho/MDD_heterogeneity/backup/data/SE_logistic_withCovs_230217_unrel_weights.txt" )
ses.sw <- fread( "/home/anscho/MDD_heterogeneity/backup/data/SE_sandwich_logistic_withCovs_230217_unrel_weights.txt" )

gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

p.t <- 0.05/(11*5*2)
ps <- results
for ( i in 14:25 ) { ps[[ i ]] <- 1*( ps[[ i ]] > p.t ) }

pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/r2N_logistic_ipsych2012_withCovs_unrel_weights_230217.pdf", 
      width=5, height=5 )
par( mgp=c(2,0.5,0), tcl=-0.25 )

rows <- 1:5
r2l.covs <- t( as.matrix( cbind( results[ rows, 6 ], results[ rows, 6 ], 
				results[ rows, 6 ], results[ rows, 6 ], 
				results[ rows, 6 ], results[ rows, 6 ] ) ) )
r2l <- t( as.matrix( results[ rows, c(7:11,13) ] ) )
r2l <- r2l - r2l.covs
p.r2l <- ps[ rows, c(2,14:19,21,23:24) ]

mids <- barplot( r2l, beside=T, col=gg_color_hue(6), names.arg=c( "ADHD","ASD","BPD","MDD","SCZ" ),
		ylim=c( 0,0.14 ), ylab=expression( Nagelkerke ~ r^2 ) ) 
legend( "top", legend=c( "PGS", "FGRS", "FH", "PGS+FGRS", "PGS+FH", "PGS+FGRS+FH" ),
		col=gg_color_hue(6), pch=15, bty='n', ncol=3, pt.cex=1.5, cex=0.66 )

for ( i in 1:5 ) {
	catchT <- F 
	if ( p.r2l$p_PGS_FGRS_FH[i] == 1 ) {
		segments( mids[6,i], r2l[6,i]+0.003, mids[4,i], r2l[6,i]+0.003 )
		segments( mids[4,i], r2l[6,i]+0.003, mids[4,i], r2l[6,i]+0.002 ) 
		segments( mids[6,i], r2l[6,i]+0.003, mids[6,i], r2l[6,i]+0.002 )
		text( mids[4,i]+(mids[6,i]-mids[4,i])/2, r2l[6,i]+0.003, "n.s.", cex=0.66, pos=3, offset=0.1 )  
		catchT <- T 
	} 
	if ( p.r2l$p_PGS_FH_FGRS[i] == 1 ) {
		segments( mids[6,i], r2l[6,i]+0.0015, mids[5,i], r2l[6,i]+0.0015 )
		segments( mids[5,i], r2l[6,i]+0.0015, mids[5,i], r2l[6,i]+0.0005 ) 
		segments( mids[6,i], r2l[6,i]+0.0015, mids[6,i], r2l[6,i]+0.0005 )
		if ( catchT == F ) {
			text( mids[4,i]+(mids[6,i]-mids[4,i])/2, r2l[6,i]+0.0015, "n.s.", cex=0.66, pos=3, offset=0.1 )
	  	}
	} 
}

dev.off()


pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/r2N_logistic_ipsych2015i_withCovs_unrel_weights_230217.pdf", 
      width=5, height=5 )
par( mgp=c(2,0.5,0), tcl=-0.25 )

rows <- 6:10
r2l.covs <- t( as.matrix( cbind( results[ rows, 6 ], results[ rows, 6 ], 
				results[ rows, 6 ], results[ rows, 6 ], 
				results[ rows, 6 ], results[ rows, 6 ] ) ) )
r2l <- t( as.matrix( results[ rows, c(7:11,13) ] ) )
r2l <- r2l - r2l.covs
p.r2l <- ps[ rows, c(2,14:19,21,23:24) ]

mids <- barplot( r2l, beside=T, col=gg_color_hue(6), names.arg=c( "ADHD","ASD","BPD","MDD","SCZ" ),
		ylim=c( 0,0.14 ), ylab=expression( Nagelkerke ~ r^2 ) ) 
#legend( "top", legend=c( "PGS", "FGRS", "FH", "PGS+FGRS", "PGS+FH", "PGS+FGRS+FH" ),
#		col=gg_color_hue(6), pch=15, bty='n', ncol=3, pt.cex=1.5, cex=0.66 )

for ( i in 1:5 ) {
	catchT <- F
	if ( p.r2l$p_PGS_FGRS_FH[i] == 1 ) {
		segments( mids[6,i], r2l[6,i]+0.003, mids[4,i], r2l[6,i]+0.003 )
		segments( mids[4,i], r2l[6,i]+0.003, mids[4,i], r2l[6,i]+0.002 ) 
		segments( mids[6,i], r2l[6,i]+0.003, mids[6,i], r2l[6,i]+0.002 )
		text( mids[4,i]+(mids[6,i]-mids[4,i])/2, r2l[6,i]+0.003, "n.s.", cex=0.66, pos=3, offset=0.1 )
		catchT <- T  
	} 
	if ( p.r2l$p_PGS_FH_FGRS[i] == 1 ) {
		segments( mids[6,i], r2l[6,i]+0.0015, mids[5,i], r2l[6,i]+0.0015 )
		segments( mids[5,i], r2l[6,i]+0.0015, mids[5,i], r2l[6,i]+0.0005 ) 
		segments( mids[6,i], r2l[6,i]+0.0015, mids[6,i], r2l[6,i]+0.0005 )
		if ( catchT == F ) {
			text( mids[4,i]+(mids[6,i]-mids[4,i])/2, r2l[6,i]+0.0015, "n.s.", cex=0.66, pos=3, offset=0.1 )
	  	}
	} 
}

dev.off()

## ORs

# 2012
rows <- 1:5
p.lor <- ors[ rows, 3:14 ]
p.ors <- exp( p.lor )
p.se <- ses.sw[ rows, 3:14 ]
p.u95 <- exp( p.lor + 1.96*p.se )
p.l95 <- exp( p.lor - 1.96*p.se )

pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/OR_logistic_ipsych2012_withCovs_unrel_weights_230217.pdf", width=9, height=6 )

cols <- c( rep( gg_color_hue(6)[1], 4 ), rep( gg_color_hue(6)[2], 4 ), rep( gg_color_hue(6)[3], 4 ) )
ltys <- c( 1,3,4,5, 1,2,4,5, 1,2,3,5 ) 
par( mar=c(2,3,1,1), mgp=c(2,0.5,0), tcl=-0.25 )
plot( 1, 1, type='n',
	xaxt='n', xlim=c(0,65), ylim=c(0.8,1.85),
	ylab="OR per s.d.", xlab='' )
polygon( c( 13,13,26,26 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )
polygon( c( 39,39,52,52 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )

for ( i in 1:5 ) {
	points( (i+12*(i-1)):((i-1)+12*(i)), p.ors[i,], 
		col=cols, pch=18, cex=1 )
	segments( (i+12*(i-1)):((i-1)+12*(i)), as.numeric( p.u95[i,] ), 
			(i+12*(i-1)):((i-1)+12*(i)), as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5, lty=ltys )
	segments( (i+12*(i-1)):((i-1)+12*(i))-0.25, as.numeric( p.l95[i,] ), 
			(i+12*(i-1)):((i-1)+12*(i))+0.25, as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5 )
	segments( (i+12*(i-1)):((i-1)+12*(i))-0.25, as.numeric( p.u95[i,] ), 
			(i+12*(i-1)):((i-1)+12*(i))+0.25, as.numeric( p.u95[i,] ), 
			col=cols, lwd=1.5 )
}

#abline( v=c(13,26,39,52), col='grey' )
abline( h=1 )
axis( 1, at=c( 6.5, 19.5, 32.5, 45.5, 58.5 ), lab=c("ADHD","ASD","BPD","MDD","SCZ") )

legend( 'topleft', legend=c( "PGS", "No adj.", "FH adj.", "FGRS","PGS adj.", "Full adj.", "FH","FGRS adj.","" ), 
	lty=c( 1,1,4,1,2,5,1,3,0 ), col=c(cols[1],1,1,cols[5],1,1,cols[9],1,1), ncol=4, bty='n', lwd=1.5 ) 
#legend( 'top', legend=rep( "",9 ), 
#	pch=0, col=c(cols[1],0,0,cols[5],0,0,cols[9],0,0), ncol=9, bty='n' ) 

dev.off()

# 2015i
rows <- 6:10
p.lor <- ors[ rows, 3:14 ]
p.ors <- exp( p.lor )
p.se <- ses.sw[ rows, 3:14 ]
p.u95 <- exp( p.lor + 1.96*p.se )
p.l95 <- exp( p.lor - 1.96*p.se )

pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/OR_logistic_ipsych2015i_withCovs_unrel_weights_230217.pdf", width=9, height=6 )

cols <- c( rep( gg_color_hue(6)[1], 4 ), rep( gg_color_hue(6)[2], 4 ), rep( gg_color_hue(6)[3], 4 ) )
ltys <- c( 1,3,4,5, 1,2,4,5, 1,2,3,5 ) 
par( mar=c(2,3,1,1), mgp=c(2,0.5,0), tcl=-0.25 )
plot( 1, 1, type='n',
	xaxt='n', xlim=c(0,65), ylim=c(0.8,1.85),
	ylab="OR per s.d.", xlab='' )
polygon( c( 13,13,26,26 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )
polygon( c( 39,39,52,52 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )

for ( i in 1:5 ) {
	points( (i+12*(i-1)):((i-1)+12*(i)), p.ors[i,], 
		col=cols, pch=18, cex=1  )
	segments( (i+12*(i-1)):((i-1)+12*(i)), as.numeric( p.u95[i,] ), 
			(i+12*(i-1)):((i-1)+12*(i)), as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5, lty=ltys )
	segments( (i+12*(i-1)):((i-1)+12*(i))-0.25, as.numeric( p.l95[i,] ), 
			(i+12*(i-1)):((i-1)+12*(i))+0.25, as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5 )
	segments( (i+12*(i-1)):((i-1)+12*(i))-0.25, as.numeric( p.u95[i,] ), 
			(i+12*(i-1)):((i-1)+12*(i))+0.25, as.numeric( p.u95[i,] ), 
			col=cols, lwd=1.5 )
}

#abline( v=c(13,26,39,52), col='grey' )
abline( h=1 )
axis( 1, at=c( 6.5, 19.5, 32.5, 45.5, 58.5 ), lab=c("ADHD","ASD","BPD","MDD","SCZ") )

legend( 'topleft', legend=c( "PGS","No adj.","FH adj.","FGRS","PGS adj.","Full adj.","FH","FGRS adj.","" ), 
	lty=c( 1,1,4,1,2,5,1,3,0 ), col=c(cols[1],1,1,cols[5],1,1,cols[9],1,1), ncol=4, bty='n', lwd=1.5 ) 
#legend( 'top', legend=rep( "",9 ), 
#	pch=0, col=c(cols[1],0,0,cols[5],0,0,cols[9],0,0), ncol=9, bty='n' ) 

dev.off()

## ORs / Main Figure

# 2012
rows <- 1:5
p.lor <- ors[ rows, c(3:4,7:8) ]
p.ors <- exp( p.lor )
p.se <- ses.sw[ rows, c(3:4,7:8) ]
p.u95 <- exp( p.lor + 1.96*p.se )
p.l95 <- exp( p.lor - 1.96*p.se )

pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/OR_logistic_ipsych2012_withCovs_unrel_weights_230217_MainText.pdf",
	 width=5, height=5 )

cols <- c( rep( gg_color_hue(6)[1], 2 ), rep( gg_color_hue(6)[2], 2 ) )
ltys <- c( 1,2, 1,2 ) 

par( mar=c(2,3,1,1), mgp=c(2,0.5,0), tcl=-0.25 )
plot( 1, 1, type='n',
	xaxt='n', xlim=c(0,25), ylim=c(0.8,1.85),
	ylab="OR per s.d.", xlab='' )
polygon( c( 5,5,10,10 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )
polygon( c( 15,15,20,20 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )

for ( i in 1:5 ) {
	points( (i+4*(i-1)):((i-1)+4*(i)), p.ors[i,], 
		col=cols, pch=18, cex=1 )
	segments( (i+4*(i-1)):((i-1)+4*(i)), as.numeric( p.u95[i,] ), 
			(i+4*(i-1)):((i-1)+4*(i)), as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5, lty=ltys )
	segments( (i+4*(i-1)):((i-1)+4*(i))-0.25, as.numeric( p.l95[i,] ), 
			(i+4*(i-1)):((i-1)+4*(i))+0.25, as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5 )
	segments( (i+4*(i-1)):((i-1)+4*(i))-0.25, as.numeric( p.u95[i,] ), 
			(i+4*(i-1)):((i-1)+4*(i))+0.25, as.numeric( p.u95[i,] ), 
			col=cols, lwd=1.5 )
}

#abline( v=c(13,26,39,52), col='grey' )
abline( h=1 )
axis( 1, at=c( 2.5, 7.5, 12.5, 17.5, 22.5 ), lab=c("ADHD","ASD","BPD","MDD","SCZ") )

legend( 'bottomleft', legend=c( "PGS", "FGRS", "No adj.", "PGS/FGRS adj." ), 
	lty=c( 1, 1, 1, 2 ), col=c( cols[1], cols[3], 1, 1 ), ncol=2, bty='n', lwd=1.5 ) 
#legend( 'top', legend=rep( "",9 ), 
#	pch=0, col=c(cols[1],0,0,cols[5],0,0,cols[9],0,0), ncol=9, bty='n' ) 

dev.off()

# 2015i
rows <- 6:10
p.lor <- ors[ rows, c(3:4,7:8) ]
p.ors <- exp( p.lor )
p.se <- ses.sw[ rows, c(3:4,7:8) ]
p.u95 <- exp( p.lor + 1.96*p.se )
p.l95 <- exp( p.lor - 1.96*p.se )

pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/OR_logistic_ipsych2015i_withCovs_unrel_weights_230217_MainText.pdf", 
	width=5, height=5 )

cols <- c( rep( gg_color_hue(6)[1], 2 ), rep( gg_color_hue(6)[2], 2 ) )
ltys <- c( 1,2, 1,2 ) 

par( mar=c(2,3,1,1), mgp=c(2,0.5,0), tcl=-0.25 )
plot( 1, 1, type='n',
	xaxt='n', xlim=c(0,25), ylim=c(0.8,1.85),
	ylab="OR per s.d.", xlab='' )
polygon( c( 5,5,10,10 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )
polygon( c( 15,15,20,20 ), c( 0.8,1.85,1.85,0.8 ), col='grey97', border=NA )

for ( i in 1:5 ) {
	points( (i+4*(i-1)):((i-1)+4*(i)), p.ors[i,], 
		col=cols, pch=18, cex=1 )
	segments( (i+4*(i-1)):((i-1)+4*(i)), as.numeric( p.u95[i,] ), 
			(i+4*(i-1)):((i-1)+4*(i)), as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5, lty=ltys )
	segments( (i+4*(i-1)):((i-1)+4*(i))-0.25, as.numeric( p.l95[i,] ), 
			(i+4*(i-1)):((i-1)+4*(i))+0.25, as.numeric( p.l95[i,] ), 
			col=cols, lwd=1.5 )
	segments( (i+4*(i-1)):((i-1)+4*(i))-0.25, as.numeric( p.u95[i,] ), 
			(i+4*(i-1)):((i-1)+4*(i))+0.25, as.numeric( p.u95[i,] ), 
			col=cols, lwd=1.5 )
}

#abline( v=c(13,26,39,52), col='grey' )
abline( h=1 )
axis( 1, at=c( 2.5, 7.5, 12.5, 17.5, 22.5 ), lab=c("ADHD","ASD","BPD","MDD","SCZ") )

legend( 'bottomleft', legend=c( "PGS", "FGRS", "No adj.", "PGS/FGRS adj." ), 
	lty=c( 1, 1, 1, 2 ), col=c( cols[1], cols[3], 1, 1 ), ncol=2, bty='n', lwd=1.5 ) 
#legend( 'top', legend=rep( "",9 ), 
#	pch=0, col=c(cols[1],0,0,cols[5],0,0,cols[9],0,0), ncol=9, bty='n' ) 

dev.off()


