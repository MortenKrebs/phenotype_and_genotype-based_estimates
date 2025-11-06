rm(list = ls())

## Packages

library(data.table)
library(ggplot2)
library(knitr)

## Load data

load( '~/MDD_heterogeneity/backup/data/psychData_MDK_230127.RData' )

## Fixed Parameters

EurPCs=F
N_cohort=1472762 # Pedersen2017

#c( "SCZ", "BP", "MDD", "ASD", "ADHD" )
# ASD not (0.0023+0.0087)/2  ?
life_time_risk <- NULL
life_time_risk$F9200 <- 0.03235
life_time_risk$F9200i <- 0.03235
life_time_risk$F8100 <- 0.01965
life_time_risk$F8100i <- 0.01965
life_time_risk$F3100 <- 0.0158
life_time_risk$F3300 <- 0.12285
life_time_risk$F2100 <- 0.01745

pNames <- c( "ADHD", "ADHD", "ASD", "ASD", "BPD", "MDD", "SCZ" )

results <- NULL
for ( c in c( 'i2012', 'i2015' ) ) {
  for ( p in c( 'F9200', 'F9200i', 'F8100', 'F8100i', 'F3100', 'F3300', 'F2100' ) ) {
    
    if ( !p=='F2100' ) {
	i <- which( names( dat ) == paste0( c, "_", p ) ) 
	data <- dat[[ i ]][!is.na(PC1)&!is.na(FGRS)&!is.na(PGS_eu)]
    } else {
	i <- which( names( dat ) == paste0( c, "_", p, "_noDanes" ) )
	data <- dat[[ i ]][!is.na(PC1)&!is.na(FGRS)&!is.na(PGS)] 
	data$PGS_eu <- data$PGS 
    }    
    
    if ( EurPCs ) {
    	data <- data[,which(!names(data) %in% paste0("PC",1:10)),with=F]
    	names(data)[which(names(data) %in% paste0("PC",1:10,"_e"))] <- paste0("PC",1:10)
    }

    # R2 on the probit liability scale using a probit model 
    P <- data[ ,mean( dx,na.rm = T ) ]
    P <- data[ dx==1, sum(1/prob) ] / data[ ,sum(1/prob) ]
    K <- as.numeric( life_time_risk[ p ] )
    wv <- (1-P)*K/(P*(1-K)) 		#weighting factor
    wt <- data$dx + 1
    wt[ wt==2 ] <- wv 			#weighting array 
    wt.o <- wt
    wt=wt/data$prob
    wt=wt/max(wt)
    #wt=wt/quantile(wt,.99)
    #wt[wt>1]=1 

    # weighted probit models / non-integer #successes in a binomial glm!
    pmv0 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
                 weights=wt, family=binomial(probit) ) #weighted probit model
    pmv1 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   scale( PGS_eu ), weights=wt, family=binomial(probit) ) 
    pmv2 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   scale( FGRS ), weights=wt, family=binomial(probit) ) 
    pmv3 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   scale( FH ), weights=wt, family=binomial(probit) ) 
    pmv4 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   scale( FGRS ) + scale( PGS_eu ), weights=wt, family=binomial(probit) )
    pmv5 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   scale( FH ) + scale( PGS_eu ), weights=wt, family=binomial(probit) )
    pmv6 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   scale( FGRS ) + scale( FH ), weights=wt, family=binomial(probit) )
    pmv7 <- glm( data=data, dx~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   scale( FGRS ) + scale( FH ) + scale( PGS_eu ), weights=wt, family=binomial(probit) )
    
    vr <- runif( length( data$dx ), 0, 1 ) 	#uniform random values
    vsel <- pmv0$linear.predictors[ vr<wt ]	#select controls and a proportion (wv) of cases
    R2_prob0 <- var(vsel)/(var(vsel)+1) 
    vsel <- pmv1$linear.predictors[ vr<wt ]	#select controls and a proportion (wv) of cases
    R2_prob1 <- var(vsel)/(var(vsel)+1) 
    vsel <- pmv2$linear.predictors[ vr<wt ] 	#select controls and a proportion (wv) of cases
    R2_prob2 <- var(vsel)/(var(vsel)+1) 
    vsel <- pmv3$linear.predictors[ vr<wt ]	#select controls and a proportion (wv) of cases
    R2_prob3 <- var(vsel)/(var(vsel)+1)
    vsel <- pmv4$linear.predictors[ vr<wt ]	#select controls and a proportion (wv) of cases
    R2_prob4 <- var(vsel)/(var(vsel)+1)
    vsel <- pmv5$linear.predictors[ vr<wt ]	#select controls and a proportion (wv) of cases
    R2_prob5 <- var(vsel)/(var(vsel)+1)
    vsel <- pmv6$linear.predictors[ vr<wt ]	#select controls and a proportion (wv) of cases
    R2_prob6 <- var(vsel)/(var(vsel)+1)
    vsel <- pmv7$linear.predictors[ vr<wt ]	#select controls and a proportion (wv) of cases
    R2_prob7 <- var(vsel)/(var(vsel)+1)
    
    stat.test <-  data.frame( rbind(
      c( "I", data.frame( anova( pmv0,pmv1,test = "LRT" ) )[1,] ),	
      c( "PGS", data.frame( anova( pmv0,pmv1,test = "LRT" ) )[2,] ),
      c( "FGRS", data.frame( anova( pmv0,pmv2,test = "LRT" ) )[2,] ),
      c( "FH", data.frame( anova( pmv0,pmv3,test = "LRT" ) )[2,] ),
      c( "PGS_FGRS", data.frame( anova( pmv1,pmv4,test = "LRT" ) )[2,] ),
      c( "PGS_FH", data.frame( anova( pmv1,pmv5,test = "LRT" ) )[2,] ),
      c( "FGRS_PGS", data.frame( anova( pmv2,pmv4,test = "LRT" ) )[2,] ),
      c( "FGRS_FH", data.frame( anova( pmv2,pmv6,test = "LRT" ) )[2,] ),
      c( "FH_PGS", data.frame( anova( pmv3,pmv5,test = "LRT" ) )[2,] ),
      c( "FH_FGRS", data.frame( anova( pmv3,pmv6,test = "LRT" ) )[2,] ),
      c( "PGS_FGRS_FH", data.frame( anova( pmv4,pmv7,test = "LRT" ) )[2,] ),
      c( "PGS_FH_FGRS", data.frame( anova( pmv5,pmv7,test = "LRT" ) )[2,] ),
      c( "FGRS_FH_PGS", data.frame( anova( pmv6,pmv7,test = "LRT" ) )[2,] )
    ) )
    stat.test$p <- pchisq( as.numeric(stat.test$Deviance), as.numeric(stat.test$Df),
                           lower.tail = F )
    #				, log.p = T )/log(10)					## Why?
    #	stat.test$e <- floor( stat.test$p )						## Why?
    #	stat.test$p <- 10^(stat.test$p -stat.test$e)					## Why?
    
    stat.test[,3] <- round( as.numeric( stat.test[,3] ),2 )
    stat.test[,5] <- round( as.numeric( stat.test[,5] ),2 )
    stat.test[,6] <- signif( as.numeric( stat.test[,6] ),3 )
    stat.test[,7] <- signif( as.numeric( stat.test[,7] ),3 )
    
    n <- length( data$dx )
    nCase <- sum( data$dx == 1 )
    nControl <- sum( data$dx == 0 )
    
    cor.wt.1 <- cov.wt( cbind( data$PGS_eu, data$FGRS ), wt=wt, cor=T )$cor[ 1,2 ]
    cor.wt.2 <- cov.wt( cbind( data$PGS_eu, data$FH ), wt=wt, cor=T )$cor[ 1,2 ]
    cor.wt.3 <- cov.wt( cbind( data$FGRS, data$FH ), wt=wt, cor=T )$cor[ 1,2 ]

    fgrs.res <- lm( data=data, scale( FGRS )~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 
                 	, weights=wt )$residual
    fh.res <- lm( data=data, scale( FH )~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 
                 	, weights=wt )$residual
    pgs.res <- lm( data=data, scale( PGS_eu )~age+gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 
                 	, weights=wt )$residual

    pcor.wt.1 <- cov.wt( cbind( fgrs.res, pgs.res ), wt=wt, cor=T )$cor[ 1,2 ]
    pcor.wt.2 <- cov.wt( cbind( fh.res, pgs.res ), wt=wt, cor=T )$cor[ 1,2 ]
    pcor.wt.3 <- cov.wt( cbind( fgrs.res, fh.res ), wt=wt, cor=T )$cor[ 1,2 ]

    results <- rbind( results, c( c, p, n, nCase, nControl, 
                                  R2_prob0, R2_prob1, R2_prob2, R2_prob3, 
				  R2_prob4, R2_prob5, R2_prob6, R2_prob7,
                                  stat.test$p[ 2:13 ], 
				  cor.wt.1, cor.wt.2, cor.wt.3,
				  pcor.wt.1, pcor.wt.2, pcor.wt.3 ) )
    
    fwrite( stat.test, 
            file=paste0( "~/MDD_heterogeneity/backup/data/", c, "_", 
                         p, "_modelStats_weightedProbit_withCovs_unrel_weights_230127.txt"  ), 
            row.names=F, col.names=T, sep='\t' )
    
  }
}

results <- data.table( results )

for ( i in 6:25 ) {
  results[[i]] <- as.numeric( results[[i]] )
}

names( results ) <- c( 'cohort', 'trait', 'n', 'n_case', 'n_control', 
                       'r2l_covs', 'r2l_PGS', 'r2l_FGRS', 'r2l_FH', 
                       'r2l_PGS_FGRS', 'r2l_PGS_FH', 'r2l_FGRS_FH', 
                       'r2l_PGS_FGRS_FH',
                       paste0( "p_", stat.test$V1[2:13] ), 
		       'PGS_FGRS_cor_wt', 'PGS_FH_cor_wt', 'FGRS_FH_cor_wt', 
		       'PGS_FGRS_pcor_wt', 'PGS_FH_pcor_wt', 'FGRS_FH_pcor_wt' )

write.table( results, "~/MDD_heterogeneity/backup/data/R2Liab_weightedProbit_withCovs_unrel_weights_230127.txt", 
             row.names=F, col.names=T, quote=F )

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

p.t <- 0.05/(6*5*2)
ps <- results
for ( i in 14:25 ) { ps[[ i ]] <- 1*( ps[[ i ]] > p.t ) }


pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/r2l_ipsych2012_withCovs_unrel_weights_230127.pdf", width=5, height=5 )
par( mgp=c(2,0.5,0), tcl=-0.25 )

r2l.covs <- t( as.matrix( cbind( results[ c(2,4:7), 6 ], results[ c(2,4:7), 6 ], 
                                 results[ c(2,4:7), 6 ], results[ c(2,4:7), 6 ], 
                                 results[ c(2,4:7), 6 ], results[ c(2,4:7), 6 ] ) ) )
r2l <- t( as.matrix( results[ c(2,4:7), c(7:11,13) ] ) )
r2l <- r2l - r2l.covs
p.r2l <- ps[ c(2,4:7), c(2,14:19,21,23:24) ]

mids <- barplot( r2l, beside=T, col=gg_color_hue(6), names.arg=c( "ADHD","ASD","BPD","MDD","SCZ" ),
                 ylim=c( 0,0.075 ), ylab=expression( r[Liab.]^2 ) ) 
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


pdf( "/faststorage/jail/project/MDD_heterogeneity/backup/data/r2l_ipsych2015i_withCovs_unrel_weights_230127.pdf", width=5, height=5 )
par( mgp=c(2,0.5,0), tcl=-0.25 )

r2l.covs <- t( as.matrix( cbind( results[ c(9,11:14), 6 ], results[ c(9,11:14), 6 ], 
                                 results[ c(9,11:14), 6 ], results[ c(9,11:14), 6 ], 
                                 results[ c(9,11:14), 6 ], results[ c(9,11:14), 6 ] ) ) )
r2l <- t( as.matrix( results[ c(9,11:14), c(7:11,13) ] ) )
r2l <- r2l - r2l.covs
p.r2l <- ps[ c(9,11:14), c(2,14:19,21,23:24) ]

mids <- barplot( r2l, beside=T, col=gg_color_hue(6), names.arg=c( "ADHD","ASD","BPD","MDD","SCZ" ),
                 ylim=c( 0,0.075 ), ylab=expression( r[Liab.]^2 ) ) 
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

