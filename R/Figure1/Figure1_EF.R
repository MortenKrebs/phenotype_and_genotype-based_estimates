rm( list = ls() )

## Packages

library( data.table )
library( corrplot )
library( knitr )

## Functions

gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

## Load Data

load( file='/home/anscho/MDD_heterogeneity/backup/data/corMats_221031.RData' )
load( file='/home/anscho/MDD_heterogeneity/backup/data/psychData_AJS_221031.RData' )

## Pearson - PGS / FGRS 
t <- cor.p.fgrs.i2012[ c(3,6:8,10), c(1,3,5,7,9) ]
row.names( t ) <- colnames( t ) <- c( "ADHD", "ASD", "BPD", "MDD", "SCZ" )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_PGS_FGRS_2012_221101_eurUnrel.pdf", width=5, height=5 )

	corrplot( t,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.pos='n',
		mar=c(0,4,4,0),
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 1, col.lim = c(-.1,1),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 2, at=5:1, col=gg_color_hue(6)[1], las=2 )
	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 3, at=1:5, line=-1, col=gg_color_hue(6)[2], las=2 )
	mtext( c( 'PGS' ), 2, at=3, line=2, col=gg_color_hue(6)[1] )
	mtext( c( 'FGRS' ), 3, at=3, line=2, col=gg_color_hue(6)[2] )

dev.off()

t <- cor.p.fgrs.i2015i[ c(3,6:8,10), c(1,3,5,7,9) ]
row.names( t ) <- colnames( t ) <- c( "ADHD", "ASD", "BPD", "MDD", "SCZ" )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_PGS_FGRS_2015i_221101_eurUnrel.pdf", width=5, height=5 )

	corrplot( t,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.pos='n',
		mar=c(0,4,4,0),
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 1, col.lim = c(-.1,1),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 2, at=5:1, col=gg_color_hue(6)[1], las=2 )
#	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 3, at=1:5, line=-1, col=gg_color_hue(6)[2] )
	mtext( c( 'PGS' ), 2, at=3, line=2, col=gg_color_hue(6)[1] )
#	mtext( c( 'FGRS' ), 3, at=3, line=0, col=gg_color_hue(6)[2] )

dev.off()

## Full Data Correlations

cor.p.all.2012 <- round( cor( corDat.2012[ ,c( 4,7:9,11,12,14,16,18,20,22:26 ) ] ), 3 )
row.names( cor.p.all.2012 ) <- colnames( cor.p.all.2012 ) <- rep( c( "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 3 )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_ALL_2012_221101_eurUnrel.pdf", width=8, height=8 )

	corrplot( cor.p.all.2012,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.col=c( rep( gg_color_hue(6)[1], 5 ), rep( gg_color_hue(6)[2], 5 ), rep( gg_color_hue(6)[3], 5 ) ),
		tl.cex=1,
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 0.75, col.lim = c(-.1,1),
		mar=c( 0,4,1,0),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( "PGS", 2, line=1, at=12.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 2, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 2, line=1, at=2.5, col=gg_color_hue(6)[3] ) 
	mtext( "PGS", 3, line=1, at=2.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 3, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 3, line=1, at=12.5, col=gg_color_hue(6)[3] ) 

dev.off()


cor.p.all.2015i <- round( cor( corDat.2015i[ ,c( 4,7:9,11,12,14,16,18,20,22:26 ) ] ), 3 )
row.names( cor.p.all.2015i ) <- colnames( cor.p.all.2015i ) <- rep( c( "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 3 )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_ALL_2015i_221101_eurUnrel.pdf", width=8, height=8 )

	corrplot( cor.p.all.2015i,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.col=c( rep( gg_color_hue(6)[1], 5 ), rep( gg_color_hue(6)[2], 5 ), rep( gg_color_hue(6)[3], 5 ) ),
		tl.cex=1,
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 0.75, col.lim = c(-.1,1),
		mar=c( 0,4,1,0),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( "PGS", 2, line=1, at=12.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 2, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 2, line=1, at=2.5, col=gg_color_hue(6)[3] ) 
	mtext( "PGS", 3, line=1, at=2.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 3, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 3, line=1, at=12.5, col=gg_color_hue(6)[3] ) 

dev.off()

## European with Relatives


## Pearson - PGS / FGRS 
t <- cor.p.fgrs.i2012[ c(2,5,7,8,10), c(1,3,5,7,9) ]
row.names( t ) <- colnames( t ) <- c( "ADHD", "ASD", "BPD", "MDD", "SCZ" )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_PGS_FGRS_2012_221101_eur.pdf", width=5, height=5 )

	corrplot( t,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.pos='n',
		mar=c(0,4,4,0),
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 1, col.lim = c(-.1,1),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 2, at=5:1, col=gg_color_hue(6)[1], las=2 )
	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 3, at=1:5, line=-1, col=gg_color_hue(6)[2], las=2 )
	mtext( c( 'PGS' ), 2, at=3, line=2, col=gg_color_hue(6)[1] )
	mtext( c( 'FGRS' ), 3, at=3, line=2, col=gg_color_hue(6)[2] )

dev.off()

t <- cor.p.fgrs.i2015i[ c(2,5,7,8,10), c(1,3,5,7,9) ]
row.names( t ) <- colnames( t ) <- c( "ADHD", "ASD", "BPD", "MDD", "SCZ" )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_PGS_FGRS_2015i_221101_eur.pdf", width=5, height=5 )

	corrplot( t,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.pos='n',
		mar=c(0,4,4,0),
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 1, col.lim = c(-.1,1),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 2, at=5:1, col=gg_color_hue(6)[1], las=2 )
#	mtext( c( 'ADHD','ASD','BPD','MDD','SCZ' ), 3, at=1:5, line=-1, col=gg_color_hue(6)[2] )
	mtext( c( 'PGS' ), 2, at=3, line=2, col=gg_color_hue(6)[1] )
#	mtext( c( 'FGRS' ), 3, at=3, line=0, col=gg_color_hue(6)[2] )

dev.off()

## Full Data Correlations

cor.p.all.2012 <- round( cor( corDat.2012[ ,c( 3,6,8:9,11,12,14,16,18,20,22:26 ) ] ), 3 )
row.names( cor.p.all.2012 ) <- colnames( cor.p.all.2012 ) <- rep( c( "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 3 )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_ALL_2012_221101_eur.pdf", width=8, height=8 )

	corrplot( cor.p.all.2012,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.col=c( rep( gg_color_hue(6)[1], 5 ), rep( gg_color_hue(6)[2], 5 ), rep( gg_color_hue(6)[3], 5 ) ),
		tl.cex=1,
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 0.75, col.lim = c(-.1,1),
		mar=c( 0,4,1,0),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( "PGS", 2, line=1, at=12.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 2, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 2, line=1, at=2.5, col=gg_color_hue(6)[3] ) 
	mtext( "PGS", 3, line=1, at=2.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 3, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 3, line=1, at=12.5, col=gg_color_hue(6)[3] ) 

dev.off()


cor.p.all.2015i <- round( cor( corDat.2015i[ ,c( 3,6,8:9,11,12,14,16,18,20,22:26 ) ] ), 3 )
row.names( cor.p.all.2015i ) <- colnames( cor.p.all.2015i ) <- rep( c( "ADHD", "ASD", "BPD", "MDD", "SCZ" ), 3 )

pdf( "/home/anscho/MDD_heterogeneity/backup/data/CorPlot_ALL_2015i_221101_eur.pdf", width=8, height=8 )

	corrplot( cor.p.all.2015i,
		method="color",
		#p.mat = res1$p,
		#insig = "label_sig",
		#sig.level = c(.001, .01, .05),
		pch.cex = 0.8,
		pch.col = "red",
		tl.col=c( rep( gg_color_hue(6)[1], 5 ), rep( gg_color_hue(6)[2], 5 ), rep( gg_color_hue(6)[3], 5 ) ),
		tl.cex=1,
		outline=TRUE, addCoef.col = "black", number.digits = 3, number.cex = 0.75, col.lim = c(-.1,1),
		mar=c( 0,4,1,0),
		col=colorRampPalette(c("darkred","red","white", "deepskyblue","darkblue"))(100),
	)
	mtext( "PGS", 2, line=1, at=12.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 2, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 2, line=1, at=2.5, col=gg_color_hue(6)[3] ) 
	mtext( "PGS", 3, line=1, at=2.5, col=gg_color_hue(6)[1] ) 
	mtext( "FGRS", 3, line=1, at=7.5, col=gg_color_hue(6)[2] ) 
	mtext( "FH", 3, line=1, at=12.5, col=gg_color_hue(6)[3] ) 

dev.off()











