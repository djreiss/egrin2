## NEW STRATEGY FOR DREAM5 DATA:
##   NOTE THIS IS DIFFERENT FROM zzz_eco4.R in that we DO vary parameters; we subset the data even further (randomly)
##       and we MUST USE CMONKEY VERSION 4.9.0!!!

## use raw 5-fold cv -- remove the 1st 161 conditions for iter 1,11,21, etc; remove
## conditions 162-323 for iters 2,12,22, etc...


if ( length( grep( "--iter=", commandArgs(), fixed=T ) > 0 ) ) {
  iter <- commandArgs(); iter <- grep( "--iter=", iter, fixed=T, val=T )
  iter <- strsplit( iter, "=", fixed=T )[[ 1 ]][ 2 ]
}

if ( exists( "iter" ) ) {
  debug.on()
  require( cMonkey )

  if ( TRUE ) { ## Load the data. Using data from DREAM5 competition!
    ##   http://homes.esat.kuleuven.be/~kmarchal/Supplementary_Information_Lemmens_2008/Index.html#conditional
    load('format_data.RData')
    rm(list=ls()[!ls()%in%c('ratios.normed','iter')])
    ratios <- ratios.normed; rm(ratios.normed)
    colnames(ratios) <- paste( 'C', 1:ncol(ratios), sep='' )

    remove.seq <- round(seq(1,805,length=6))
    names(remove.seq) <- as.character(0:5)
    iter <- as.integer(iter)
    print(iter)
    conditions.excluded <- remove.seq[as.character(iter%%5)]:(remove.seq[as.character(iter%%5+1)])
    print( conditions.excluded )
    ratios <- ratios[ ,-(conditions.excluded) ]
  }

  rnd.seed <- iter ##as.integer( substr( gsub( '[-:. ]', "", as.character( Sys.time() ) ), 12, 20 ) )
  set.seed( rnd.seed )
  
  n.iter <- 2000
  parallel.cores <- 1; parallel.cores.motif <- 2
  require( parallel )
  options( cores=2 ); options( mc.cores=2 )
  ratios <- ratios[ ,sample( 1:ncol( ratios ), sample( 180:280, 1 ) ) ]
  k.clust <- sample( 350:550, 1 )
  motif.upstream.scan <- c( sample( (-50):0, 1 ), sample( 150:250, 1 ) )
  motif.upstream.search <- c( sample( (-20):0, 1 ), sample( 100:200, 1 ) )
  net.weights <- c( string=runif( 1 ) * 0.5 + 0.2, operons=runif( 1 ) * 0.5 + 0.2 )
  if ( sample( 1:2, 1 ) == 1 ) net.weights <- net.weights[ -which( names( net.weights ) == 'string' ) ]
  if ( sample( 1:2, 1 ) == 1 ) net.weights <- net.weights[ -which( names( net.weights ) == 'operons' ) ]
  maxw <- sample( 12:30, 1 )
  meme.cmd <- "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 24 -mod zoops -nostatus -text -cons $none"
  meme.cmd <- gsub( "maxw 24", sprintf( "maxw %d", maxw ), meme.cmd )
  n.motifs <- sample( 1:3, 1 )
  row.scaling <- sample( 4:8, 1 ) ## default 6

  e <- cmonkey.init( organism="eco", ratios=ratios, plot.iters=0 ) ##k.clust=k.clust, net.weights=net.weights, 
  ##n.clust.per.row=5, n.clust.per.col=250 )

  cmonkey( e, dont.init=T )
  save.image( file=sprintf( "output/zzz_eco_%03d.RData", as.integer( iter ) ) )
  cat( 'DONE:', iter, date(), '\n' )
  
} else {
    
  require( parallel )
  options( cores=4 ); options( mc.cores=4 )
  if ( Sys.getenv('HOST') == 'massive' ) { options( cores=8 ); options( mc.cores=8 ) }
  if ( ! exists( 'range.to.run' ) ) range.to.run <- 1:1000
  mclapply( range.to.run, function(iter) {
    cat( iter, "\n" )
    if ( file.exists( sprintf( "output/zzz_eco_%03d.Rout", iter ) ) ) next
    cmd <- sprintf( "/bin/nice -9 R CMD BATCH --no-save --no-restore --iter=%d zzz_eco5.R output/zzz_eco_%03d.Rout",
                   iter, iter )
    print( cmd )
    print( date() )
    system( cmd )
  }, mc.preschedule=F )
  
}

debug.off()
stop()

