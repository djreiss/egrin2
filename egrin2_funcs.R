all.motifs.to.meme.file <- function( ks=1:k.clust, ms=NA, seq.type=names(mot.weights)[1],
                                    e.value.cutoff=100, resid.cutoff=0.8, to.file=TRUE ) {
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )

  cluster.motif.lines <- function( k, m=NA ) { ## Write out a single bicluster's (usually 2) motifs to MEME format
      if ( k %% 100 == 0 ) print(k)
      lines <- character()
      memeOut <- meme.scores[[ seq.type ]][[ k ]]
      if ( is.null( memeOut ) || memeOut == "" ) return( lines )
      memeOut <- meme.scores[[ seq.type ]][[ k ]]$meme.out
      if ( is.null( memeOut ) ) return( lines )
      if ( ! is.infinite( resid.cutoff ) && clusterStack[[ k ]]$resid > resid.cutoff ) return( lines )
      ##max.motifs <- max( max.motifs, length( memeOut ) )
      for ( i in 1:length( memeOut ) ) {
          if ( ! is.na( m ) && i != m ) next
          if ( memeOut[[ i ]]$e.value > e.value.cutoff ) next
          pssm <- memeOut[[ i ]]$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
          lines <- c( lines, "",
                 sprintf( "MOTIF bic_%03d_%02d_%.3f_%.3e", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                     sprintf( "BL   MOTIF bic_%03d_%02d_%.3f_%.3e width=0 seqs=0", k, i,
                             clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                     sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
                             memeOut[[ i ]]$sites, memeOut[[ i ]]$e.value ) )
          lines <- c( lines, apply( pssm, 1, function( i )
                                   sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
      }
      lines
  }

  ##lines.t <- lapply( ks, cluster.motif.lines )
  if ( is.null( ms ) || is.na( ms ) ) ms <- rep( NA, length(ks) )
  lines.t <- lapply( 1:length(ks), function(i) cluster.motif.lines( ks[i], ms[i] ) )
  if ( length( ks ) > 1 ) lines.t <- c( lines, do.call( 'c', lines.t ) )
  else lines.t <- c( lines, lines.t[[1]] )
  if ( to.file ) {
      tfile <- my.tempfile( "tomtom_t_", )
      cat( lines.t, file=tfile, sep="\n" ) ## Write out all motifs
      return( tfile )
  }
  return( lines.t )
}

read.in.fimo <- function( fimo.file, p.cutoff ) {
    fimo.out <- system( sprintf( 'bunzip2 -c %s | awk \'($6<=%s){print}\'',
                                fimo.file, as.character( p.cutoff ) ), intern=T )
    fimo.out <- as.data.frame( do.call( rbind, lapply( fimo.out, function( i ) strsplit( i, '\t' )[[ 1 ]] ) ) )
    colnames( fimo.out ) <- strsplit( system( sprintf( 'bunzip2 -c %s | head -1', fimo.file ), intern=T ), '\t' )[[ 1 ]]
    fimo.out$Start <- as.integer( as.character( fimo.out$Start ) )
    fimo.out$Stop <- as.integer( as.character( fimo.out$Stop ) )
    fimo.out$`Log-odds` <- as.numeric( as.character( fimo.out$`Log-odds` ) )
    fimo.out$`p-value` <- as.numeric( as.character( fimo.out$`p-value` ) )
    tmp <- do.call( rbind, strsplit( as.character( fimo.out$Motif ), '_' ) )
    fimo.out$bic <- as.integer( tmp[ ,4 ] )
    fimo.out$mot <- as.integer( tmp[ ,5 ] )
    fimo.out$Strand <- substr( tmp[ ,1 ], 1, 1 )
    rm( tmp )
    fimo.out$Motif <- NULL
    fimo.out <- as.data.table( fimo.out )
    setkey( fimo.out, bic, mot, Seq, Start )
    fimo.out
}

get.clusterStack <- function( rd ) {
    cat( 'CLUSTERSTACKS', rd, '\n' )
    filename <- sprintf( '%s/%s_clusterStacks.RData', output.dir, gsub('.RData', '', rd) )
    if ( file.exists( filename ) ) { load( filename ); return( clusterStack ) }
        
    load( sprintf( '%s/%s', output.dir, rd ) )
    clusterStack <- e$clusterStack
    save( clusterStack, file=filename )
    clusterStack
}

get.all.clusterStacks <- function() {
    all.rdata.clusterStacks <- lapply( rdatas, function( rd ) {
        get.clusterStack( rd )
    } )
    names( all.rdata.clusterStacks ) <- rdatas
    all.rdata.clusterStacks
}

## Get meme.out structure each motif in each bicluster in each cmonkey run (contains pssms)
## save individually
get.meme.outs <- function( rd, seq.type='upstream meme' ) {
    cat( 'MEME OUTS', rd, '\n' )
    filename <- sprintf( '%s/%s_meme_outs.RData', output.dir, gsub('.RData', '', rd) )
    if ( file.exists( filename ) ) { load( filename ); return( meme.outs ) }
    
    load( sprintf( '%s/%s', output.dir, rd ) )
    meme.outs <- lapply( e$meme.scores[[ seq.type ]], function( out ) {
        out$pv.ev <- NULL; out$prev.run <- NULL; out } )
    meme.outs$all.pv <- NULL
    ##print(object.size(meme.outs))
    save( meme.outs, file=filename )
    meme.outs
}

get.all.meme.outs <- function( seq.type='upstream meme' ) {
    all.rdata.meme.outs <- lapply( rdatas, function( rd ) {
        get.meme.outs( rd, seq.type )
    } )
    names( all.rdata.meme.outs ) <- rdatas
    all.rdata.meme.outs
}

get.meme.files <- function( rd ) {
    cat( 'MEME FILES', rd, '\n' )
    filename <- sprintf( '%s/%s_meme_files.RData', output.dir, gsub('.RData', '', rd) )
    if ( file.exists( filename ) ) { load( filename ); return( mots.files ) }
    
    load( sprintf( '%s/%s', output.dir, rd ) )
    e$all.motifs.to.meme.file <- all.motifs.to.meme.file; environment(e$all.motifs.to.meme.file) <- e
    mots.files <- lapply( 1:e$k.clust, function(k)
                         lapply( 1:max(unlist(e$n.motifs)), function(m) {
                             mf <- e$all.motifs.to.meme.file(ks=k, ms=m, to.file=F, e.val=Inf, resid=Inf)
                             mf <- gsub( 'bic_', sprintf( 'bic_%s_', basename(rd) ), mf )
                         } ) )
    save( mots.files, file=filename )
    mots.files
}

get.all.meme.files <- function() {
    all.rdata.meme.files <- mclapply( rdatas, function( rd ) {
        get.meme.files( rd )
    } )
    names( all.rdata.meme.files ) <- rdatas
    all.rdata.meme.files
}
