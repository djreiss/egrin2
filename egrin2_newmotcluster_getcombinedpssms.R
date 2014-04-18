######################################################################
######################################################################
######################################################################
#########################################################################
## Compute combined pssms for each motif cluster. Here we need to run tomtom
##  between all of the motifs in each motif cluster.
## It will also plot them if an individual cmonkey run "e" is in the environment

cf.names <- sapply( strsplit( gsub( '.RData', '', basename( rdatas ) ), '_' ), '[', 3 )

e.value.cutoff=100; resid.cutoff=0.8; dist.meth="ed"; q.thresh=0.5; min.overlap=4; q.pseudo=0; t.pseudo=0
cmd <- "./progs/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"

if ( ! exists( 'all.rdata.meme.outs' ) ) all.rdata.meme.outs <- get.all.meme.outs()

get.aligned.pssm <- function( tt.out, seq.type='upstream meme', n.gene.weight=F, return.aligned.pssms=F ) {
    meme.let <- c( "A", "C", "G", "T" )

    mot.names <- unique( c( as.character( tt.out$`Query ID` ), as.character( tt.out$`Target ID` ) ) )
    ## unfortunate use of '_' to separate bic_mot id's (and they exist in filenames) so some trickery here:
    ## Should be last 4 values are bic, mot, resid, e-val, so everything else if RData filename:
    print( mot.names[ 1 ] )
    tmp1 <- strsplit( gsub( 'bic_', '', mot.names[ 1 ] ), "_" )[[ 1 ]]
    rd <- paste( tmp1[ 1:(length(tmp1)-4) ], collapse='_' )
    ##tmp.env <- new.env(); load( sprintf( '%s/%s', output.dir, rd ), envir=tmp.env )

    ## Get first motif, then take all alignments to it, and construct combined pssm
    tmp1 <- as.integer( tmp1 )
    tmp.mot <- all.rdata.meme.outs[[ rd ]][[ tmp1[ length(tmp1)-3 ] ]]$meme.out[[ tmp1[ length(tmp1)-2 ] ]]
    ##tmp.mot <- tmp.env$e$meme.scores[[ seq.type ]][[ tmp1[ length(tmp1)-3 ] ]]$meme.out[[ tmp1[ length(tmp1)-2 ] ]]
    ##rm( tmp.env )
    
    pssm <- orig.pssm <- tmp.mot$pssm
    colnames( pssm ) <- meme.let
    if ( n.gene.weight ) pssm <- pssm * tmp.mot$sites
    if ( return.aligned.pssms ) aligned.pssms <- list()

    max.width <- max( nchar( c( as.character( tt.out$`Query consensus` ),
                               as.character( tt.out$`Target consensus` ) ) ) )
    for ( jj in 1:max.width ) pssm <- rbind( rep( 0, 4 ), pssm, rep( 0, 4 ) )
    if ( return.aligned.pssms ) aligned.pssms[[ mot.names[ 1 ] ]] <- pssm
    first.ind <- which( apply( pssm, 1, function( j ) any( j != 0 ) ) )[ 1 ]
    orig.width <- nrow( orig.pssm )
    tttt <- subset( tt.out, `Query ID` == mot.names[ 1 ] | `Target ID` == mot.names[ 1 ] )

    for ( m in mot.names[ 2:length( mot.names ) ] ) {
        ttt <- unique( subset( tttt, `Query ID` == m | `Target ID` == m ) )
        if ( nrow( ttt ) <= 0 ) next
        else if ( nrow( ttt ) > 1 ) ttt <- ttt[ 1, ]

        print( m )
        tmp1 <- strsplit( gsub( 'bic_', '', m ), "_" )[[ 1 ]]
        rd <- paste( tmp1[ 1:(length(tmp1)-4) ], collapse='_' )
        ##tmp.env <- new.env(); load( sprintf( '%s/%s', output.dir, rd ), envir=tmp.env )
        tmp1 <- as.integer( tmp1 )
        tmp.mot <- all.rdata.meme.outs[[ rd ]][[ tmp1[ length(tmp1)-3 ] ]]$meme.out[[ tmp1[ length(tmp1)-2 ] ]]
        ##tmp.mot <- tmp.env$e$meme.scores[[ seq.type ]][[ tmp1[ length(tmp1)-3 ] ]]$meme.out[[ tmp1[ length(tmp1)-2 ] ]]
        ##rm( tmp.env )
        
        pssm2 <- tmp.mot$pssm
        if ( n.gene.weight ) pssm2 <- pssm2 * tmp.mot$sites
        if ( ttt$Orientation == "-" ) {
            if ( ttt$`Target ID` == mot.names[1] ) offset <- first.ind + orig.width - ttt$`Optimal offset` - nrow( pssm2 )
            else offset <- first.ind - ttt$`Optimal offset`
            pssm2 <- pssm2[ ,4:1 ][ nrow( pssm2 ):1, ]
        } else {
            if ( ttt$`Query ID` == mot.names[1] ) offset <- first.ind - ttt$`Optimal offset`
            else offset <- first.ind + ttt$`Optimal offset`
        }
        pssm[ offset:( offset + nrow( pssm2 ) - 1 ), ] <- pssm[ offset:( offset + nrow( pssm2 ) - 1 ), ] + pssm2
        if ( return.aligned.pssms ) {
            tmp.pssm <- pssm * 0
            tmp.pssm[ offset:( offset + nrow( pssm2 ) - 1 ), ] <- pssm2
            aligned.pssms[[ m ]] <- tmp.pssm
        }
    }

    first.ind2 <- which( apply( pssm, 1, function( j ) any( j != 0 ) ) )[ 1 ]
    pssm <- pssm[ -( 1:( first.ind2 - 1 ) ), ]
    last.ind2 <- which( apply( pssm, 1, function( j ) all( j == 0 ) ) )[ 1 ]
    if ( ! is.na( last.ind2 ) ) pssm <- pssm[ -( ( last.ind2 - 1 ):nrow( pssm ) ), ]
    pssm <- ( pssm + 1e-9 ) / ( max( pssm, na.rm=T ) + 1e-9 )
    combined.pssm <- pssm

    if ( return.aligned.pssms ) {
        for ( m in names( aligned.pssms ) ) {
            pssm <- aligned.pssms[[ m ]]
            pssm <- pssm[ -( 1:( first.ind2 - 1 ) ), ]
            pssm <- pssm[ -( ( last.ind2 - 1 ):nrow( pssm ) ), ]
            aligned.pssms[[ m ]] <- pssm
        }
        return( list( combined=combined.pssm, aligned=aligned.pssms ) )
    }
    combined.pssm
}

if ( ! exists( 'all.rdata.meme.files' ) ) {
    all.rdata.meme.files <- get.all.meme.files()
}

##for ( i in 1:length( clusts ) ) {
new.clusts <- lapply( clusts[ 1:attr(clusts,'mc.length')], function( clust ) {
    ##clust <- clusts[[ i ]]
    print( length( clust ) )
    if ( ! is.null( attr( clust, 'aligned.pssm' ) ) ) return( clust )
    tmp <- do.call( rbind, strsplit( clust, '_' ) )
    tmp <- cbind( tmp[,1], paste( tmp[,2], tmp[,3], sep='_' ) )
    tmp <- tapply( tmp[,2], tmp[,1], c )
    mots.files <- lapply( names(tmp), function( t ) {
        #print( t )
        rd <- rdatas[ which( cf.names == t ) ]
        #load( sprintf( '%s/%s', output.dir, rd ) )
        #e$all.motifs.to.meme.file <- all.motifs.to.meme.file; environment(e$all.motifs.to.meme.file) <- e
        tt <- strsplit( tmp[[t]], '_' )
        ks <- as.integer(sapply(tt,'[', 1))
        ms <- as.integer(sapply(tt,'[', 2))
        #mots.file <- e$all.motifs.to.meme.file(ks=ks, ms=ms, to.file=F, e.val=Inf, resid=Inf)
        #mots.file <- gsub( 'bic_', sprintf( 'bic_%s_', basename(rd) ), mots.file )
        mots.files <- lapply( 1:length(ks), function(i) all.rdata.meme.files[[ rd ]][[ ks[i] ]][[ ms[i] ]] )
        mots.file <- c( mots.files[[1]], unlist( lapply( mots.files[2:length(mots.files)], '[', -(1:8) ) ) )
        mots.file
    } )

    mots.file <- c( mots.files[[1]], unlist( lapply( mots.files[2:length(mots.files)], '[', -(1:8) ) ) )
    tfile <- tempfile( "tomtom_t_", )
    cat( mots.file, file=tfile, sep="\n" ) ## Write out all motifs
    cmd.1 <- sprintf( cmd, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tfile )
    cmd.1 <- paste( cmd.1, "-query", tfile )
    print( cmd.1 )
    tout.orig <- system( cmd.1, intern=T ) ## note this has not been filtered to remove self-self comparisons
    tout <- do.call( rbind, strsplit( tout.orig, "\t" ) )
    colnames( tout ) <- gsub( '#', '', tout[ 1, ,drop=F ] ); tout <- tout[ -1, ,drop=F ]
    tout <- as.data.frame( tout[ tout[ ,1 ] != tout[ ,2 ], ,drop=F ] ) ## de-symmetrize the output
    tout$`p-value` <- as.numeric( as.character( tout$`p-value` ) )
    tout$`q-value` <- as.numeric( as.character( tout$`q-value` ) )
    tout$`Optimal offset` <- as.integer( as.character( tout$`Optimal offset` ) )
    tout <- tout[ order( tout$`q-value`, tout$`p-value` ), ]
    cat( "TOUT:", dim( tout ), "\n" )

    tt.out <- subset( tout, `p-value` <= 0.01 )
    attr( clust, 'tt.out' ) <- tmp$tt.out
    pssm <- get.aligned.pssm( tt.out )
    attr( clust, 'aligned.pssm' ) <- pssm
    if ( exists( 'e' ) ) e$viewPssm( pssm, main=paste( length(clust),
                         length(unique(c(as.character(tt.out$`Query ID`),as.character(tt.out$`Target ID`)))) ) )
    clust
} ) ##, mc.preschedule=F )
