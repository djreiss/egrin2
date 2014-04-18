######################################################################
######################################################################
######################################################################
#################################################################
## Compute per-run fraction of times each motif is in a coding region, save the results

p.cutoff <- 1e-6

do.coding.fracs <- function(f1) {
    cat('CODING FRACS', f1, '\n')
    f1a <- gsub( '.RData', '', basename(f1) )
    out.file <- sprintf("./%s/%s_coding_fracs.RData", output.dir, f1a)
    fimo.file <- sprintf("./%s/%s_fimo_out.bz2", output.dir, f1a)
    if ( file.exists( out.file ) ) { load( out.file ); return( coding.fracs ) }
    if ( ! file.exists( fimo.file ) ) { print( "SKIPPING (2)" ); return() }
    load( sprintf( '%s/%s', output.dir, f1 ) )
    fimo.out <- read.in.fimo( fimo.file, p.cutoff ) ##read.delim( bzfile(fimo.file) )
    
    coding.seqs <- lapply( names( e$genome.info$genome.seqs ), function( n )
                          rep( FALSE, nchar( e$genome.info$genome.seqs[ n ] ) ) )
    names( coding.seqs ) <- names( e$genome.info$genome.seqs )
    tmp <- subset( e$genome.info$feature.tab, type %in% c( 'CDS', 'rRNA', 'tRNA' ) )
    tmp$where <- as.character( tmp$contig )
    tmp$Start <- as.integer( as.character( tmp$start_pos ) )
    tmp$Stop <- as.integer( as.character( tmp$end_pos ) )
    for ( i in 1:nrow( tmp ) ) {
        ttmp <- tmp[ i, ]
        if ( ! is.na( ttmp$where ) && ! is.na( ttmp$Start ) && ! is.na( ttmp$Stop ) )
            coding.seqs[[ ttmp$where ]][ ttmp$Start:ttmp$Stop ] <- TRUE
        rm( ttmp )
    }
  
    scans <- as.data.table( subset( fimo.out, `p-value` <= p.cutoff ) )

    in.coding <- rep( NA, nrow( scans ) )
    for ( cc in names( coding.seqs ) ) {
        in.coding[ scans$Seq == cc ] <- coding.seqs[[ cc ]][ scans[ scans$Seq == cc, ]$Start ] |
            coding.seqs[[ cc ]][ scans[ scans$Seq == cc, ]$Stop ]
    }
    scans$in.coding <- in.coding

    mots <- unique( paste( scans$bic, scans$mot, sep='_' ) )

    frac.in.coding <- do.call( c, lapply( 1:length( mots ), function( i ) {
        bimo <- as.integer( strsplit( mots[i], '_' )[[ 1 ]] )
        sc <- scans[ bic == bimo[1] & abs( mot ) == bimo[2] ]
        if ( nrow( sc ) <= 0 || all( is.na( sc$Start ) ) ) return( NA ) ##next
        hits <- sc$in.coding
        mean( hits, na.rm=T )
    } ) )

    names( frac.in.coding ) <- mots
    mean.coding.frac <- sum( sapply( coding.seqs, sum, na.rm=T ) ) / sum( sapply( coding.seqs, length ) )
    coding.fracs <- list( all.fracs=frac.in.coding, mean.fracs=mean.coding.frac, p.cutoff=p.cutoff )
    save( coding.fracs, file=out.file )
    coding.fracs
}    

##for ( f1 in 1:length(rdatas) ) {
coding.fracs <- mclapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
    do.coding.fracs( f1 )
}, mc.preschedule=F )

names(coding.fracs) <- rdatas
