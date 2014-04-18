######################################################################
######################################################################
######################################################################
#########################################################################
## Compute the "motif shadows" for each run -- locations covered by each motif in each run -
##   save as motif_shadows.RData

p.cutoff <- 1e-6

do.motif.shadows <- function(f1) {
    cat('MOTIF SHADOWS', f1, '\n')
    f1a <- gsub( '.RData', '', f1 )
    if ( file.exists( sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) ) ) {
        load( sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) )
        return( m )
    }
    system( sprintf('touch %s/%s_motif_shadows.RData', output.dir, f1a) ) ## placeholder
    load( sprintf( '%s/%s', output.dir, f1 ) )
    nc <- sapply( e$genome.info$genome.seqs, nchar )
    rm( e, ratios ); gc()

    f1a <- gsub( '.RData', '', f1 )
    fimo.file <- sprintf("%s/%s_fimo_out.bz2", output.dir, f1a)
    fimo.out <- read.in.fimo( fimo.file, p.cutoff )

    motifs <- unique( paste( fimo.out$bic, fimo.out$mot, sep='_' ) )

    m <- lapply( 1:length( motifs ), function( i ) {
        m.tmp <- lapply( nc, function(i)integer() ); names( m.tmp ) <- names( nc )
        bi <- as.integer( strsplit( motifs[i], '_' )[[ 1 ]][ 1 ])
        mo <- as.integer( strsplit( motifs[i], '_' )[[ 1 ]][ 2 ])
        scans <- fimo.out[ J( c( bi, bi ), c( mo, -mo ) ) ]
        if ( nrow(scans) <= 0 ) return( m.tmp )
        scans <- scans[ ! is.na( scans$Start ), ] ##posns ), ]
        if ( nrow(scans) <= 0 ) return( m.tmp ) 
        for ( iii in 1:nrow( scans ) ) {
            inds <- scans$Start[ iii ]:scans$Stop[ iii ]
            chr <- levels( scans$Seq )[ scans$Seq[ iii ] ]
            m.tmp[[ chr ]] <- c( m.tmp[[ chr ]], inds )
        }
        ##cat( i, f1a, length( motifs ), motifs[ i ], nrow( scans ), "\n" ) ##sum( m.tmp ), "\n" )
        lapply( m.tmp, unique )
    } )
    names( m ) <- motifs
    save( m, file=sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) )
    m
}    
    
##motif.shadows <-
mclapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
    do.motif.shadows(f1)
}, mc.preschedule=F )

##names( motif.shadows ) <- rdatas
