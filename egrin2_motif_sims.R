######################################################################
######################################################################
######################################################################
#########################################################################
## Next compute the overlap (similarities) of all "motif shadows" vs all others - save as motif_sims.tsv.bz2
## This takes a while as it is all runs vs. all runs

do.motif.sims <- function(f1) {
    ##m1 <- motif.shadows[[f1]]
    f1a <- gsub( '.RData', '', f1 )
    if ( ! file.exists( sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) ) ) return()
    load( sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) )
    m1 <- m; rm( m )
        
    f1 <- gsub( '.RData', '', basename(f1) )
    mclapply( i:length(rdatas), function(j) { ## allow self-self comparisons.
        f2 <- rdatas[j]    
        ##m2 <- motif.shadows[[f2]]
        f2a <- gsub( '.RData', '', f2 )
        if ( ! file.exists( sprintf('%s/%s_motif_shadows.RData', output.dir, f2a) ) ) return()
        load( sprintf('%s/%s_motif_shadows.RData', output.dir, f2a) )
        m2 <- m; rm( m )
        f2 <- gsub( '.RData', '', basename(f2) )
        print( sprintf('%s/%s_vs_%s_motif_sims.tsv.bz2', output.dir, f1, f2) )
        if ( file.exists( sprintf('%s/%s_vs_%s_motif_sims.tsv.bz2', output.dir, f1, f2) ) ) return()
        motif.sims <- data.frame()
        for ( ii in 1:length(m1) ) {
            m11 <- m1[[ii]]
            for ( jj in 1:length(m2) ) {
                m22 <- m2[[jj]]
                if ( length(m11) != length(m22) ) { print( "WARNING: WRONG LENGTH!" ) } ## length should be # of chroms
                sum1 <- 0; sum2 <- 0; sum3 <- 0
                for ( iii in 1:length(m11) ) {
                    x <- m11[[iii]]
                    y <- m22[[iii]]
                    if ( length( x ) <= 0 || length( y ) <= 0 ) next
                    tmp <- x %in% y
                    sum1 <- sum1 + sum( ! tmp )
                    sum2 <- sum2 + sum( ! ( y %in% x ) )
                    sum3 <- sum3 + length( unique( c( x, y ) ) ) 
                }
                tmp <- ( sum1 + sum2 ) / sum3 ## this is a DISTANCE, so 1 is NO SIMILARITY
                if ( tmp < 1 ) {
                    if ( tmp < 0.9 ) cat( f1, names(m1)[ii], f2, names(m2)[jj], tmp, '\n' )
                    motif.sims <- rbind( motif.sims, data.frame( m1=names(m1)[ii], m2=names(m2)[jj], dist=tmp ) )
                }
            }
        }
        write.table( motif.sims, file=sprintf('%s/%s_vs_%s_motif_sims.tsv', output.dir, f1, f2),
                    sep='\t', quote=F, row=F )
        system( sprintf('bzip2 -fv %s/%s_vs_%s_motif_sims.tsv', output.dir, f1, f2) )
    }, mc.preschedule=F )
}

for ( f1 in 1:length(rdatas) ) {
    f1 <- rdatas[f1]
    do.motif.sims(f1)
}
