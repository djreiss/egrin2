######################################################################
######################################################################
######################################################################
#########################################################################
## Let's change how we do tomtom-ing -- let's tomtom (and save the results of) each
##    output RData file (from above) against each other. That way we can do it during
##    the cMonkey running.

e.value.cutoff=100; resid.cutoff=0.8; dist.meth="ed"; q.thresh=0.5; min.overlap=4; q.pseudo=0; t.pseudo=0
cmd <- "%s/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"
## see http://meme.nbcr.net/meme/doc/tomtom.html for output format of tomtom

##for ( f1 in 1:length(rdatas) ) {
lapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
    cat('TOMTOM', f1, '\n')
    load( sprintf( '%s/%s', output.dir, f1 ) )
    e1 <- e; rm( e )
    ##sys.source( "~/scratch/biclust/cmonkey-motif-other.R", envir=e1 ) ## make sure the env has "all.motifs.to.mast.file()"
    e1$all.motifs.to.meme.file <- all.motifs.to.meme.file; environment(e1$all.motifs.to.meme.file) <- e1
    tf1 <- e1$all.motifs.to.meme.file(e.val=Inf)

    f1a <- basename( gsub( '.RData', '', f1 ) )
    
    system( sprintf('rpl bic_ %s_ %s', f1a, tf1) )

    rdatas.2 <- list.files( path=output.dir, patt=glob2rx('zzz_eco_???.RData'), full=F )
    ##for ( f2 in 1:length(rdatas.2) ) { ## all vs. all including self vs self -- note should use mclapply to parallelize
    mclapply( 1:length(rdatas.2), function(f2) {
        f2 <- rdatas.2[f2]
        cat('TOMTOM2', f2, '\n')
        f2a <- gsub( '.RData', '', basename(f2) )
        fname <- sprintf('%s/%s_vs_%s_tomtom.tsv', output.dir, f1a, f2a)
        if ( file.exists( sprintf('%s.bz2', fname) ) ) {
            print( sprintf('SKIPPING: %s', fname) ); return() } ## already did this comparison
        if ( file.exists( sprintf('%s/%s_vs_%s_tomtom.tsv.bz2', output.dir, f2a, f1a) ) ) {
            print( sprintf('SKIPPING (REVERSED): %s', fname) ); return() } ## symmetric -- already did this comparison
        load( sprintf( '%s/%s', output.dir, f2 ) )
        e2 <- e; rm( e )
        ##sys.source( "~/scratch/biclust/cmonkey-motif-other.R", envir=e2 )
        e2$all.motifs.to.meme.file <- all.motifs.to.meme.file; environment(e2$all.motifs.to.meme.file) <- e2
        tf2 <- e2$all.motifs.to.meme.file(e.val=Inf)
        system( sprintf('rpl bic_ %s_ %s', f2a, tf2) )
        if ( file.exists( sprintf('%s.bz2', fname) ) ) {
            print( sprintf('SKIPPING: %s', fname) ); return() } ## someone else is already doing this comparison
        if ( file.exists( sprintf('%s/%s_vs_%s_tomtom.tsv.bz2', output.dir, f2a, f1a) ) ) {
            print( sprintf('SKIPPING (REVERSED): %s', fname) ); return() } ## symmetric -- already did this comparison
        cmd.1 <- sprintf( cmd, e1$progs.dir, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tf2 )
        cmd.1 <- paste( cmd.1, "-query", tf1, '|bzip2 -c >', sprintf('%s.bz2', fname) )
        print( cmd.1 )
        system( cmd.1 ) ## note this has not been filtered to remove self-self comparisons

        if ( FALSE ) { ## used to be - read in results and save out processed file.
            tout <- system( cmd.1, intern=T )
            tout <- do.call( rbind, strsplit( tout, "\t" ) )
            colnames( tout ) <- tout[ 1, ,drop=F ]; tout <- tout[ -1, ,drop=F ]
            tout <- as.data.frame( tout[ tout[ ,1 ] != tout[ ,2 ], ,drop=F ] ) ## de-symmetrize the output
            cat( "TOUT:", dim( tout ), "\n" )
            ## Save the file
            write.table( tout, row.names=F, sep='\t', quote=F, file=sprintf('%s/%s_vs_%s_tomtom.tsv', output.dir, f1a, f2a) )
            system( sprintf( 'bzip2 -fv %s/%s_vs_%s_tomtom.tsv', output.dir, f1a, f2a) )
        }
    }, mc.preschedule=F )
    
    rm(e,e1,e2,tout)
} )
