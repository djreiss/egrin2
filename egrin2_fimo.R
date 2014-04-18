######################################################################
######################################################################
######################################################################
#########################################################################
## First fimo all motifs from each run against the genome and save the output...

cmd <- "./progs/fimo --max-stored-scores 9999999 --max-seq-length 1e8 --text --verbosity 2 %s %s" ## | bzip2 -c >%s.bz2"

do.fimo <- function(f1) {
    cat('FIMO', f1, '\n')
    f1a <- gsub( '.RData', '', basename(f1) )
    out.file <- sprintf("./%s/%s_fimo_out", output.dir, f1a)
    if ( file.exists( out.file ) || file.exists( sprintf( "%s.bz2", out.file ) ) ) { print( "SKIPPING" ); return() }
    load( sprintf( '%s/%s', output.dir, f1 ) )
    e$all.motifs.to.meme.file <- all.motifs.to.meme.file; environment(e$all.motifs.to.meme.file) <- e

    dir.create( sprintf( './%s/fimo_out', output.dir ) )
    seqs.file <- e$my.tempfile("fimo_seqs", tmpdir=sprintf( "./%s/fimo_out", output.dir ) )
    writeLines( paste( paste(">", names(e$genome.info$genome.seqs), sep=""),
                      e$genome.info$genome.seqs, sep="\n"), con=seqs.file )
    
    inds <- c( seq( 1, e$k.clust, by=20 ), e$k.clust ) ## fimo can't handle more than about 100 motifs!!??
    for ( i in 1:( length( inds ) - 1 ) ) {
        mots.file <- e$all.motifs.to.meme.file(ks=inds[i]:inds[i+1], e.val=Inf, resid=Inf)
        system( sprintf('rpl bic_ %s_ %s', f1a, mots.file) )
        if ( i == 1 && ( file.exists( out.file ) || file.exists( sprintf( "%s.bz2", out.file ) ) ) ) {
            print( "SKIPPING" ); return() }
        cmd1 <- sprintf(cmd, mots.file, seqs.file)
        if ( i == 1 ) cmd1 <- sprintf( "%s >%s", cmd1, out.file )
        else cmd1 <- sprintf( "%s >>%s", cmd1, out.file )
        print(cmd1)
        out <- system(cmd1, intern=T)
    }
    ## Need to remove the headers that are inside the output file, other than the first one
    tfile <- e$my.tempfile('fimo_tmp', tmpdir=sprintf('./%s/fimo_out', output.dir))
    system( sprintf( 'head -1 %s >%s', out.file, tfile ) )
    system( sprintf( 'grep -v Motif %s >>%s', out.file, tfile ) )
    system( sprintf( 'mv -fv %s %s', tfile, out.file ) )
    system( sprintf( 'bzip2 -vf %s', out.file ) )
}    

##for ( f1 in 1:length(rdatas) ) {
mclapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
    do.fimo( f1 )
}, mc.preschedule=F )
