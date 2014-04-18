######################################################################
######################################################################
######################################################################
#########################################################################
## Do the motif clustering based on positions in the fimo.outs (computed above).
##  note this is a new method, not using tomtom, but computing all the fimo scans is
##  much faster so we should make this work!
## Next cluster the motifs based on the overlap (similarities) of all "motif shadows"
##   using mcl

if ( ! exists( 'distance.weight.cutoff' ) ) distance.weight.cutoff <- 0.95   ## 0.99 ??
if ( ! exists( 'mcl.I' ) ) mcl.I <- 4.5
if ( ! exists( 'mcl.pi' ) ) mcl.pi <- 2.0
if ( ! exists( 'n.cutoff' ) ) n.cutoff <- 10

## Force it to be stopped after 15 minutes, and killed 1m later if still running...
mcl.cmd <- paste( 'timeout -k 1m 15m ./progs/mcl-10-201/local/bin/mcl %s/new.mcltemp -o %s --abc',
                 ' -I %.1f -v all -te 3 -S 200000 -scheme 7 --analyze=y -pi %.1f 2>&1' )

cat( "Going to run mcl now...\n" )
fs <- list.files( path=output.dir, patt=glob2rx('*_motif_sims.tsv.bz2'), full=F )
for ( f in fs ) {
    append <- '>>'
    if ( f == fs[1] ) append <- '>'
    tmp <- strsplit( f, '_' )[[ 1 ]]
    system( sprintf( "bunzip2 -cv %s/%s | awk '($3<=%.2f){print \"%s_\"$1,\"%s_\"$2,1-$3}' %s %s/new.mcltemp", 
                    output.dir, f, distance.weight.cutoff, tmp[3], tmp[7], append, output.dir ) )
    system( sprintf( 'wc -l %s/new.mcltemp', output.dir ) )
}

outfile <- sprintf( "%s/new.mcltemp.I%s.pi%s", output.dir, gsub('.','',sprintf("%.1f",mcl.I),fixed=T),
                   gsub('.','',sprintf("%.1f",mcl.pi),fixed=T) )
mcl.cmd.in <- sprintf( mcl.cmd, output.dir, outfile, mcl.I, mcl.pi )
print( mcl.cmd.in )
mcl.out <- system( mcl.cmd.in, intern=T, ignore.stderr=F )
##if ( ! improve.the.clusters ) unlink( 'new.mcltemp' )
print( mcl.out )

system( sprintf( 'gzip -fv %s', outfile ) )
clusts <- strsplit( readLines( gzfile( sprintf( "%s.gz", outfile ) ) ), "\t" )
cat( "GOT", length( clusts ), "motif clusters with", length( unlist( clusts ) ), "motifs.\n" )
clusts <- clusts[ sapply( clusts, length ) >= 3 ] ## eliminate tiny clusters
cat( "GOT", sum( sapply( clusts, length ) >= 3 ), "motif clusters (length > 3) with",
    length( unlist( clusts[ sapply( clusts, length ) >= 3 ] ) ), "motifs.\n" )
cat( "GOT", sum( sapply( clusts, length ) >= 10 ), "motif clusters (length > 10) with",
    length( unlist( clusts[ sapply( clusts, length ) >= 10 ] ) ), "motifs.\n" )
mc.length <- max( which( sapply( clusts, length ) >= n.cutoff ) ) ## ignore small clusters
if ( is.infinite( mc.length ) ) mc.length <- max( which( sapply( clusts, length ) >= 3 ) )
attr( clusts, 'mcl.cmd' ) <- mcl.cmd.in
attr( clusts, 'mcl.out' ) <- mcl.out
attr( clusts, 'mc.length' ) <- mc.length
save( clusts, file=sprintf( "%s/new_motif_shadows_%.1f_%.1f_%.1f_clusts.RData",
                  output.dir, distance.weight.cutoff, mcl.I, mcl.pi ) )

## Compute bad motif clusters by seeing which of its motifs have a high 'coding.frac'.
## Assumes this script was run with 'coding.fracs' in 'tasks' so 'coding.fracs' exists.
if ( exists( 'coding.fracs' ) ) {
    cf.names <- sapply( strsplit( gsub( '.RData', '', basename( names( coding.fracs ) ) ), '_' ), '[', 3 )
    bad.clusts <- character()
    for ( i in 1:attr(clusts,'mc.length') ) {
        cat('BAD MOTIF CLUSTS', i, '\n')
        c <- clusts[[ i ]]
        tmp <- strsplit( c, '_' )
        cfs <- sapply( tmp, function( t ) {
            cf <- coding.fracs[[ which( cf.names == t[1] ) ]]
            cf2 <- cf$all.fracs[[ paste( t[2:3], collapse='_' ) ]] ##<= cf$mean.fracs - 0.01
            cf2
        } )
        mn <- mean( cfs, na.rm=T )
        sd <- sd( cfs, na.rm=T )
        if ( is.na( mn ) || is.na( sd ) || mn + sd * 2 > coding.fracs[[1]]$mean.fracs )
            bad.clusts <- c( bad.clusts, paste( 'MOTC', i, sep='_' ) )
        ##cat( i, paste( 'MOTC', i, sep='_' ) %in% bad.clusts, "\n" )
    }

    attr( clusts, 'bad.clusts' ) <- bad.clusts
}

save( clusts, file=sprintf( "%s/new_motif_shadows_%.1f_%.1f_%.1f_clusts.RData",
                  output.dir, distance.weight.cutoff, mcl.I, mcl.pi ) )
