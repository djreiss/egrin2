require( data.table )
require( parallel )
debug.on()

## if ON.CLUSTER == TRUE then we are running it on an SGE cluster
if ( ! exists( 'ON.CLUSTER' ) ) ON.CLUSTER <- FALSE

if ( ! exists( 'tasks' ) || ! exists( 'rdatas' ) ) {
    ## If this script was run from zzz_eco4_tomtom_onSGE.R then mc.cores has already been set.
    options( mc.cores=12 )
}

## Let's do all of the tomtom-ing, and fimo-ing, etc. on an individual RData basis, rather than after
## integrating all RData files (from different cMonkey runs) into a single big cMonkey run...

if ( ! exists( 'tasks' ) )
    tasks <- c( 'fimo', 'coding.fracs', 'motif.shadows', 'nwinf', 
               'regdb', 'motif.regdb.aupr',  ## compares MOTIFS to regdb (eco only!)
               'nwinf.regdb.aupr',           ## compares Inf. to regdb   (eco only!)
               'motif.sims', 'newmotcluster',
               'corems', 'tomtom' )

if ( ! exists( 'output.dir' ) ) output.dir <- 'output/'

if ( ! exists( 'rdatas' ) ) { ## allows us to just run it on one file
    rdatas <- list.files( path=output.dir, patt=glob2rx('zzz_eco_???.RData'), full=F )
}

if ( ! exists( 'rdatas2' ) ) {  ## second list of rdatas with which to compare motifs against (double loop) (tasks 'motif.sims' or 'tomtom')
    ## by default if rdatas is already provided and it is one file, then assume we're comparing the one file against all others (for running on cluster)
    if ( length( rdatas ) == 1 ) rdatas2 <- list.files( path=output.dir, patt=glob2rx('zzz_eco_???.RData'), full=F )
    else rdatas2 <- rdatas
}

if ( ! exists( 'do.plot' ) ) do.plot <- FALSE ## plot aupr statistics? (tasks 'motif.regdb.aupr' and 'nwinf.regdb.aupr')
if ( ! exists( 'gold.option' ) ) gold.option <- 'dream5' ## other options: 'strong+weak', 'strong', 'weak', 'confirmed'

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
    fimo.out <- as.data.frame( do.call.rbind( lapply( fimo.out, function( i ) strsplit( i, '\t' )[[ 1 ]] ) ) )
    colnames( fimo.out ) <- strsplit( system( sprintf( 'bunzip2 -c %s | head -1', fimo.file ), intern=T ), '\t' )[[ 1 ]]
    fimo.out$Start <- as.integer( as.character( fimo.out$Start ) )
    fimo.out$Stop <- as.integer( as.character( fimo.out$Stop ) )
    fimo.out$`Log-odds` <- as.numeric( as.character( fimo.out$`Log-odds` ) )
    fimo.out$`p-value` <- as.numeric( as.character( fimo.out$`p-value` ) )
    tmp <- do.call.rbind( strsplit( as.character( fimo.out$Motif ), '_' ) )
    fimo.out$bic <- as.integer( tmp[ ,4 ] )
    fimo.out$mot <- as.integer( tmp[ ,5 ] )
    fimo.out$Strand <- substr( tmp[ ,1 ], 1, 1 )
    rm( tmp )
    fimo.out$Motif <- NULL
    fimo.out <- as.data.table( fimo.out )
    setkey( fimo.out, bic, mot, Seq, Start )
    fimo.out
}

get.all.clusterStacks <- function() {
    all.rdata.clusterStacks <- lapply( rdatas, function( rd ) {
        cat( 'CLUSTERSTACKS', rd, '\n' )
        filename <- sprintf( '%s/%s_clusterStacks.RData', output.dir, gsub('.RData', '', rd) )
        if ( file.exists( filename ) ) { load( filename ); return( clusterStack ) }
        
        load( sprintf( '%s/%s', output.dir, rd ) )
        clusterStack <- e$clusterStack
        save( clusterStack, file=filename )
        clusterStack
    } )
    names( all.rdata.clusterStacks ) <- rdatas
    all.rdata.clusterStacks
}

## Get meme.out structure each motif in each bicluster in each cmonkey run (contains pssms)
## save individually
get.all.meme.outs <- function( seq.type='upstream meme' ) {
    all.rdata.meme.outs <- lapply( rdatas, function( rd ) {
        cat( 'MEME OUTS', rd, '\n' )
        filename <- sprintf( '%s/%s_meme_outs.RData', output.dir, gsub('.RData', '', rd) )
        if ( file.exists( filename ) ) { load( filename ); return( meme.outs ) }
        
        load( sprintf( '%s/%s', output.dir, rd ) )
        meme.outs <- lapply( e$meme.scores[[ seq.type ]], function( out ) {
            out$pv.ev <- NULL; out$prev.run <- NULL; out } )
        meme.outs$all.pv <- NULL
        print(object.size(meme.outs))
        save( meme.outs, file=filename )
        meme.outs
    } )
    names( all.rdata.meme.outs ) <- rdatas
    all.rdata.meme.outs
}

if ( length(tasks) <= 0 ) stop()

######################################################################
######################################################################
######################################################################
#########################################################################
## First fimo all motifs from each run against the genome and save the output...

if ( 'fimo' %in% tasks ) {
cmd <- "./progs/fimo --max-stored-scores 9999999 --max-seq-length 1e8 --text --verbosity 2 %s %s" ## | bzip2 -c >%s.bz2"

##for ( f1 in 1:length(rdatas) ) {
mclapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
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
    
    inds <- c( seq( 1, e$k.clust, by=100 ), e$k.clust ) ## fimo can't handle more than about 100 motifs!!??
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
}, mc.preschedule=F )

##rm(list=ls()); gc()
}

######################################################################
######################################################################
######################################################################
#################################################################
## Compute per-run fraction of times each motif is in a coding region, save the results

if ( 'coding.fracs' %in% tasks ) {
p.cutoff <- 1e-6

##for ( f1 in 1:length(rdatas) ) {
coding.fracs <- mclapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
    cat('CODING FRACS', f1, '\n')
    f1a <- gsub( '.RData', '', basename(f1) )
    out.file <- sprintf("./%s/%s_coding_fracs.RData", output.dir, f1a)
    if ( file.exists( out.file ) ) { load( out.file ); return( coding.fracs ) }
    fimo.file <- sprintf("./%s/%s_fimo_out.bz2", output.dir, f1a)
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
}, mc.preschedule=F )

names(coding.fracs) <- rdatas
}    

######################################################################
######################################################################
######################################################################
#########################################################################
## Compute the "motif shadows" for each run -- locations covered by each motif in each run -
##   save as motif_shadows.RData

if ( 'motif.shadows' %in% tasks ) {
p.cutoff <- 1e-4
    
mclapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
    cat('MOTIF SHADOWS', f1, '\n')
    f1a <- gsub( '.RData', '', f1 )
    if ( file.exists( sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) ) ) {
        load( sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) )
        return()
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
}, mc.preschedule=F )
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Next compute the overlap (similarities) of all "motif shadows" vs all others - save as motif_sims.tsv.bz2
## This takes a while as it is all runs vs. all runs

if ( 'motif.sims' %in% tasks ) {
p.cutoff <- 1e-6

apply.func <- lapply
if ( ! ON.CLUSTER ) apply.func <- mclapply

for ( f1 in rdatas ) {
    f1 <- gsub( '.RData', '', basename(f1) )
    if ( file.exists( sprintf('%s/%s_motif_shadows.RData', output.dir, f1) ) ) {
        load( sprintf('%s/%s_motif_shadows.RData', output.dir, f1) )
        m1 <- m; rm( m )
    } else {
        next
    }
    apply.func( rdatas2, function(f2) { ## allow self-self comparisons.
        f2 <- gsub( '.RData', '', basename(f2) )
        print( sprintf('%s/%s_vs_%s_motif_sims.tsv.bz2', output.dir, f1, f2) )
        if ( file.exists( sprintf('%s/%s_vs_%s_motif_sims.tsv.bz2', output.dir, f1, f2) ) ) return()
        if ( file.exists( sprintf('%s/%s_motif_shadows.RData', output.dir, f2) ) ) {
            load( sprintf('%s/%s_motif_shadows.RData', output.dir, f2) )
            m2 <- m; rm( m )
        } else {
            return()
        }
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
        if ( ON.CLUSTER ) stop()
    } ) ##, mc.preschedule=F )
}
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Do the motif clustering based on positions in the fimo.outs (computed above).
##  note this is a new method, not using tomtom, but computing all the fimo scans is
##  much faster so we should make this work!
## Next cluster the motifs based on the overlap (similarities) of all "motif shadows"
##   using mcl
if ( 'newmotcluster' %in% tasks ) {
    
if ( ! exists( 'distance.weight.cutoff' ) ) distance.weight.cutoff <- 0.95   ## 0.99 ??
if ( ! exists( 'mcl.I' ) ) mcl.I <- 4.5
if ( ! exists( 'mcl.pi' ) ) mcl.pi <- 2.0
if ( ! exists( 'n.cutoff' ) ) n.cutoff <- 10
if ( ! exists( 'exclude.bad.mots' ) ) exclude.bad.mots <- TRUE ## don't include "bad motifs"

## Force it to be stopped after 30 minutes, and killed 1m later if still running...
mcl.cmd <- paste( 'timeout -k 1m 30m ./progs/mcl-10-201/local/bin/mcl %s -o %s --abc',
                 ' -I %.1f -v all -te 3 -S 200000 -scheme 7 --analyze=y -pi %.1f 2>&1' )

cat( "Going to run mcl now...\n" )
fs <- list.files( path=output.dir, patt=glob2rx('*_motif_sims.tsv.bz2'), full=F )
dir.create( 'zzztemp' )
mcltempfile <- tempfile( tmpdir='./zzztemp' )
for ( f in fs ) {
    cat( f, distance.weight.cutoff, '\n' )
    append <- '>>'
    if ( f == fs[1] ) append <- '>'
    tmp <- strsplit( f, '_' )[[ 1 ]]
    system( sprintf( "pbzip2 -dc %s/%s | awk '($3<=%.9f){print \"%s_\"$1,\"%s_\"$2,1-$3}' %s %s", #new.mcltemp", 
                    output.dir, f, distance.weight.cutoff, tmp[3], tmp[7], append, mcltempfile ) )
    ##system( sprintf( 'wc -l %s', mcltempfile ) )
}

if ( exclude.bad.mots && exists( 'coding.fracs' ) ) {
    cfs <- do.call( rbind, lapply(names(coding.fracs), function(rd) {
        if ( is.null( coding.fracs[[rd]] ) ) return(NULL)
        tmp <- strsplit( gsub('.RData','',rd), '_' )[[ 1 ]]
        data.table( paste(tmp[3],names(coding.fracs[[rd]]$all.fracs),sep='_'), coding.fracs[[rd]]$all.fracs)
    } ) )
    cfs <- subset( cfs, V2 <= coding.fracs[[1]]$mean.fracs + 0.02 )
    tmp <- fread( mcltempfile ) ##sprintf('%s/new.mcltemp', output.dir) )
    tmp <- subset( tmp, V1 != V2 & V1 %in% cfs$V1 & V2 %in% cfs$V1 )
    write.table( tmp, sep=' ', col=F, row=F, quote=F, file=mcltempfile ) ##sprintf('%s/new.mcltemp',output.dir) )
}

print(mcltempfile);system(sprintf("wc -l %s",mcltempfile))

outfile <- sprintf( "%s/new.mcltemp.I%s.pi%s.dc%.9f", output.dir, gsub('.','',sprintf("%.1f",mcl.I),fixed=T),
                   gsub('.','',sprintf("%.1f",mcl.pi),fixed=T), distance.weight.cutoff )
mcl.cmd.in <- sprintf( mcl.cmd, mcltempfile, outfile, mcl.I, mcl.pi )
print( mcl.cmd.in )
mcl.out <- system( mcl.cmd.in, intern=T, ignore.stderr=F )
##if ( ! improve.the.clusters ) unlink( 'new.mcltemp' )
print( mcl.out )

unlink( mcltempfile )
if ( ! file.exists( outfile ) ) stop() ## must have taken longer than 15 minutes.

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
save( clusts, file=sprintf( "%s/new_motif_shadows_%.1f_%.1f_%.9f_clusts.RData",
                  output.dir, mcl.I, mcl.pi, distance.weight.cutoff ) )

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

save( clusts, file=sprintf( "%s/new_motif_shadows_%.1f_%.1f_%.9f_clusts.RData",
                  output.dir, mcl.I, mcl.pi, distance.weight.cutoff ) )
}    

######################################################################
######################################################################
######################################################################
#########################################################################
## Compute combined pssms for each motif cluster. Here we need to run tomtom
##  between all of the motifs in each motif cluster.
## It will also plot them if an individual cmonkey run "e" is in the environment

if ( 'newmotcluster.getcombinedpssms' %in% tasks && exists( 'clusts' ) ) {

cf.names <- sapply( strsplit( gsub( '.RData', '', basename( rdatas ) ), '_' ), '[', 3 )

e.value.cutoff=100; resid.cutoff=0.8; dist.meth="ed"; q.thresh=0.5; min.overlap=4; q.pseudo=0; t.pseudo=0
cmd <- "./progs/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"
max.clust.size <- 800

if ( ! exists( 'all.rdata.meme.outs' ) ) all.rdata.meme.outs <- get.all.meme.outs()

get.combined.pssm <- function( tt.out, seq.type='upstream meme', n.gene.weight=F, return.aligned.pssms=F ) {
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
## Get meme file output for each motif in each bicluster in each cmonkey run (a hierarchical list)
## save individually
all.rdata.meme.files <- mclapply( rdatas, function( rd ) {
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
} )
names( all.rdata.meme.files ) <- rdatas
}

##for ( i in 1:length( clusts ) ) {
new.clusts <- mclapply( clusts[ 1:attr(clusts,'mc.length')], function( clust ) {
    ##clust <- clusts[[ i ]]
    print( length( clust ) )
    if ( length( clust ) > max.clust.size || ! is.null( attr( clust, 'aligned.pssm' ) ) ) return( clust )
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
    attr( clust, 'tt.out' ) <- tt.out
    pssm <- get.combined.pssm( tt.out )
    attr( clust, 'combined.pssm' ) <- pssm
    #if ( exists( 'e' ) ) e$viewPssm( pssm, main=paste( length(clust),
    #                     length(unique(c(as.character(tt.out$`Query ID`),as.character(tt.out$`Target ID`)))) ) )
    clust
}, mc.preschedule=F )

for ( i in 1:length(new.clusts) ) if (class(new.clusts[[i]]) == 'try-error') new.clusts[[i]] <- clusts[[i]]
for ( n in names( attributes( clusts ) ) ) attr( new.clusts, n ) <- attr( clusts, n )
}


######################################################################
######################################################################
######################################################################
#########################################################################
## Now here is the comparison vs. regulondb part -- this compares motif locations with
##   binding sites in regulondb. Do it individually for each run and save the results
## NOTE THIS IS SPECIFIC TO E. COLI ONLY !!!

if ( 'regdb' %in% tasks ) {

binding.site.set.file <- sprintf( '%s/BindingSiteSet.txt', getwd() )
if ( ! file.exists( binding.site.set.file ) )
    system( 'wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt' )
    
p.cutoff <- 1e-6
require( data.table )

##for ( f1 in 1:length(rdatas) ) {
regdb.out <- mclapply( 1:length(rdatas), function(f1) {
    f1 <- basename(rdatas[f1])
    cat('REGDB', f1, '\n')
    f1a <- gsub( '.RData', '', f1 )
    if ( file.exists( sprintf('%s/%s_with_regdb_comparison.RData', output.dir, f1a) ) ) {
        load( sprintf('%s/%s_with_regdb_comparison.RData', output.dir, f1a) )
        cat(f1,length(combi.stats),'\n')
        return( combi.stats )
    }
    system( sprintf('touch %s/%s_with_regdb_comparison.RData', output.dir, f1a) ) ## placeholder
    load( sprintf( '%s/%s', output.dir, f1 ) )
    rm( ratios ); gc()

    fimo.file <- sprintf("%s/%s_fimo_out.bz2", output.dir, f1a)

    ## Pre-read in the fimo.out computed above. This is mostly taken from compare_regulondb_single.R ...
    fimo.out <- read.in.fimo( fimo.file, p.cutoff )

    cwd <- setwd('cmEvalEco')
    try(capture.output(try(source("compare_regulondb_single.R",local=T))))

    use.combi <- TRUE
    source('compare_regdb_regulons.R',local=T)
    combi.stats <- gsub( '"', '', gsub( '[1] ', '', capture.output( source("print_stats.R",local=T) ), fixed=T ) )
    print( combi.stats )

    use.combi <- FALSE
    source('compare_regdb_regulons.R',local=T)
    non.combi.stats <- gsub( '"', '', gsub( '[1] ', '', capture.output( source("print_stats.R",local=T) ), fixed=T ) )
  
    setwd(cwd)
    
    save( combi.stats, non.combi.stats, eco.hits, eco.hits2, file=sprintf('%s/%s_with_regdb_comparison.RData',
                                                                 output.dir, f1a) )
    return( combi.stats )
}, mc.preschedule=F )
names( regdb.out ) <- rdatas
n.mots.hit <- as.integer( sapply(strsplit(sapply(regdb.out,'[',12),' '),'[',13) )
hist( n.mots.hit, breaks=20 )
}


######################################################################
######################################################################
######################################################################
## The stuff below is required for AUPR functions (nwInf AND motif)

## Updated AUPR function based on Sriram's matlab code
## Assumes input is tf, target, weight and gold is tf, target
## Only uses Inf. predictions that have abs(weight) >= 0.1 ... see ../FINAL_compare_regulondb.R
get.aupr <- function( net, gold, weight.cut=0.1, n.cut=100000, plot.it=T, precs=0.25, title='', ... ) {
    require(caTools) ## for trapz

    net <- net[ rev( rank( abs(net$weight), ties='random' ) ), ]
    if ( nrow( net ) > n.cut ) net <- net[ 1:n.cut, ]
    net <- net[ abs(net$weight) >= weight.cut, ]

    if ( 'weight' %in% colnames(gold) ) {
        gold <- gold[ rev( rank( abs(gold$weight), ties='random' ) ), ]
        gold <- gold[ abs(net$weight) >= weight.cut, ]
    }
    if ( nrow( gold ) > n.cut ) gold <- gold[ 1:n.cut, ]
    
    netwo = paste( as.character(net$tf), as.character(net$target), sep='//' ) 
    regulator = as.character(net$tf)
    targets = as.character(net$target)

    reg1 <- as.character(gold$tf)
    tg1 <- as.character(gold$target)
    allnet = paste( as.character(gold$tf), as.character(gold$target), sep='//' ) 

    ## overlap different network sizes with regulondb
    ## leng gives the network size
    leng = unique(c(seq(10,length(netwo),length=20),10,50,100,200,250,500,1000,1500,2000,2500,5000,
        floor(length(netwo)/10),floor(length(netwo)/4),floor(length(netwo)/2),length(netwo)))
    leng <- leng[leng <= length(netwo)] #<- NA ##'';
    leng = sort(leng);
    resultss = matrix(NA, nrow=length(leng), ncol=3) #zeros(1,4);
    le = length(allnet)
    for (jj in 1:length(leng)) {
        ix = regulator[1:leng[jj]] %in% reg1 #ismember(regulator(1:leng(jj)),reg1);
        ix1 = targets[1:leng[jj]] %in% tg1 #ismember(targets(1:leng(jj)),tg1)
        ix2 = (ix & ix1);
        netwo2 = netwo[1:leng[jj]]
        netwo1 = netwo2[ix2];  # gives the subset that has the same tf/targets as regulondb
        resultss[jj,3-1] = sum(netwo1%in%allnet, na.rm=T)/sum(ix2, na.rm=T) ##ismember(netwo1,allnet))/sum(ix2); % precision
        resultss[jj,4-1] = sum(allnet%in%netwo1, na.rm=T)/le ##ismember(allnet,netwo1))/le; % recall
        resultss[jj,2-1] = leng[jj]; 
        ##resultss[jj,1] = ii;
        ##c = c + 1;
    }
    colnames(resultss) <- c( 'n', 'precision', 'recall' )

    resultss[is.na(resultss)] <- 0 ## is this right?
    
  ## calculating auc
  ## get precision recall for each network 
  X = c(1,resultss[,2],0);
  Y = c(0,resultss[,3],1);
  IX = order(Y) # sort values
  X = X[IX]; Y = Y[IX];
  AUC = trapz(Y,X); # area under curve
  attr( resultss, 'AUC' ) <- AUC
  attr( resultss, 'n' ) <- length(netwo)

  ## get # correct predictions at precision of 25%
  tmp <- sapply(precs, function(prec) approx(resultss[,'precision'],resultss[,'recall'],prec)$y) * length(allnet)
  names(tmp) <- as.character(precs)  
  attr( resultss, 'nPredAt' ) <- tmp
    # plotting it
  if ( plot.it ) plot(Y,X,typ='l',lwd=2, xlab='Recall', ylab='Precision',
                      main=sprintf('%s AUC=%.3f; Npred=%d', title, AUC, as.integer(round(tmp[1]))), ...) 
  
  resultss
}

combine.networks <- function( net1, net2 ) { ## add the weights if the edges are teh same
    tmp1 <- net1$weight
    names(tmp1) <- paste( net1$tf, net1$target, sep='//' )

    tmp2 <- net2$weight
    names(tmp2) <- paste( net2$tf, net2$target, sep='//' )
    
    ttmp <- names(tmp1) %in% names(tmp2)
    tmp3 <- tmp1[ names(tmp1)[ttmp] ] + tmp2[ names(tmp1)[ttmp] ]   ## Add together the weights???
    tmp3 <- c( tmp3, tmp1[ ! ttmp ] )
    tmp3 <- c( tmp3, tmp2[ ! names(tmp2) %in% names(tmp1) ] )

    out <- do.call.rbind( strsplit( names(tmp3), '//' ) ) ##do.call( rbind, strsplit( names(tmp3), '//' ) )
    out <- data.table( tf=out[,1], target=out[,2], weight=tmp3 )
    out <- out[ order( abs( out$weight ), decreasing=T ), ]
    out
}

get.b.numbers <- function( gene.names ) {
    gene.names <- tolower( as.character( gene.names ) )
    if ( 'flhdc' %in% gene.names ) gene.names <- c( gene.names, 'flhd', 'flhc' )
    if ( 'hipab' %in% gene.names ) gene.names <- c( gene.names, 'hipa', 'hipb' )
    if ( 'reldb' %in% gene.names ) gene.names <- c( gene.names, 'rele', 'relb' )
    if ( 'rcsab' %in% gene.names ) gene.names <- c( gene.names, 'rcsa', 'rcsb' )
    if ( 'h-ns' %in% gene.names ) gene.names <- c( gene.names, 'hns' )
    if ( 'mqsr-mqsa' %in% gene.names ) gene.names <- c( gene.names, 'mqsr', 'mqsa' )
    if ( 'rcsb-bglj' %in% gene.names ) gene.names <- c( gene.names, 'rcsb', 'bglj' )
    if ( 'relb-rele' %in% gene.names ) gene.names <- c( gene.names, 'rele', 'relb' )
    if ( 'maze-mazf' %in% gene.names ) gene.names <- c( gene.names, 'maze', 'mazf' )
    if ( 'ihf' %in% gene.names ) gene.names <- c( gene.names, 'ihfa', 'ihfb' )

    names <- fread( 'EcoData032614-204851.txt', header=T )
    names$Gene <- tolower( as.character( gsub( "'", '', names$Gene ) ) )
    names.syns <- strsplit( tolower( names$Syn ), ', ', fixed=T )
    gene.lookup <- character()
    for ( n in gene.names ) {
        tmp <- subset( names, Gene == n )
        if ( nrow( tmp ) == 1 ) gene.lookup[ n ] <- tmp$bnum
        else if ( nrow( tmp ) == 0 ) {
            tmp <- which( sapply( names.syns, function(i) n %in% i ) )
            if ( length( tmp ) == 1 ) {
                tmp <- names[ tmp, ]
                if ( ! tmp$Gene %in% gene.names ) gene.lookup[ n ] <- tmp$bnum
            } else {
                cat( "HERE2", n, nrow(tmp), '\n' )
            }
        } else {
            cat( "HERE1", n, nrow(tmp), '\n' )
        }
    }
    gene.lookup
}

## NOTE: to use gold2 (or gold2.weak or gold2.strong) set 'gold <- gold2.weak', for example.
## To run new auprs for gold2 (or others) make a new 'output.dir', link to the files in the original 'output.dir'
##    but remove the links to the "vsRegDB" files where the regdb comparisons are saved.
if ( ! exists( 'gold' ) ) {
get.gold.standard <- function( type='dream5' ) {
    if ( type == 'dream5' ) {
        gold <- fread( 'DREAM5_NetworkInference_Evaluation/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network3.tsv',
                      head=F )
        load( 'format_data.RData' )
        rm( ratios, ratios.normed, names2, names1, names, gene.names, gene.lookup, gene.lookup2, names.syns )
        gold$tf <- gene.lookup.fin[ as.character( gold$V1 ) ]
        gold$target <- gene.lookup.fin[ as.character( gold$V2 ) ]
        gold <- subset( gold, ! is.na( tf ) & ! is.na( target ) ) ## only removes 2 edges
    }
    else {
        ## Get the full regdb set, to use strong, strong+weak, weak evidence. See:
        ## http://regulondb.cs.purdue.edu/evidenceclasification
        ## I copied this table and saved it to 'evidence.tsv' via LibreOffice
        ##system( 'wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt' )
        gold2 <- fread( 'network_tf_gene.txt', head=F, skip=34 )
        gold2$tf <- get.b.numbers( unique( tolower( as.character( gold2$V1 ) ) ) )[ tolower( as.character( gold2$V1 ) ) ]
        gold2$target <- get.b.numbers( unique( tolower( as.character( gold2$V2 ) ) ) )[ tolower( as.character( gold2$V2 ) ) ]
        gold2 <- subset( gold2, ! is.na( tf ) & ! is.na( target ) ) ## removes 578 edges

        evidence <- fread( 'evidence.tsv', head=F )
        evidence.lookup <- as.character( evidence$V4 ); names( evidence.lookup ) <- as.character( evidence$V2 )
        evidences <- strsplit( gsub( '[', '', gsub( ']', '', as.character( gold2$V4 ), fixed=T ), fixed=T ), ', ' )
        isstrong <- sapply( evidences, function(i) any( evidence.lookup[i] == 'Strong' ) )
        isweak <- sapply( evidences, function(i) all( evidence.lookup[i] == 'Weak' ) )
        if ( type == 'strong+weak' ) gold <- gold2
        else if ( type == 'strong' ) gold <- gold2[ which(isstrong), ]
        else if ( type == 'weak' ) gold <- gold2[ which(!isstrong), ]

        stop()
        ## Use the "BindingSiteSet.txt" file which has complete listing and more complete evidence listings
        ##system( 'wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt' )
        ##bss=fread("BindingSiteSet_8.2.txt",sep='\t',skip=42,head=F)

        ## SELEX datasets
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/ShimadaT_TFBSs_CRP_21673794.txt' ) ## weak
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/ShimadaT_TFBSs_H-NS_21883529_gSELEX.txt' ) # strong
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/ShimadaT_TFBSs_H-NS_21883529_gSELEX_GEA.txt' ) # str
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/ShimadaT_TFBSs_LeuO_21883529_gSELEX.txt' ) # strong
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/ShimadaT_TFBSs_LeuO_21883529_gSELEX_GEA.txt' ) # str

        ## RNA-seq datasets
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/High_throughput_transcription_initiation_mapping_with_5_tri_or_monophosphate_enrichment_v3.0.txt' ) ## strong
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/Promoter_from_454_Dataset.txt' ) ## weak
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/Promoter_from_RACE_Dataset.txt' ) ## strong
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/StorzG_TSS_Table_LB_0.4.txt' ) ## strong
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/StorzG_TSS_Table_LB_2.0.txt' ) ## strong
        system( 'wget http://regulondb.ccg.unam.mx/menu/download/high_throughput_datasets/files/StorzG_TSS_Table_M63_0.4.txt' ) ## strong
        
        
    }
    gold
}

if ( gold.option != 'NONE' ) gold <- get.gold.standard( gold.option )

}

######################################################################
######################################################################
######################################################################
## Compute the motif-based networks and get AUPR measures

if ( 'motif.regdb.aupr' %in% tasks ) {

if ( ! exists( 'all.rdata.meme.outs' ) ) all.rdata.meme.outs <- get.all.meme.outs()

max.motifs <- Inf
hit.p.cutoff <- 1e-4
distance.cutoff <- Inf ##0.99

## Need to compute conversion from AaaX to b### number.
## try the table gotten from here: http://www.ecogene.org/?q=ecodownload/dbtable
## needed a bit of editing to get it to load -- replace the 'b#' with 'bnum' and the last tab in the first line

mot.networks <- list()
for ( f1 in 1:length(rdatas) ) {
    f1 <- basename(rdatas[f1])
    print(f1)
    f1a <- gsub( '.RData', '', f1 )
    if ( file.exists( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) ) ) {
        load( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) )
        mot.networks[[ f1 ]] <- network
        next
    }
    if ( ! file.exists( sprintf('%s/%s_with_regdb_comparison.RData', output.dir, f1a) ) ) next
    load( sprintf('%s/%s_with_regdb_comparison.RData', output.dir, f1a) )

    ms <- all.rdata.meme.outs[[f1]] ##e$meme.scores[[ seq.type ]]
    e.vals <- do.call( rbind, lapply( 1:max(unlist(sapply(ms,'[[','k'))), function( i ) {
        out <- c( NA, NA, NA )
        if ( is.null( ms[[ i ]] ) || ms[[ i ]] == '' ) return( out )
        mm <- ms[[ i ]]$meme.out
        for ( j in 1:length( mm ) ) if ( ! is.null( mm[[ j ]]$e.value ) ) out[ j ] <- mm[[ j ]]$e.value
        out
    } ) )

    eval.cutoff <- Inf
    if ( ! is.infinite( max.motifs ) ) eval.cutoff <- sort( e.vals[!is.na(e.vals)] )[ max.motifs ]
    tmp <- subset( eco.hits2, distance <= distance.cutoff ) ## we used 0.99 for the cutoff for the ensemble analysis
    tmp <- tmp[ order( tmp$distance ), ]
    b.m <- t( apply( sapply( strsplit( as.character( tmp$motif ), '_' ), '[', 2:3 ), 2, as.integer ) )
    tmp <- tmp[ e.vals[ b.m ] <= eval.cutoff, ]
    tmp <- subset( tmp, ! duplicated( tmp$motif ) )  ## this gives tfs to motifs match

    b.nums <- get.b.numbers(unique(tolower(as.character(tmp$eco.tf))))
    tmp$eco.tf2 <- b.nums[ tolower( as.character( tmp$eco.tf ) ) ]
    
    ## now need to get motifs to gene match
    network <- data.table()
    for ( i in 1:nrow(tmp) ) {
        m <- as.character(tmp$motif[i])
        ttmp <- as.integer(strsplit(m,'_')[[1]][2:3])  ## k, mot
        hits <- ms[[ ttmp[1] ]]$meme.out[[ ttmp[2] ]]$posns
        if ( nrow(hits) <= 0 || is.null(hits) ) { cat( m, 'NO GENES\n' ); next }
        hits <- subset( hits, p.value <= hit.p.cutoff )
        if ( nrow(hits) <= 0 ) { cat( m, 'NO GENES\n' ); next }
        genes <- as.character(hits$gene)
        tf <- tmp$eco.tf2[i]
        if (is.na(tf)) {
            tf <- as.character(tmp$eco.tf[i])
            if ( tf == 'FlhDC' ) tf <- b.nums[ c( 'flhd', 'flhc' ) ]
            if ( tf == 'HipAB' ) tf <- b.nums[ c( 'hipa', 'hipb' ) ]
            if ( tf == 'RelEB' ) tf <- b.nums[ c( 'rele', 'relb' ) ]
            if ( tf == 'RelB-RelE' ) tf <- b.nums[ c( 'rele', 'relb' ) ]
            if ( tf == 'MazE-MazF' ) tf <- b.nums[ c( 'maze', 'mazf' ) ]
            if ( tf == 'RcsAB' ) tf <- b.nums[ c( 'rcsa', 'rcsb' ) ]
            if ( tf == 'H-NS' ) tf <- b.nums[ c( 'hns' ) ]
            if ( tf == 'MqsR-MqsA' ) tf <- b.nums[ c( 'mqsr', 'mqsa' ) ]
            if ( tf == 'RcsB-BglJ' ) tf <- b.nums[ c( 'rcsb', 'bglj' ) ]
            if ( tf == 'IHF' ) tf <- b.nums[ c( 'ihfa', 'ihfb' ) ]
            if ( all( is.na( tf ) ) ) { cat( tf, as.character(tmp$eco.tf[i]), 'no hit\n' ); next }
            tf <- tf[ ! is.na(tf) ]
        }
        net <- expand.grid( tf, genes )
        network <- rbind( network, data.table( tf=net[,1], target=net[,2], weight=1-tmp$distance[i] ) )
    }
    mot.networks[[ f1 ]] <- network
}

##for (f1 in names(mot.networks)) {
mot.auprs <- list()
for ( f1 in 1:length(rdatas) ) {
    f1 <- basename(rdatas[f1])
    network <- mot.networks[[f1]]
    if ( is.null( network ) ) next
    f1a <- gsub( '.RData', '', f1 )
    if ( file.exists( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) ) ) {
        load( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) )
        mot.auprs[[f1]] <- attr(aupr, 'AUC')
    } else {
        aupr <- get.aupr( network, gold, plot=F, weight.cut=1e-5 )
        save( network, aupr, file=sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) )
    }
    cat( f1, nrow(network), attr(aupr, 'AUC'), attr(aupr,'nPredAt')[1], "\n" )
}

if ( TRUE ) {
## merge all prediction sets into a single weighted network
## for fun, let's compute AUPR stats after we add each network
if ( ! exists( 'do.plot' ) ) do.plot <- FALSE
if ( do.plot ) pdf( sprintf( '%s/motif_cumulative.pdf', output.dir ) )
big.net <- unique(mot.networks[[1]])
aupr <- get.aupr( big.net, gold, plot=do.plot, weight.cut=1e-5 )
summ.stats <- data.table()
for ( i in 2:length(mot.networks) ) {
    aupr <- get.aupr( big.net, gold, plot=do.plot, weight.cut=1e-5 )
    cat( names(mot.networks)[i], nrow(big.net), attr(aupr, 'AUC'), attr(aupr,'nPredAt')[1], "\n" )
    summ.stats <- rbind( summ.stats, data.table( file=names(mot.networks)[i], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'nPredAt')[1] ) )
    tmp <- try( combine.networks( big.net, unique(mot.networks[[i]]) ) )
    if ( ! 'try-error' %in% class(tmp) ) big.net <- unique(tmp)
}
if ( do.plot ) {
    plot( summ.stats$auc )
    dev.off()
}

## let's do it for 12 different randomized orderings of adding the networks
summ.stats2 <- mclapply( 1:12, function(rnd) {
    ord <- sample( 1:length(mot.networks) )
    if ( i == 1 ) ord <- 1:length(mot.networks) ## use one ordering that's the same as the run order
    big.net <- unique(mot.networks[[ord[1]]])
    summ.stats <- data.table()
    for ( i in 2:length(mot.networks) ) {
        aupr <- get.aupr( big.net, gold, plot=F, weight.cut=1e-5 )
        cat( rnd, i, names(mot.networks)[ord[i]], nrow(big.net), attr(aupr, 'AUC'),
            attr(aupr,'nPredAt')[1], "\n" )
        summ.stats <- rbind( summ.stats, data.table( file=names(mot.networks)[ord[i]], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'nPredAt')[1] ) )
        tmp <- try( combine.networks( big.net, unique(mot.networks[[ord[i]]]) ) )
         if ( ! 'try-error' %in% class(tmp) ) big.net <- unique(tmp)
    }
    summ.stats
}, mc.preschedule=F )

summ.stats3 <- sapply(summ.stats2,function(i)if(class(i)!='try-error') i[,auc] else NULL)
if (is.list(summ.stats3)) summ.stats3 <- do.call(cbind, summ.stats3[ ! sapply(summ.stats3, is.null) ] )
if ( do.plot ) {
    pdf( sprintf( '%s/motif_cumulative2.pdf', output.dir ) )
    ##matplot( summ.stats3, typ='l', xlab='# cMonkey runs', ylab='AUPR' )

    all.aucs <- unlist( mot.auprs )
    matplot( summ.stats3, typ='l', xlab='# cMonkey runs', ylab='AUPR', main='EGRIN 2.0 GRE Performance',
            xlim=c(0, 110), ylim=range(c(all.aucs,summ.stats3)) )
    boxplot( all.aucs,add=T, at=108, boxwex=10 )
    text( 104, mean(all.aucs), 'Indivual:', adj=c(1, 0.5) )
    
    boxplot( t( summ.stats3 ), xlab='# cMonkey runs', ylab='AUPR' )

    ## Try this: how significant is the increase between a given iter and the last one (to see if
    ##    it was worth running more runs than that)
    pvals <- sapply( 1:(nrow(summ.stats3)-1), function(i)
                    t.test(summ.stats3[i,],summ.stats3[nrow(summ.stats3),])$p.value )
    plot( -log10(pvals) )
    dev.off()
}

save( summ.stats, summ.stats2, summ.stats3, file=sprintf( '%s/motif_cumulative2.RData', output.dir ) )

rm( tmp, big.net ); gc()
}

}

######################################################################
######################################################################
######################################################################
#########################################################################
## OK, let's do the Inferelator inference on each cMonkey run...
## NOTE THIS HAS ECOLI SPECIFIC CODE, uses the DREAM5 TFs list

if ( 'nwinf' %in% tasks ) {
require( cMonkeyNwInf )
load( 'format_data.RData' ) ## my formatted DREAM5 data -- get the TFs only
rm( ratios, ratios.normed, names2, names1, names, gene.names, gene.lookup, gene.lookup2, gene.lookup.fin, names.syns )

for ( f1 in 1:length(rdatas) ) {
    f1 <- basename(rdatas[f1])
    cat('INFERELATOR', f1, '\n')
    f1a <- gsub( '.RData', '', f1 )
    if ( file.exists( sprintf('%s/%s_nwInf.RData', output.dir, f1a) ) ) {
        next
        #load( sprintf('%s/%s_nwInf.RData', output.dir, f1a) )
        #return( coeffs )
    }
    system( sprintf('touch %s/%s_nwInf.RData', output.dir, f1a) ) ## placeholder
    load( sprintf( '%s/%s', output.dir, f1 ) )
    rm( ratios ); gc()
    coeffs <- try( runnit.wrapper.halo( e, predictors=tfs, cv.choose="min+2se", tf.groups=Inf, alpha=0.8,
                                n.boot=1, tau=0, r.cutoff=Inf, r.filter=Inf, weighted=T, funcs=NA, aic.filter=NA,
                                plot=F ) )
    if ( 'try-error' %in% class( coeffs ) ) { cat( f1, 'ERROR!!!', '\n' ); next }
    save( coeffs, file=sprintf('%s/%s_nwInf.RData', output.dir, f1a) )
    for ( k in 1:length(coeffs) ) coeffs[[k]]$plot.info <- NULL
    save( coeffs, file=sprintf('%s/%s_nwInf_sm.RData', output.dir, f1a) )
    rm(predictor.mats,conditions.excluded,coeffs); gc()
}
}

######################################################################
######################################################################
######################################################################
## Compare Inferelator predictions to RegulonDB
## system('wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt')
## for MSB paper let's use the DREAM5 gold standard network

if ( 'nwinf.regdb.aupr' %in% tasks ) {

##if ( ! exists( 'all.rdata.clusterStacks' ) ) all.rdata.clusterStacks <- get.all.clusterStacks()

##for ( f1 in 1:length(rdatas) ) {
coefs <- lapply( 1:length(rdatas), function(f1) {
    f1 <- basename(rdatas[f1])
    f1a <- gsub( '.RData', '', f1 )
    if ( ! file.exists( sprintf('%s/%s_nwInf_sm.RData', output.dir, f1a) ) ) return( list() )
    if ( file.exists( sprintf('%s/%s_nwInf_sm_vsRegDB.RData', output.dir, f1a) ) ) {
        load( sprintf('%s/%s_nwInf_sm_vsRegDB.RData', output.dir, f1a) )
        cat( f1, nrow(coef), attr(aupr, 'AUC'), attr(aupr,'nPredAt')[1], "\n" )
        return( list( coef=coef, aupr=aupr ) )
    }
    load( sprintf( '%s/%s_clusterStacks.RData', output.dir, f1a ) )
    load( sprintf('%s/%s_nwInf_sm.RData', output.dir, f1a) )
    
    ## Convert the tfs -> bicluster network to a tf -> gene network
    coef <- do.call( rbind, lapply( coeffs, function(i) {
        coe <- c( i$coeffs, i$possibly.regulates ) ## include possibly.regulates!!!
        if ( length( coe ) <= 1 ) return( NULL )
        tmp <- as.data.table( expand.grid( names(coe), clusterStack[[i$k]]$rows ) )
        setnames( tmp, c( 'tf', 'target' ) )
        tmp$weight <- coe[ tmp$tf ]
        tmp } ) )
    coef <- coef[ as.character(tf) != as.character(target) ]
    coef <- coef[ order( abs( coef$weight ), decreasing=T ), ]

    dupes <- duplicated( coef[ ,c('tf', 'target'), with=F ] ) ## add together weights for duplicated interactions
    duped <- coef[ dupes, ]
    coef <- coef[ !dupes, ]
    coef <- combine.networks( coef, duped )
    
    aupr <- get.aupr( coef, gold, plot=F )
    save( coef, aupr, file=sprintf('%s/%s_nwInf_sm_vsRegDB.RData', output.dir, f1a) )
    cat( f1, nrow(coef), attr(aupr, 'AUC'), attr(aupr,'nPredAt')[1], "\n" )
    return( list( coef=coef, aupr=aupr ) )
}, mc.preschedule=F )
names( coefs ) <- rdatas 

if ( TRUE ) {
## merge all prediction sets into a single weighted network
## for fun, let's compute AUPR stats after we add each network
## TBD: add the beta, or abs(beta) ?? For paper, we did just beta...
if ( ! exists( 'do.plot' ) ) do.plot <- FALSE
if ( do.plot ) pdf( sprintf( '%s/nwInf_cumulative.pdf', output.dir ) )
big.net <- unique(coefs[[1]]$coef)
##big.net$weight <- abs(big.net$weight)
aupr <- get.aupr( big.net, gold, plot=do.plot )
summ.stats <- data.table( file=names(coefs)[1], nr=nrow(big.net),
                         auc=attr(aupr, 'AUC'), npred=attr(aupr,'nPredAt')[1] )
for ( i in 2:length(coefs) ) {
    aupr <- get.aupr( big.net, gold, plot=do.plot )
    cat( names(coefs)[i], nrow(big.net), attr(aupr, 'AUC'), attr(aupr,'nPredAt')[1], "\n" )
    summ.stats <- rbind( summ.stats, data.table( file=names(coefs)[i], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'nPredAt')[1] ) )
    net <- unique(coefs[[i]]$coef)
    ##net$weight <- abs(net$weight)
    tmp <- try( combine.networks( big.net, net ) )
    if ( ! 'try-error' %in% class(tmp) ) big.net <- unique(tmp)
}
if ( do.plot ) {
    plot( summ.stats$auc )
    dev.off()
}

## let's do it for 12 different randomized orderings of adding the networks
summ.stats2 <- mclapply( 1:12, function(rnd) {
    ord <- sample( 1:length(coefs) )
    if ( rnd == 1 ) ord <- 1:length(coeffs) ## use one ordering that's the same as the run order
    big.net <- coefs[[ord[1]]]$coef
    summ.stats <- data.table()
    for ( i in 2:length(coefs) ) {
        aupr <- get.aupr( big.net, gold, plot=F )
        cat( rnd, i, names(coefs)[ord[i]], nrow(big.net), attr(aupr, 'AUC'), attr(aupr,'nPredAt')[1], "\n" )
        summ.stats <- rbind( summ.stats, data.table( file=names(coefs)[ord[i]], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'nPredAt')[1] ) )
        tmp <- try( combine.networks( big.net, coefs[[ord[i]]]$coef ) )
        if ( ! 'try-error' %in% class(tmp) ) big.net <- tmp
        rm( tmp ); gc()
    }
    summ.stats
}, mc.preschedule=F )
    
summ.stats3 <- sapply(summ.stats2,function(i)if(class(i)!='try-error') i[,auc] else NULL)
if (is.list(summ.stats3)) summ.stats3 <- do.call(cbind, summ.stats3[ ! sapply(summ.stats3, is.null) ] )
##summ.stats3 <- sapply(summ.stats2,function(i)i[,auc])
if ( do.plot ) {
    pdf( sprintf( '%s/nwInf_cumulative2.pdf', output.dir ) )
    all.aucs <- sapply( coefs, function(i) attr(i$aupr, 'AUC') )
    matplot( summ.stats3, typ='l', xlab='# cMonkey runs', ylab='AUPR', main='EGRIN 2.0 Inferelator Performance',
            xlim=c(0, 110), ylim=range(c(all.aucs,summ.stats3)) )
    boxplot( all.aucs,add=T, at=108, boxwex=10 )
    text( 104, mean(all.aucs), 'Indivual:', adj=c(1, 0.5) )
    boxplot( t( summ.stats3 ), xlab='# cMonkey runs', ylab='Inferelator AUPR' )
    dev.off()
}
save( summ.stats, summ.stats2, summ.stats3, file=sprintf( '%s/nwInf_cumulative2.RData', output.dir ) )

rm( tmp, big.net ); gc()
}
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Get clusterStacks for each run, save individually.
## Compute gene-gene co-reg network from these
if ( 'corems' %in% tasks ) {
    
if ( ! exists( 'all.rdata.clusterStacks' ) ) all.rdata.clusterStacks <- get.all.clusterStacks()

gene.tab <- table( unlist( lapply( all.rdata.clusterStacks, lapply, '[[', 'rows' ) ) )
genes <- sort( names( gene.tab ) )

cond.tab <- table( unlist( lapply( all.rdata.clusterStacks, lapply, '[[', 'cols' ) ) )
conds <- sort( names( cond.tab ) )

gene.gene <- matrix( 0, nrow=length(genes), ncol=length(genes) )
rownames(gene.gene) <- colnames(gene.gene) <- genes

for ( run in all.rdata.clusterStacks ) {
    for ( clust in run ) {
        if ( length( clust$rows ) <= 1 ) next
        rows <- clust$rows ##t( combn( clust$rows, 2 ) )
        gene.gene[ rows, rows ] <- gene.gene[ rows, rows ] + 1
    }
    print(range(gene.gene))
}

stop()

## git clone https://github.com/scalefreegan/Corems.git
## cd Corems; unzip weighted_community_codes.zip
## g++ -o adjmat2wpairs adjmat2wpairs.cpp; g++ -o cluster_communities cluster_communities.cpp
## g++ -o getting_communities getting_communities.cpp
## NEED BOOST LIB FOR compute_tanimoto

R.m <- gene.gene
diag(R.m) <- 0
R.m.2 <- t(apply(R.m,1,function(x) if (sum(x,na.rm=T)>0) return(x/sum(x,na.rm=T)) else return(x/Inf) ))

cwd <- setwd("Corems")
source("analyze_corems.R")
source("processEGRIN.R")
source("processEGRIN.R")
source("analyze_corems.R")
source("main.R")
setwd(cwd)
load("format_data.RData"); ratios=ratios.normed; rm(ratios.normed,names,names1,names2,names.syns); gc() ## for ratios
gBg=R.m.2
rm(R.m,R.m.2); gc()
o <- new.env(parent = baseenv())
o$parameters <- lapply(params,function(i){eval(as.symbol(i))})
names(o$parameters) <- params
system(paste("mkdir",OUTDIR,sep=" "))
system(paste("mkdir filehash"))
gBg.backbone <- multiscaleBackbone(gBg, multicore=F)
writeEdgeList(gBg.backbone)
cwd=setwd(OUTDIR)
system("ln -s ../Corems/adjmat2wpairs; ln -s ../Corems/cluster_communities; ln -s ../Corems/compute_tanimoto; ln -s ../Corems/getting_communities")
runCoremDetection(numGenes = dim(gBg)[1])
o$link.community.threshold <- chooseCutoff()
graphics.off()
setwd(cwd)
o$gg <- list(gBg=gBg, gBg.backbone=gBg.backbone)
rm(gBg, gBg.backbone); gc()
o$corems <- list()
o$corems$all <- loadCorems(o$link.community.threshold)
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Let's change how we do tomtom-ing -- let's tomtom (and save the results of) each
##    output RData file (from above) against each other. That way we can do it during
##    the cMonkey running.
## Note if ON.CLUSTER == TRUE then we are running it on a cluster - do one iteration
##    and then stop.

if ( 'tomtom' %in% tasks ) {

e.value.cutoff=100; resid.cutoff=0.8; dist.meth="ed"; q.thresh=0.5; min.overlap=4; q.pseudo=0; t.pseudo=0
cmd <- "%s/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"
## see http://meme.nbcr.net/meme/doc/tomtom.html for output format of tomtom

for ( f1 in 1:length(rdatas) ) {
##lapply( 1:length(rdatas), function(f1) {
    f1 <- rdatas[f1]
    cat('TOMTOM', f1, '\n')
    load( sprintf( '%s/%s', output.dir, f1 ) )
    e1 <- e; rm( e )
    ##sys.source( "~/scratch/biclust/cmonkey-motif-other.R", envir=e1 ) ## make sure the env has "all.motifs.to.mast.file()"
    e1$all.motifs.to.meme.file <- all.motifs.to.meme.file; environment(e1$all.motifs.to.meme.file) <- e1
    tf1 <- e1$all.motifs.to.meme.file(e.val=Inf)

    f1a <- basename( gsub( '.RData', '', f1 ) )
    
    system( sprintf('rpl bic_ %s_ %s', f1a, tf1) )

    ##for ( f2 in 1:length(rdatas2) ) { ## all vs. all including self vs self -- note should use mclapply to parallelize
    apply.func <- lapply
    if ( ! ON.CLUSTER ) apply.func <- mclapply
    apply.func( 1:length(rdatas2), function(f2) {
        f2 <- rdatas2[f2]
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
        rm( e2 )
        system( sprintf('rpl bic_ %s_ %s', f2a, tf2) )
        if ( file.exists( sprintf('%s.bz2', fname) ) ) {
            print( sprintf('SKIPPING: %s', fname) ); return() } ## someone else is already doing this comparison
        if ( file.exists( sprintf('%s/%s_vs_%s_tomtom.tsv.bz2', output.dir, f2a, f1a) ) ) {
            print( sprintf('SKIPPING (REVERSED): %s', fname) ); return() } ## symmetric -- already did this comparison
        cmd.1 <- sprintf( cmd, e1$progs.dir, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tf2 )
        cmd.1 <- paste( cmd.1, "-query", tf1, '|bzip2 -c >', sprintf('%s.bz2', fname) )
        print( cmd.1 )
        print( sprintf('RUNNING: %s', fname) )
        system( cmd.1 ) ## note this has not been filtered to remove self-self comparisons
        if ( ON.CLUSTER ) stop()

        if ( FALSE ) { ## used to be - read in results and save out processed file.
            tout <- system( cmd.1, intern=T )
            tout <- do.call( rbind, strsplit( tout, "\t" ) )
            colnames( tout ) <- tout[ 1, ,drop=F ]; tout <- tout[ -1, ,drop=F ]
            tout <- as.data.frame( tout[ tout[ ,1 ] != tout[ ,2 ], ,drop=F ] ) ## de-symmetrize the output
            cat( "TOUT:", dim( tout ), "\n" )
            ## Save the file
            write.table( tout, row.names=F, sep='\t', quote=F, file=sprintf('%s/%s_vs_%s_tomtom.tsv', output.dir, f1a, f2a) )
            system( sprintf( 'bzip2 -fv %s/%s_vs_%s_tomtom.tsv', output.dir, f1a, f2a) )
            rm( tout )
        }
    } ) ##, mc.preschedule=F )
    
    rm(e)
}

##rm(list=ls()); gc()
stop()
}

if ( 'sqlite' %in% tasks ) {

## Write out all vital cMonkey info to a set of tsv's; then optionally convert them to a sqlite database file.
cmonkey.run.to.sqlite <- function(e, fname, to.sqlite=F) {
  require( data.table )
  dir.create( fname )
  
  ## PARAMETERS
  tmp <- NULL
  for ( param in ls(e$cmonkey.params) ) {
      print(param)
      val <- e$cmonkey.params[[param]]
      if ( is.null(val) || length(val) <= 0 ) next
      if ( is.matrix(val) ) val <- dim(val)
      if ( is.list(val) && length(val) == 1 ) val <- val[[1]]
      if ( is.numeric(val) ) val <- round(val,dig=4)
      if ( param == 'session.info' ) {
          val['TERMCAP'] <- ''
          val <- val[nchar(val) > 0]
          val <- paste(names(val), val, sep='=', collapse=',')
      } else {
          if ( length(val) > 1 ) val <- paste(unlist(val),collapse=',')
      }
      val <- gsub( '[\n\t]', '', val )
      tmp <- rbind( tmp, data.table( param=param, val=as.character(val) ) )
  }
  write.table( tmp, file=sprintf('%s/%s',fname,'params.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'params.tsv'), wait=F )

  ## STATS
  stats <- e$stats
  for ( i in colnames(stats) ) if ( is.numeric(stats[[i]]) ) stats[[i]] <- round(stats[[i]], dig=4)
  write.table( stats, file=sprintf('%s/%s',fname,'stats.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
  system( sprintf( 'bzip2 -fv %s/%s',fname,'stats.tsv'), wait=F )
  rm( stats )
  
  ## DATA - ratios, string, operons, annotations
  #tmp <- data.table( round(e$ratios$ratios,dig=3) )
  #tmp <- data.table( gene=rownames(e$ratios$ratios), tmp )
  tmp <- which( e$ratios$ratios < Inf, arr=T )
  tmp <- data.table( gene=rownames(e$ratios$ratios)[ tmp[,1] ], cond=colnames(e$ratios$ratios)[ tmp[,2] ],
                    ratio=(round(e$ratios$ratios,dig=3))[tmp] )
  write.table( tmp, file=sprintf('%s/%s',fname,'ratios.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'ratios.tsv'), wait=F )
  tmp <- data.table( e$genome.info$feature.tab )
  write.table( tmp, file=sprintf('%s/%s',fname,'feature_tab.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'feature_tab.tsv'), wait=F )
  tmp <- data.table( e$genome.info$feature.names )
  setcolorder( tmp, c(2,1,3) )
  write.table( tmp, file=sprintf('%s/%s',fname,'feature_names.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'feature_names.tsv'), wait=F )
  tmp <- data.table( e$genome.info$operons )
  setcolorder( tmp, c(2,1) )
  write.table( tmp, file=sprintf('%s/%s',fname,'operons.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'operons.tsv'), wait=F )
  tmp <- data.table( gene=names(e$genome.info$all.upstream.seqs[[1]]), sequence=e$genome.info$all.upstream.seqs[[1]] )
  write.table( tmp, file=sprintf('%s/%s',fname,'upstream_seqs.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'upstream_seqs.tsv'), wait=F )
  if ( ! is.null( e$networks$operons ) ) {
    tmp <- data.table( e$networks$operons )
    write.table( tmp, file=sprintf('%s/%s',fname,'operon_network.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
    rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'operon_network.tsv'), wait=F )
  }
  if ( ! is.null( e$networks$string ) ) {
    tmp <- data.table( e$networks$string )
    write.table( tmp, file=sprintf('%s/%s',fname,'string_network.tsv'), quote=F, sep='\t', row.names=F, col.names=T, append=F )
    rm( tmp ); system( sprintf( 'bzip2 -fv %s/%s',fname,'string_network.tsv'), wait=F )
  }
  
  ## BICLUSTERS - info, genes, conditions, motifs
  bzcon1 <- sprintf('%s/bicluster.tsv',fname)
  bzcon2 <- sprintf('%s/bicluster_genes.tsv',fname)
  bzcon3 <- sprintf('%s/bicluster_conds.tsv',fname)
  ##bzcon4 <- sprintf('%s/bicluster_motifs.tsv',fname)
  wrote <- FALSE
  for ( k in 1:e$k.clust ) {
    ##bb <- paste('BIC',k,sep='_')
    if ( k %% 100 == 0 ) print(k)
    tab1 <- tab2 <- tab3 <- tab4 <- NULL
    clust <- e$clusterStack[[k]] ##out$get.bicluster.info(bb)[[1]]
    if ( ! is.null(clust) ) {
      ##fname <- names(e$fnames.to.cluster[which(e$fnames.to.cluster==k)])
      ##if ( is.null( fname ) ) fname <- ''
      ##k.orig <- clust$k
      tab1 <- data.table( bic=k, nrow=clust$nrow, ncol=clust$ncol, resid=round(clust$resid,dig=3),
                         pclust=round(clust$p.clust,dig=3), eval=min(clust$e.val,na.rm=T) ) ##, k_orig=k.orig, fname=fname )
      if ( ! is.null( tab1 ) ) write.table( tab1, bzcon1, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    }
    genes <- clust$rows
    if ( ! is.null( genes ) ) tab2 <- data.table( bic=k, gene=unique(genes) )
    if ( ! is.null( tab2 ) ) write.table( tab2, bzcon2, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    conds <- clust$cols
    if ( ! is.null( conds ) ) tab3 <- data.table( bic=k, cond=unique(conds) )
    if ( ! is.null( tab3 ) ) write.table( tab3, bzcon3, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    ##mots <- out$get.motifs(bicluster=bb)[[1]]
    ##if ( ! is.null( mots ) ) tab4 <- data.table( bic=k, mot=gsub('MOT_','',unique(mots)) )
    ##if ( ! is.null( tab4 ) ) write.table( tab4, bzcon4, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    wrote <- TRUE
  }
  system( sprintf( 'bzip2 -fv -9 %s', bzcon1 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon2 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon3 ), wait=F )
  ##system( sprintf( 'bzip2 -fv -9 %s', bzcon4 ), wait=F )

  ## MOTIFS - info; meme positions; mast positions; pssms
  bzcon1 <- sprintf('%s/motif.tsv',fname)
  bzcon2 <- sprintf('%s/motif_meme_posn.tsv',fname)
  bzcon3 <- sprintf('%s/motif_mast_posn.tsv',fname)
  bzcon4 <- sprintf('%s/motif_pssm.tsv',fname)
  wrote <- FALSE
  for ( k in 1:e$k.clust ) {
    ##bic <- paste('BIC',k,sep='_')
    mots <- e$meme.scores[[ 1 ]][[ k ]]$meme.out ##out$get.motifs(bicluster=bic)[[1]]
    pv <- e$meme.scores[[ 1 ]][[ k ]]$pv.ev[[2]]
    for ( mm in 1:length(mots) ) {
      if ( k %% 100 == 0 ) print(paste(k,mm))
      m <- mots[[ mm ]]
      minf <- m ##out$get.motif.info(m)[[1]]
      if ( ! is.null( minf ) ) {
        ##cf <- out$coding.fracs$all.fracs[m]
        tab <- data.table( bic=k, mot=mm, width=minf$width, llr=minf$llr, eval=minf$e.value, sites=minf$sites ) ##,
                          ##coding=cf, good=as.integer((cf < out$coding.fracs$mean.fracs - 0.01)) )
        if ( ! is.null( tab ) ) write.table( tab, bzcon1, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        tab <- NULL; if ( nrow(minf$posns) > 0 ) tab <- as.data.table( cbind( bic=k, mot=mm, minf$posns ) )
        if ( ! is.null( tab ) ) write.table( tab, bzcon2, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        tab <- NULL; if ( ! is.null(pv) && nrow(pv) > 0 ) {
            tmp <- subset( pv, abs(mots) == mm )
            if ( nrow(tmp) > 0 ) tab <- as.data.table( cbind( bic=k, mot=mm, tmp ) )
        }
        if ( ! is.null( tab ) ) write.table( tab, bzcon3, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        pssm <- minf$pssm; colnames( pssm ) <- e$col.let; pssm <- as.data.table( pssm )
        tab <- data.table( cbind( bic=k, mot=mm, ind=1:nrow(pssm), round(pssm,dig=3) ) )
        if ( ! is.null( tab ) ) write.table( tab, bzcon4, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        wrote <- TRUE
      }
    }
  }
  system( sprintf( 'bzip2 -fv -9 %s', bzcon1 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon2 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon3 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon4 ), wait=F )


  ## INFERELATOR COEFFS (TBD)


  ## TRICK TO READ THESE IN FAST (using data.table package):
  ## fread(paste(readLines(bzfile('bicluster.tsv.bz2')),collapse='\n')) ## never mind - this crashes. Need to bunzip2 it first, then use fread()
  ## or use sqldf
  ## read.csv.sql(bzfile('filehash/DATABASES/bicluster.tsv.bz2'),sep='\t')
  ## so this:
  ## tf <- e$my.tempfile(); system(sprintf('bunzip2 -c filehash/DATABASES/motif_promoter_counts.tsv.bz2 >%s', tf))
  ## tab <- fread(tf); unlink(tf)
  ## or, e.g.:
  ## tf <- file(tf); tmp <- sqldf("select ind,count from tf where motc=0",file.format=list(sep='\t',head=T),verbose=T)
  ## or
  ## df=read.csv.sql(filter='bunzip2 -c filehash/DATABASES/motif.tsv.bz2')

  if ( to.sqlite ) {
  ## once these are created, can do e.g.:
  ## sqldf('select * from motif limit 10', dbname='filehash/DATABASES/motif.db')
  ## or:
  ## sqldf('select * from motif_scans where bic=10000', dbname='filehash/DATABASES/motif_scans.db')
  ## or:
  ## sqldf(dbname='filehash/DATABASES/motif_scans.db')  ## open a persistent connection
  ##    sqldf('create index pval_index on motif_scans(p_value)')
  ##    tmp=sqldf('select * from motif_scans where p_value<1e-7')
  ## sqldf() ## close the connection
  files <- list.files(path=fname, patt=glob2rx('*.tsv.bz2'), full=T )
  unlink( sprintf( '%s/%s.db', fname, 'cmonkey' ) )
  for ( f in files ) {
      print(f)
      typematch <- c( integer='integer', character='character', numeric='real' )
      tf <- tempfile(); system(sprintf('bunzip2 -c %s | head -1000000 >%s', f, tf))
      tmp <- fread(tf, nrows=1000000, sep='\t')
      print(head(tmp,5))
      classes <- sapply( tmp, class )
      if ( any( classes == 'character' ) ) {
          max.char <- sapply( tmp, function(i) max(nchar(i)) )[which(classes == 'character')]
          for ( i in max.char ) typematch[sprintf('character(%d)',i)] <- sprintf('character(%d)',i)
          classes[which(classes == 'character')] <- sprintf('character(%d)',max.char)
          classes[which(max.char>1024)] <- 'blob'
      }
      rm(tmp)
      print(classes)
      system(sprintf('bunzip2 -cv %s | tail -n +2 >%s', f, tf)) ## skip first (header) line
      system(sprintf('wc -l %s',tf))
      table.name <- gsub('.tsv.bz2','',basename(f),fixed=T)
      names(classes) <- gsub( '[.+-/|]', '_', names(classes) )
      index.classes <- classes ##classes[classes != 'numeric']
      scrp <- tempfile()
      ## cat( sprintf( 'create table %s (%s, primary key(%s));\n', table.name,
      ##              paste(paste(names(classes), typematch[classes]), collapse=','),
      ##              paste(names(index.classes),collapse=',') ), file=scrp, append=F )
      ## NOTE primary key not needed. Sqlite automatically generates a "rowid" column used as primary key
      cat( sprintf( 'create table %s (%s);\n', table.name,
                   paste(paste(names(classes), typematch[classes]), collapse=',') ), file=scrp, append=F )
      cat( '.separator "\\t"\n', file=scrp, append=T )
      cat( sprintf( '.import %s %s\n', tf, table.name ), file=scrp, append=T )
      cat( sprintf( 'create index %s_index on %s(%s);\n', table.name, table.name,
                   paste(names(index.classes),collapse=',') ), file=scrp, append=T )
      cat( '.quit\n', file=scrp, append=T )
      ##unlink( sprintf( '%s/%s.db', fname, table.name ) )
      system( sprintf( 'cat %s', scrp ) )
      system( sprintf( 'sqlite3 %s/cmonkey.db <%s', fname, scrp ) ) ## put them all in a single cmonkey.db db file
      unlink(tf); unlink(scrp)
  }
  ## cat( paste( 'Script to attach all databases:\n',
  ##              paste( sprintf( "attach database '%s/%s.db' as %s;\n", fname,
  ##                             gsub('.tsv.bz2','',basename(files),fixed=T),
  ##                             gsub('.tsv.bz2','',basename(files),fixed=T) ),
  ##                    collapse='' ) ) )
}
}

for ( f1 in rdatas ) {
    f1 <- gsub( '.RData', '', basename(f1) )
    print(f1)
    if ( file.exists( sprintf('%s/%s/cmonkey.db', output.dir, f1) ) ) next
    load( sprintf( '%s/%s.RData', output.dir, f1 ) )
    cmonkey.run.to.sqlite( e, sprintf('%s/%s', output.dir, f1), to.sqlite=T )
    rm( e ); gc()
}

}
