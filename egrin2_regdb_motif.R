######################################################################
######################################################################
######################################################################
#########################################################################
## Now here is the comparison vs. regulondb part -- this compares motif locations with
##   binding sites in regulondb. Do it individually for each run and save the results
## NOTE THIS IS SPECIFIC TO E. COLI ONLY !!!

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

    load( sprintf('%s/%s_motif_shadows.RData', output.dir, f1a) )
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
