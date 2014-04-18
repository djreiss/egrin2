######################################################################
######################################################################
######################################################################
#########################################################################
## OK, let's do the Inferelator inference on each cMonkey run...
## NOTE THIS HAS ECOLI SPECIFIC CODE, uses the DREAM5 TFs list

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
