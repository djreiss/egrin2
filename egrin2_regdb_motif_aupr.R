
max.motifs <- Inf
hit.p.cutoff <- 1e-4
distance.cutoff <- Inf ##0.99

## Need to compute conversion from AaaX to b### number.
## try the table gotten from here: http://www.ecogene.org/?q=ecodownload/dbtable
## needed a bit of editing to get it to load -- replace the 'b#' with 'bnum' and the last tab in the first line

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
            } ##else {
                ##cat( "HERE2", n, nrow(tmp), '\n' )
            ##}
        } ##else {
          ##  cat( "HERE1", n, nrow(tmp), '\n' )
        ##}
    }
    gene.lookup
}

##mot.networks <- list()
##for ( f1 in 1:length(rdatas) ) {
mot.networks <- mclapply( 1:length(rdatas), function(f1) {
    f1 <- basename(rdatas[f1])
    print(f1)
    f1a <- gsub( '.RData', '', f1 )
    if ( file.exists( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) ) ) { ## already done
        load( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) )
        return( network )
    } 
    if ( ! file.exists( sprintf('%s/%s_with_regdb_comparison.RData', output.dir, f1a) ) ) return(NULL)
    load( sprintf('%s/%s_with_regdb_comparison.RData', output.dir, f1a) )

    ms <- get.meme.outs( f1 )
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

    b.nums <- get.b.numbers(as.character(tmp$eco.tf))
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
    if ( is.null( network ) ) return( NULL )
    aupr <- get.aupr( network, gold, plot=F, weight.cut=1e-5 )
    cat( f1, nrow(network), attr(aupr, 'AUC'), attr(aupr,'npred_at_prec25'), "\n" )
    save( network, aupr, file=sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) )
    ##mot.networks[[ f1 ]] <- network
    network
}, mc.preschedule=F )
names( mot.networks ) <- rdatas

if ( FALSE ) {
## merge all prediction sets into a single weighted network
## for fun, let's compute AUPR stats after we add each network
## mot.networks <- list()
## for ( f1 in 1:length(rdatas) ) {
##     f1 <- basename(rdatas[f1])
##     f1a <- gsub( '.RData', '', f1 )
##     if ( ! file.exists( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) ) ) next
##     load( sprintf('%s/%s_motifs_vsRegDB.RData', output.dir, f1a) )
##     mot.networks[[f1]] <- network
##     cat( f1, nrow(network), attr(aupr, 'AUC'), attr(aupr,'npred_at_prec25'), "\n" )
## }

pdf( 'motif_cumulative.pdf' )
big.net <- mot.networks[[1]]
summ.stats <- data.table()
for ( i in 2:length(mot.networks) ) {
    aupr <- get.aupr( big.net, gold, plot=T, weight.cut=1e-5 )
    cat( names(mot.networks)[i], nrow(big.net), attr(aupr, 'AUC'), attr(aupr,'npred_at_prec25'), "\n" )
    summ.stats <- rbind( summ.stats, data.table( file=names(mot.networks)[i], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'npred_at_prec25') ) )
    tmp <- mot.networks[[i]] ##; tmp$weight <- tmp$weight / max(tmp$weight,na.rm=T)
    tmp <- try( combine.networks( big.net, tmp ) )
    if ( ! 'try-error' %in% class(tmp) ) big.net <- tmp
}
plot( summ.stats$auc )
dev.off()
}

if ( FALSE ) {
## let's do it for 12 different randomized orderings of adding the networks
summ.stats2 <- mclapply( 1:30, function(rnd) {
    ord <- sample( 1:length(mot.networks) )
    big.net <- mot.networks[[ord[1]]]
    summ.stats <- data.table()
    for ( i in 2:length(mot.networks) ) {
        aupr <- get.aupr( big.net, gold, plot=F, weight.cut=1e-5 )
        cat( rnd, i, names(mot.networks)[ord[i]], nrow(big.net), attr(aupr, 'AUC'),
            attr(aupr,'npred_at_prec25'), "\n" )
        summ.stats <- rbind( summ.stats, data.table( file=names(mot.networks)[ord[i]], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'npred_at_prec25') ) )
        tmp <- mot.networks[[ord[i]]] ##; tmp$weight <- tmp$weight / max(tmp$weight,na.rm=T)
        tmp <- try( combine.networks( big.net, tmp ) )
        if ( ! 'try-error' %in% class(tmp) ) big.net <- tmp
    }
    summ.stats
}, mc.preschedule=F )

summ.stats3 <- sapply(summ.stats2,function(i)if(class(i)!='try-error') i[,auc] else NULL)
if (is.list(summ.stats3)) summ.stats3 <- do.call(cbind, summ.stats3[ ! sapply(summ.stats3, is.null) ] )
pdf( 'motif_cumulative2.pdf' )
matplot( summ.stats3, typ='l', xlab='# cMonkey runs', ylab='AUPR' )
boxplot( t( summ.stats3 ), xlab='# cMonkey runs', ylab='AUPR' )
dev.off()
save( summ.stats2, file='motif_cumulative2.RData' )
}

