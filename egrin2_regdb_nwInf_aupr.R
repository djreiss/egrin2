######################################################################
######################################################################
######################################################################
## Compare Inferelator predictions to RegulonDB
## system('wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt')
## for MSB paper let's use the DREAM5 gold standard network

##for ( f1 in 1:length(rdatas) ) {
coefs <- mclapply( 1:length(rdatas), function(f1) {
    f1 <- basename(rdatas[f1])
    f1a <- gsub( '.RData', '', f1 )
    if ( ! file.exists( sprintf('%s/%s_nwInf_sm.RData', output.dir, f1a) ) ) return( list() )
    if ( file.exists( sprintf('%s/%s_nwInf_sm_vsRegDB.RData', output.dir, f1a) ) ) {
        load( sprintf('%s/%s_nwInf_sm_vsRegDB.RData', output.dir, f1a) )
        cat( f1, nrow(coef), attr(aupr, 'AUC'), attr(aupr,'npred_at_prec25'), "\n" )
        return( list( coef=coef, aupr=aupr ) )
    }
    ##load( sprintf( '%s/%s_clusterStacks.RData', output.dir, f1a ) )
    clusterStack <- get.clusterStack( f1 )
    load( sprintf('%s/%s_nwInf_sm.RData', output.dir, f1a) )
    
    ## Convert the tfs -> bicluster network to a tf -> gene network
    coef <- do.call( rbind, lapply( coeffs, function(i) {
        coe <- c( i$coeffs, i$possibly.regulates ) ## include possibly.regulates!!!
        tmp <- as.data.table( expand.grid( names(coe), clusterStack[[i$k]]$rows ) )
        colnames( tmp ) <- c( 'tf', 'target' )
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
    cat( f1, nrow(coef), attr(aupr, 'AUC'), attr(aupr,'npred_at_prec25'), "\n" )
    return( list( coef=coef, aupr=aupr ) )
}, mc.preschedule=F )
names( coefs ) <- rdatas 

if ( FALSE ) {
## merge all prediction sets into a single weighted network
## for fun, let's compute AUPR stats after we add each network
pdf( 'nwInf_cumulative.pdf' )
big.net <- coefs[[1]]$coef
summ.stats <- data.table()
for ( i in 2:length(coefs) ) {
    aupr <- get.aupr( big.net, gold, plot=T )
    cat( names(coefs)[i], nrow(big.net), attr(aupr, 'AUC'), attr(aupr,'npred_at_prec25'), "\n" )
    summ.stats <- rbind( summ.stats, data.table( file=names(coefs)[i], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'npred_at_prec25') ) )
    tmp <- try( combine.networks( big.net, coefs[[i]]$coef ) )
    if ( ! 'try-error' %in% class(tmp) ) big.net <- tmp
}
plot( summ.stats$auc )
dev.off()

## let's do it for 12 different randomized orderings of adding the networks
summ.stats2 <- mclapply( 1:30, function(rnd) {
    ord <- sample( 1:length(coefs) )
    big.net <- coefs[[ord[1]]]$coef
    summ.stats <- data.table()
    for ( i in 2:length(coefs) ) {
        aupr <- get.aupr( big.net, gold, plot=F )
        cat( rnd, i, names(coefs)[ord[i]], nrow(big.net), attr(aupr, 'AUC'), attr(aupr,'npred_at_prec25'), "\n" )
        summ.stats <- rbind( summ.stats, data.table( file=names(coefs)[ord[i]], nr=nrow(big.net),
                                                auc=attr(aupr, 'AUC'), npred=attr(aupr,'npred_at_prec25') ) )
        tmp <- try( combine.networks( big.net, coefs[[ord[i]]]$coef ) )
        if ( ! 'try-error' %in% class(tmp) ) big.net <- tmp
    }
    summ.stats
}, mc.preschedule=F )

summ.stats3 <- sapply(summ.stats2,function(i)i[,auc])
pdf( 'nwInf_cumulative2.pdf' )
matplot( summ.stats3, typ='l', xlab='# cMonkey runs', ylab='AUPR' )
boxplot( t( summ.stats3 ), xlab='# cMonkey runs', ylab='AUPR' )
dev.off()
save( summ.stats2, file='nwInf_cumulative2.RData' )
}
