######################################################################
######################################################################
######################################################################
## The stuff below is required for AUPR functions (nwInf AND motif)

## Updated AUPR function based on Sriram's matlab code
## Assumes input is tf, target, weight and gold is tf, target
## Only uses Inf. predictions that have abs(weight) >= 0.1 ... see ../FINAL_compare_regulondb.R
get.aupr <- function( net, gold, weight.cut=0.1, plot.it=T, ... ) {
    net <- net[ rev( rank( abs(net$weight), ties='random' ) ), ]
    net <- net[ abs(weight) >= weight.cut, ]
    if ( nrow( net ) > 100000 ) net <- net[ 1:100000, ]

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
  if ( plot.it ) plot(Y,X,typ='l',lwd=2, xlab='Recall', ylab='Precision', ...) # plotting it
  attr( resultss, 'AUC' ) <- AUC
  attr( resultss, 'n' ) <- length(netwo)

  ## get # correct predictions at precision of 25%
  tmp <- approx(resultss[,'precision'],resultss[,'recall'],0.25)$y * length(allnet)
  attr( resultss, 'npred_at_prec25' ) <- tmp
  
  resultss
}

combine.networks <- function( net1, net2 ) { ## add the weights if the edges are teh same
    tmp1 <- net1$weight
    names(tmp1) <- paste( net1$tf, net1$target, sep='//' )

    tmp2 <- net2$weight
    names(tmp2) <- paste( net2$tf, net2$target, sep='//' )

    ## Add together the weights???
    tmp3 <- tmp1[ names(tmp1)[names(tmp1) %in% names(tmp2)] ] + tmp2[ names(tmp1)[names(tmp1) %in% names(tmp2)] ]
    tmp3 <- c( tmp3, tmp1[ ! names(tmp1) %in% names(tmp2) ], tmp2[ ! names(tmp2) %in% names(tmp1) ] )

    out <- do.call( rbind, strsplit( names(tmp3), '//' ) )
    out <- data.table( tf=out[,1], target=out[,2], weight=tmp3 )
    out <- out[ order( abs( out$weight ), decreasing=T ), ]
    out
}

gold <- fread( 'DREAM5_NetworkInference_Evaluation/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network3.tsv', head=F )
load( 'format_data.RData' )
rm( ratios, ratios.normed, names2, names1, names, gene.names, gene.lookup, gene.lookup2, names.syns )
gold$tf <- gene.lookup.fin[ as.character( gold$V1 ) ]
gold$target <- gene.lookup.fin[ as.character( gold$V2 ) ]
gold <- subset( gold, ! is.na( tf ) & ! is.na( target ) ) ## only removes 2 edges
require(caTools) ## for trapz
