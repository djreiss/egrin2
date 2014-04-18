######################################################################
######################################################################
######################################################################
#########################################################################
## Get clusterStacks for each run, save individually.
## Compute gene-gene co-reg network from these
    
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
