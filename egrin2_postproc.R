debug.on()
require( parallel )
options( mc.cores=12 )
require( data.table )

## Let's do all of the tomtom-ing, and fimo-ing, etc. on an individual RData basis, rather than after
## integrating all RData files (from different cMonkey runs) into a single big cMonkey run...

if ( ! exists( 'tasks' ) )
    tasks <- c( 'fimo', 'coding.fracs', 'motif.shadows',
               'regdb', 'regdb.aupr',  ## compares MOTIFS to regdb; compute auprs
               'nwinf', 'nwinf.regdb',  ## run nwInf; compares Inf. to regdb; compute auprs
               'motif.sims', 'newmotcluster',
               'corems', 'tomtom' )
if ( ! exists( 'output.dir' ) ) output.dir <- 'output/'

if ( ! exists( 'rdatas' ) ) { ## allows us to just run it on one file
    rdatas <- list.files( path=output.dir, patt=glob2rx('zzz_eco_???.RData'), full=F )
}

source( 'egrin2/egrin2_funcs.R' )

if ( length(tasks) <= 0 ) stop()

######################################################################
######################################################################
######################################################################
#########################################################################
## First fimo all motifs from each run against the genome and save the output...

if ( 'fimo' %in% tasks ) {
    source( 'egrin2/egrin2_fimo.R' )
}

######################################################################
######################################################################
######################################################################
#################################################################
## Compute per-run fraction of times each motif is in a coding region, save the results

if ( 'coding.fracs' %in% tasks ) {
    source( 'egrin2/egrin2_coding_fracs.R' )
}    

######################################################################
######################################################################
######################################################################
#########################################################################
## Compute the "motif shadows" for each run -- locations covered by each motif in each run -
##   save as motif_shadows.RData

if ( 'motif.shadows' %in% tasks ) {
    source( 'egrin2/egrin2_motif_shadows.R' )
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Next compute the overlap (similarities) of all "motif shadows" vs all others - save as motif_sims.tsv.bz2
## This takes a while as it is all runs vs. all runs

if ( 'motif.sims' %in% tasks ) {
    source( 'egrin2/egrin2_motif_sims.R' )
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
    source( 'egrin2/egrin2_newmotcluster.R' )
}    

######################################################################
######################################################################
######################################################################
#########################################################################
## Compute combined pssms for each motif cluster. Here we need to run tomtom
##  between all of the motifs in each motif cluster.
## It will also plot them if an individual cmonkey run "e" is in the environment

if ( 'newmotcluster.getcombinedpssms' %in% tasks && exists( 'clusts' ) ) {
    source( 'egrin2/egrin2_newmotcluster_getcombinedpssms.R' )    
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Now here is the comparison vs. regulondb part -- this compares motif locations with
##   binding sites in regulondb. Do it individually for each run and save the results
## NOTE THIS IS SPECIFIC TO E. COLI ONLY !!!

if ( 'regdb' %in% tasks ) {
    source( 'egrin2/egrin2_regdb_motif.R' )    
    n.mots.hit <- as.integer( sapply(strsplit(sapply(regdb.out,'[',12),' '),'[',13) )
    hist( n.mots.hit, breaks=20 )
}

######################################################################
######################################################################
######################################################################
## Compute the motif-based networks and get AUPR measures

if ( 'regdb.aupr' %in% tasks ) {
    source( "egrin2/egrin2_aupr_funcs.R" )
    
    source( "egrin2/egrin2_regdb_motif_aupr.R" )
}

######################################################################
######################################################################
######################################################################
#########################################################################
## OK, let's do the Inferelator inference on each cMonkey run...
## NOTE THIS HAS ECOLI SPECIFIC CODE, uses the DREAM5 TFs list

if ( 'nwinf' %in% tasks ) {
    source( "egrin2/egrin2_nwInf.R" )
}

######################################################################
######################################################################
######################################################################
## Compare Inferelator predictions to RegulonDB
## system('wget http://regulondb.ccg.unam.mx/menu/download/datasets/files/network_tf_gene.txt')
## for MSB paper let's use the DREAM5 gold standard network

if ( 'nwinf.regdb' %in% tasks ) {
    source( "egrin2/egrin2_aupr_funcs.R" )
    
    source( "egrin2/egrin2_regdb_nwInf_aupr.R" )
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Get clusterStacks for each run, save individually.
## Compute gene-gene co-reg network from these

if ( 'corems' %in% tasks ) {
    source( 'egrin2/egrin2_corems.R' )
}

######################################################################
######################################################################
######################################################################
#########################################################################
## Let's change how we do tomtom-ing -- let's tomtom (and save the results of) each
##    output RData file (from above) against each other. That way we can do it during
##    the cMonkey running.

if ( 'tomtom' %in% tasks ) {
    source( 'egrin2/egrin2_tomtom.R' )
}
