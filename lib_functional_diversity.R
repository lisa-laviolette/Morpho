#' Compute functional diversity indices
#'
#' @param coord morphological space (morphs x axes)
#' @param weight weight matrix (date x morphs)
#' @param verb print information about the progress of computation
multidimFD4  <-  function(coord, weight, verb=TRUE) {

  # library required for indices computation
  library("geometry")
  library("ape")
  # library("foreach")
  # library("doParallel")

  # registerDoParallel(cores=4)

  # saving name of current working directory
  # current_wd <- getwd()


  ##############################################################################
  # info on study case


###############################

  # number and names of axes
  nm_axes <- colnames(coord)
  nb_axes <- length(nm_axes)

  # number and names of assemblages
  nm_asb <- row.names(weight)
  nb_asb <- length(nm_asb)

  # matrix to store results
  indices <- c( "Nb_sp", "Tot_weight", c("FRic","FDiv","FEve") )
  FD <- matrix(NA, nb_asb, length(indices), dimnames=list(nm_asb,indices))

  ##############################################################################
  # preliminary computation at the species pool level

  #######################################

  # convex hull volume of the species pool
  FRic_pool <- convhulln(coord,"FA")$vol

  ##############################################################################

  for (k in nm_asb) {

    # FD[k, "assemblage"] <- k

    ###########################################################
    # preparing data

    # names, number, weight and coordinates of of species present
    weight_k <- weight[k,]
    nm_sp_k <- row.names(coord)[which(weight_k>0)]
    nb_sp_k <- length(nm_sp_k)
    weight_sp_k <- weight[k,nm_sp_k]
    coord_sp_k <- coord[nm_sp_k,]
    if(nb_sp_k==1) { coord_sp_k <- matrix(coord_sp_k,nrow=1,dimnames=list(nm_sp_k,colnames(coord)) ) } # matrix object

    # names of species absent
    nm_sp_absent_k <- names(which(weight[k,]==0))

    # total weight
    FD[k, "Tot_weight"] <- sum(weight_sp_k)

    #relative weight
    rel_weight_sp_k <- weight_sp_k/sum(weight_sp_k)

    # species richness
    FD[k, "Nb_sp"] <- nb_sp_k

    ###########################################################
    # multivariate indices
    # indices based on vertices and volume of convex hull, only if more species than number of axes
    ########################
    # Functional richness = convex hull volume
    FD[k, "FRic"] <- round(convhulln(coord_sp_k,"FA")$vol/FRic_pool,6)

    ########################
    # Functional Divergence

    # identity of vertices
    vert0 <- convhulln(coord_sp_k,"Fx TO 'vert.txt'")
    vert1 <- scan("vert.txt",quiet=T)
    vertices_k <- (vert1+1)[-1]

    # coordinates of the center of gravity of the vertices (B)
    B_k <- apply(coord_sp_k[vertices_k,],2,mean)

    # Euclidean dstance to B (dB)
    dB_k <- apply(coord_sp_k, 1, function(x) { (sum((x-B_k)^2) )^0.5} )

    # mean of dB values and deviations to mean
    meandB_k <- mean(dB_k)
    devdB_k <- dB_k-meandB_k

    # abundance-weighted mean deviation
    abdev_k <-  rel_weight_sp_k*devdB_k
    ababsdev_k <-  rel_weight_sp_k*abs(devdB_k)

    FD[k, "FDiv"] <- round( (sum(abdev_k)+meandB_k) / (sum(ababsdev_k)+meandB_k) ,6)


    # end of if more species than number of axes

    ##########################
    # Functional Evenness

    # inter-species Euclidean distance
    distT_k <- dist(coord_sp_k, method="euclidian")

    # topology of Minimum Spanning Tree and conversion of the 'mst' matrix into 'dist' class
    linkmst_k <- mst(distT_k)
    mstvect_k <- as.dist(linkmst_k)

    # pairwise cumulative relative abundances and conversion into 'dist' class
    ab2_k <- matrix(0,nrow=nb_sp_k,ncol=nb_sp_k)
    for (q in 1:nb_sp_k)
      for (r in 1:nb_sp_k)
        ab2_k[q,r] <- rel_weight_sp_k[q]+rel_weight_sp_k[r] # end of q,r
    ab2vect_k <- as.dist(ab2_k)

    # EW index for the (S-1) segments
    EW_k <- rep(0,nb_sp_k-1)
    flag <- 1
    for (m in 1:((nb_sp_k-1)*nb_sp_k/2))
    {if (mstvect_k[m]!=0) {EW_k[flag] <- distT_k[m]/(ab2vect_k[m]) ; flag <- flag+1}}  # end of m

    # PEW index and comparison with 1/S-1
    minPEW_k <- rep(0,nb_sp_k-1)  ;  OdSmO_k <- 1/(nb_sp_k-1)
    for (l in 1:(nb_sp_k-1))
      minPEW_k[l] <- min( (EW_k[l]/sum(EW_k)) , OdSmO_k )  # end of l

    # FEve
    FD[k, "FEve"] <- round( ( (sum(minPEW_k))- OdSmO_k) / (1-OdSmO_k ) ,6)


    # printing step achieved
    if (verb==TRUE) print(paste("FD of assemblage '",k,"' computed",sep="") )
  }# end of working on assemblage k

  # loop on assemblages for computing and plotting functional diversity indices


  # End of indices computation
  ######################################################################################################################

  ######################################################################################################################

  # returning to current working directory
  # setwd(current_wd)

  # returning results
  return(FD)

}# end of function multidimFD
########################################################################################################################################
########################################################################################################################################
