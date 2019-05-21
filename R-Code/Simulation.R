simulateCellTypes = function(
    referenceCellTypes,
    markerGenes,
    numSamples,
    seed = NULL,
    verbose = TRUE
){

    # Split reference into marker genes and non marker genes
    if(verbose){message("Processing input files ...")}
    markerCounts = referenceCellTypes[
        rownames(referenceCellTypes) %in% markerGenes,]
    
    if(nrow(markerCounts) == 0){
        stop("Non of the marker genes present in reference cell-types!")
    }
    otherCounts = referenceCellTypes[
        !rownames(referenceCellTypes) %in% markerGenes,]

    # Estimate parameters for marker genes and other genes separately
    if(verbose){message("Estimating parameters ...")}

    paramMarkers = estimateParameters(
        markerCounts,
        simMarkers = TRUE
    )
    paramOthers = estimateParameters(
        otherCounts,
        simMarkers = FALSE
    )

    # Simulate counts
    if(verbose){message("Simulating counts ...")}
    if(is.null(seed)) {
        seed = sample(1:1000000, size = 1)
    }
    #set.seed(seed)

    sim_Counts_Markers = simulateCounts(
        simParam = paramMarkers,
        nSamples = numSamples,
        nGenes = nrow(markerCounts),
        simSeed = seed,
        simMarkers = TRUE)
    sim_Counts_Markers[!paramMarkers$detectG,] = as.matrix(markerCounts[!paramMarkers$detectG,])

    sim_Counts_Others = simulateCounts(
        simParam = paramOthers,
        nSamples = numSamples,
        nGenes = nrow(otherCounts),
        simSeed = seed,
        simMarkers = TRUE
    )
    sim_Counts_Others[!paramOthers$detectG,] = as.matrix(otherCounts[!paramOthers$detectG,])

    # Rename genes of simulated data, join and create sample IDs
    if(verbose){message("Preparing output ...")}
    rownames(sim_Counts_Markers) = rownames(markerCounts)
    rownames(sim_Counts_Others) = rownames(otherCounts)

    res = rbind(sim_Counts_Markers, sim_Counts_Others)
    colnames(res) = paste0("simu_", 1:numSamples)

    if(verbose){message("Done!")}
    return(res)
}

estimateParameters = function(
    countData,
    simMarkers = TRUE
){
    # Set parameters
    sigma = 1.96

    # Kick out empty samples and keep only expressed genes
    totalS = ncol(countData)
    totalG = nrow(countData)
    fullS = colSums(countData, na.rm = TRUE) > 0
    detectG = rowMeans(countData, na.rm = TRUE) > 0
    countData = countData[detectG, fullS]

    nsamples = dim(countData)[2]
    counts0 = countData == 0
    nn0 = rowSums(!counts0)

    # the negative binomial
    mu = rowSums(countData) / ncol(countData)
    s2 = rowSums((countData - mu) ^ 2) / ncol(countData)
    size = mu ^ 2 / (s2 - mu + 1e-04)
    size = ifelse(size > 0, size, NA)
    p0 = (nsamples - nn0) / nsamples
    mu = mu[!is.na(size)]
    p0 = p0[!is.na(size)]
    remove = rownames(countData)[is.na(size)]
    detectG[names(detectG) %in% remove] = FALSE
    size = size[!is.na(size)]
    phi.g = 1 / size
    phi.c = mean(phi.g)
    ldisp = log2(phi.g)
    lsize = log2(size)
    lmu = log2(mu + 1)

    estG = length(mu)
    estS = length(ncol(countData))

    # meansizefit
    meansizefit = loess.sd(lsize ~ lmu, nsigma = sigma)

    # meandispfit
    meandispfit = loess.sd(ldisp ~ lmu, nsigma = sigma)

    # return object
    paramData = list(means = mu,
                     dispersion = phi.g,
                     common.dispersion = phi.c,
                     size = size,
                     p0 = p0,
                     meansizefit = meansizefit,
                     meandispfit = meandispfit,
                     estS = estS,
                     estG = estG,
                     totalS = totalS,
                     totalG = totalG,
                     detectG = detectG,
                     sigma = sigma)

    return(paramData)
}

simulateCounts = function(
    simParam,
    nSamples,
    nGenes,
    simSeed = NULL,
    simMarkers = TRUE
){
    if(simMarkers){
        lfcs = as.matrix(rep(0, nGenes))
    } else {
        lfcs = as.matrix(rep(0, nGenes))
    }

    # define NB params
    mu = simParam$means
    meansizefit = simParam$meansizefit

    # For markers use observed mean parameters, for other genes sample
    if(simMarkers) {
        present = simParam$detectG
        if(sum(present) < nGenes){
            warning("Detected one or more marker genes with no expression
                    in reference samples!")
            true.means = rep(0, nGenes)
            true.means[present] = mu
        } else {
            true.means = mu
        }
    } else {
        index = sample(1:length(mu), size = nGenes, replace = TRUE)
        true.means = mu[index]
    }

    # estimate size parameter associated with true mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule = 2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule = 2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # size factor
    all.facs = rep(1, nSamples)
    #all.facs = sample(seq(1, 2.5, by = 0.1), size = nSamples, replace = TRUE)

    # effective means
    effective.means = outer(true.means, all.facs, "*")

    mod = as.matrix(rep(1, nSamples))

    # make mean expression with beta coefficients added as defined
    # by model matrix
    mumat = log2(effective.means + 1) + lfcs %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(
        rnbinom(nSamples * nGenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
        ncol = nSamples,
        nrow = nGenes,
        dimnames = list(paste0(rownames(mumat),"_", seq_len(nGenes)),
                        NULL))

    return(counts)
}

simulateNegativeControls = function(
    nGenes,
    numSamples,
    normMean = 50,
    normSD = 500,
    seed = NULL,
    verbose = TRUE
){
    if(is.null(seed)) {
        seed = sample(1:1000000, size = 1)
    }
    #set.seed(seed)

    # Simulate counts based on randomly sampled mean and size values
    # for negative binomial distribution
    means = rnorm(nGenes, mean = normMean, sd = normSD)
    means[means < 0] = 0
    all.facs = sample(seq(1, 3, by = 0.1), size = numSamples, replace = TRUE)
    effective.means = outer(means, all.facs, "*")
    mumat = log2(effective.means + 1)
    mumat[mumat < 0] = min(log2(effective.means + 1))
    sizevec = rnorm(nGenes, mean = 0, sd = 4)

    counts = matrix(
        rnbinom(numSamples * nGenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
        ncol = numSamples,
        nrow = nGenes
    )
    colnames(counts) = paste0("Negative_Control_", c(1:numSamples))

    return(counts)
}
