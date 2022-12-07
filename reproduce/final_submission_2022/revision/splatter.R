suppressPackageStartupMessages({
  library(splatter)
  library(scater)
})

set.seed(1)
sce <- mockSCE()
params <- splatEstimate(sce)
sim <- splatSimulate(params)

# poc
params.groups <- newSplatParams()
sim <- splatSimulateGroups(params.groups,
                           batchCells=c(3000),
                           group.prob=c(0.15,0.15,0.23,0.23,0.24),
                           de.prob=c(0.4,0.2,0.15,0.15,0.15))
sim <- logNormCounts(sim)
sim <- runPCA(sim)
plotPCA(sim, colour_by = "Group")

library(zellkonverter)
writeH5AD(sim,file='/Users/ligk2e/Desktop/sim.h5ad')

# discovery and conversed mode
params.groups <- newSplatParams()
sim <- splatSimulateGroups(params.groups,
                           batchCells=c(3000),
                           group.prob=c(0.2,0.2,0.15,0.15,0.15,0.15),
                           de.prob=c(0.2,0.2,0.005,0.005,0.005,0.005),
                           de.facLoc=c(0.1,0.1,0.5,0.5,0.5,0.5))
sim <- logNormCounts(sim)
sim <- runPCA(sim)
plotPCA(sim, colour_by = "Group")

library(zellkonverter)
writeH5AD(sim,file='/Users/ligk2e/Desktop/sim_dc.h5ad')
