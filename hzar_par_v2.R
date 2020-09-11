#!/usr/bin/env Rscript

## Check for the arguments
mincenter = -115
maxcenter = 3000
maxwidth = 3115
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("No input file!!!.n", call.=FALSE)
}

## Load the packages
library(hzar)

if(require(doMC)){
  registerDoMC(5)
} else {
  registerDoSEQ()
}
#setwd("G:/OneDrive - PUCRS - BR/Doutorado/analises/hzar")

## Load molecular data from the data table.
#leopardus <- read.table("mknExLoci.txt",header=TRUE)#, sep=",")
leopardus <- read.table(args[1], header=TRUE, sep=",")

## Blank out space in memory to hold molecular analysis
if(length(apropos("^leopd$",ignore.case=FALSE)) == 0 || !is.list(leopd) ) {
  leopd <- list()
}

## Space to hold data and results
obs <- list() # Space to hold the observed data
models <- list() # Space to hold the models to fit
fitRs <- list() # Space to hold the compiled fit requests
runs <- list() # Space to hold the output data chains
analysis <- list() # Space to hold the analysed data

## Creating data object for SNP
print("##### Running hzar.doMolecularData1DPops #####")
obs <- hzar.doMolecularData1DPops(leopardus$dist,
                                  leopardus$freq,
                                  leopardus$samples,
                                  siteID=leopardus$locationID)

## Constructing different clineMetaModel objects by using a function
print("##### Running hzar.makeCline1DFreq #####")
leopd.loadlocmodel <- function(scaling,tails,
                               id=paste(scaling,tails,sep="."))
  models[[id]] <<- hzar.makeCline1DFreq(obs, scaling, tails)

#leopd.loadlocmodel("none", "none", "modelNN")
#leopd.loadlocmodel("none", "right", "modelNR")
#leopd.loadlocmodel("none", "left", "modelNL")
#leopd.loadlocmodel("none", "mirror", "modelNM")
#leopd.loadlocmodel("none", "both", "modelNB")
leopd.loadlocmodel("fixed", "none", "modelXN")
leopd.loadlocmodel("fixed", "right", "modelXR")
leopd.loadlocmodel("fixed", "left", "modelXL")
leopd.loadlocmodel("fixed", "mirror", "modelXM")
leopd.loadlocmodel("fixed", "both", "modelXB")
#leopd.loadlocmodel("free", "none", "modelFN")
#leopd.loadlocmodel("free", "right", "modelFR")
#leopd.loadlocmodel("free", "left", "modelFL")
#leopd.loadlocmodel("free", "mirror", "modelFM")
#leopd.loadlocmodel("free", "both", "modelFB")
model_name <- names(models)

## Modify all models to focus on the region where the observed data were collected.
#models <- sapply(models,
#                 hzar.model.addBoxReq,
#                 1000 , 2000, # Set yhe min and max dist
#                 simplify=FALSE)

## Compile each of the models to prepare for fitting, create an hzar.fitRequest object [demora, paralelo]
print("##### Running in loop hzar.first.fitRequest.old.ML #####")
fitRs$init <- foreach(i=1:length(models)) %dopar% {
	hzar.first.fitRequest.old.ML(models[[i]], obs, verbose = FALSE)
}
names(fitRs$init) <- model_name

## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528),
       D=c(978,544,99,596,528,124),
       E=c(544,99,596,528,124,978),
       F=c(99,596,528,124,978,544))

## Update the settings for the fitter
print("##### Setting chain parameters in loop #####")
seed <- 0
foreach(i=1:length(fitRs$init)) %do% {
  seed = seed +1

  fitRs$init[[i]]$mcmcParam$seed[[1]] <- mainSeed[[seed]]
  fitRs$init[[i]]$mcmcParam$burnin <- fitRs$init[[i]]$mcmcParam$chainLength %/% 10

  if (seed == 6) {
    seed <- 0
  }
}

## Run each model for an initial chain [demora, paralelo]
print("##### Running in loop hzar.doFit #####")
runs$init <- foreach(i=1:length(fitRs$init)) %dopar% {
  hzar.doFit(fitRs$init[[i]])
}
names(runs$init) <- model_name

## Compile a new set of fit requests using the initial chains 
print("##### Running in loop hzar.next.fitRequest #####")
fitRs$chains <- foreach(i=1:length(runs$init)) %dopar% {
  hzar.next.fitRequest(runs$init[[i]])
}
names(fitRs$chains) <- model_name

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel
print("##### Running hzar.multiFitRequest #####")
fitRs$chains <- hzar.multiFitRequest(fitRs$chains,
                                     each = 3,
                                     baseSeed = NULL)

## Just to be thorough, randomize the initial value for each fit
#print("##### Setting cline parameters in loop #####")
#foreach(i=1:length(fitRs$chains)) %do% {

  ## center
#  newcenter <- runif(1, mincenter, maxcenter)
#  fitRs$chains[[i]]$modelParam$init["center"] <- newcenter
  
  ## width
#  newwidth <- runif(1, 0, maxwidth)
#  fitRs$chains[[i]]$modelParam$init["width"] <- newwidth
  
#}

#foreach(i=31:length(fitRs$chains)) %do% { # So rodar para os scaling free
  
  ## pMIn
#  newpmin <- runif(1, 0, 1)
#  fitRs$chains[[i]]$modelParam$init["pMin"] <- newpmin
  
  ## pMax
#  newpmax <- runif(1, 0, 1)
#  fitRs$chains[[i]]$modelParam$init["pMax"] <- newpmax
  
#}

#tails <- 0 # 0-4 modelos de tail, aqui so vai ser mudado o 1,2,4
#replic <- 0 # 0-2 replicas
#foreach(k=1:length(fitRs$chains)) %do% { # rodar separadamente para cada tail
  
#  if (tails > 4) {
#    tails = 0
#  }
  
  ## deltaL and tauL
#  if (tails == 2 || tails == 4) {
#    newdelL <- runif(1, 0, maxwidth)
#    fitRs$chains[[k]]$modelParam$init["deltaL"] <- newdelL
#    newtauL <- runif(1, 0, 1)
#    fitRs$chains[[k]]$modelParam$init["tauL"] <- newtauL
#  }
  
  ## deltaR and tauR
#  if (tails == 1 || tails == 4) {
#    newdelR <- runif(1, 0, maxwidth)
#    fitRs$chains[[k]]$modelParam$init["deltaR"] <- newdelR
#    newtauR <- runif(1, 0, 1)
#    fitRs$chains[[k]]$modelParam$init["tauR"] <- newtauR
#  }
  
#  if (replic < 2) {
#    replic = replic + 1
#  } else {
#    replic = 0
#    tails = tails + 1
#  }
  
#}

## Run a chain of 3 runs for every fit request [demora, paralelo]
print("##### Running hzar.doChain.multi #####")
runs$chains <-  hzar.doChain.multi(fitRs$chains,
                                   doPar = TRUE,
                                   inOrder = TRUE,
                                   count = 3)

## Summary result parameters and log likelihood [ ver como printar em um output ]
print("##### Summary result parameters and log likelihood #####")
foreach(i=1:length(runs$chains)) %dopar% {
  if (i == 1 || i == 4 || i == 7 || i == 10 || i == 13 || i == 16 || i == 19 || i == 22 || i == 25 ||
      i == 28 || i == 31 || i == 34 || i == 37 || i == 40 || i == 43) {
    summary(do.call(mcmc.list,
                    lapply(runs$chains[i:i+2],
                           function(x) hzar.mcmc.bindLL(x[[3]])))) # Pega a terceira chain de cada uma das 3 replicas   
  }
}

## Create a model data group (hzar.dataGroup object) for each model from the initial runs [demora pouco]
print("##### Running in loop hzar.dataGroup.add #####")
analysis$initDGs <- foreach(i=1:length(runs$init)) %dopar% {
  hzar.dataGroup.add(runs$init[[i]])
}
names(analysis$initDGs) <- model_name

## Create a model data group for the null model (expected allele frequency independent of distance along cline)
print("##### Running hzar.dataGroup.null #####")
analysis$initDGs$nullModel <- hzar.dataGroup.null(obs)

## Create a hzar.obsDataGroup object from the 16 hzar.dataGroup just created
print("##### Running hzar.make.obsDataGroup #####")
analysis$oDG <- hzar.make.obsDataGroup(analysis$initDGs)
analysis$oDG <- hzar.copyModelLabels(analysis$initDGs, analysis$oDG)

## Convert all 135 runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object [demora pouco]
print("##### Running hzar.make.obsDataGroup hzar.dataGroup.add #####")
analysis$oDG <- hzar.make.obsDataGroup(lapply(runs$chains,
                                              hzar.dataGroup.add), analysis$oDG)

## Compare the 16 cline models to the null model graphically
print("##### Running hzar.plot.cline #####")
output <- paste(args[1], ".all_clines.png", sep = "")
png(width=400, height=400, filename=output, pointsize=8)
hzar.plot.cline(analysis$oDG)
dev.off()

## Do model selection based on the AICc scores
print("##### Running hzar.AICc.hzar.obsDataGroup #####")
print(analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(analysis$oDG))

## Print out the model with the minimum AICc score
print("##### Model with the minimum AICc #####")
output <- paste(args[1], ".best_model.txt", sep = "")
print(analysis$model.name <- rownames(analysis$AICcTable
                                )[[ which.min(analysis$AICcTable$AICc )]])
write(analysis$model.name, file=output)

## Extract the hzar.dataGroup object for the selected model
analysis$model.selected <- analysis$oDG$data.groups[[analysis$model.name]]

## Look at the variation in parameters for the selected model
print("##### Running hzar.getLLCutParam #####")
if (analysis$model.name != "nullModel") {
  print(hzar.getLLCutParam(analysis$model.selected,
                           names(analysis$model.selected$data.param)))
  write(analysis$model.name, file=output, append = TRUE)
} else {
  print("nullModel was the best model")
}

## Print the maximum likelihood cline for the selected model
print("##### Running hzar.get.ML.cline #####")
hzar.get.ML.cline(analysis$model.selected)

## Plot the maximum likelihood cline with 95% credible cline region for the selected model
print("##### Running hzar.plot.cline #####")
output <- paste(args[1], ".best_cline.png", sep = "")
png(width=400, height=400, filename=output, pointsize=8)
hzar.plot.cline(analysis$model.selected)
dev.off()
