# plotar resultados NGSadmix (adaptado direto do site)

setwd("")

# Get ID and pop info for each individual
pop<-scan("poplabel_v4",what="theFuck")
ind<-scan("indlabel_v4",what="theFuck")

# Read inferred admixture proportions file
q<-read.table("batch001_run003_filtered_k5.qopt")

# Plot them (ordered by population or individual name)
ord = order(pop)
ordi = order(ind)

png(filename="plots/run003_ld_ngsadmix_k5.png", width=900, height=300)
par(mar=c(5,4,1,1))
barplot(t(q)[,ordi],
        col=c(4,5,8,1,3),
        names=ind[ordi],
        las=2,
        ylab="Admixture proportions (k=5)",
        main="run003 - pruned due to ld",
        cex.axis=1, 
        cex.names=1)
dev.off()

# order according to population only
png(filename="plots/batch001_run002_ngsadmix_k3_pop.png", width=500, height=300)
par(mar=c(5,4,1,1))
barplot(t(q)[,ord],
        col=c(4,5,8,1,3),
        space=0,
        border=NA,
        xlab="Individuals from each population",
        ylab="Admixture proportions (K=3)",
        main="Batch001 - run002",
        cex.axis=1)
text(tapply(1:length(pop),pop[ord],mean),-0.05,unique(pop[ord]),xpd=T)
abline(v=cumsum(sapply(unique(pop[ord]),function(x){sum(pop[ord]==x)})),col=1,lwd=1.2)
dev.off()