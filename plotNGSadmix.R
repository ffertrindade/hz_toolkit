# Plot NGSadmix results from *.qopt files

args <- commandArgs(T)
if (length(args)!=3) {
  stop("Usage: plotNGSadmix.R qopt_file label_file k_num", call.=FALSE)
}

qopt <- args[1]
label <- args[2]
k <- args[3]
rm(args)

# Get ID and pop info for each individual
pop<-scan(label)

# Read inferred admixture proportions file
q<-read.table(qopt)

# Plot them (ordered by population or individual name)
ord = order(pop)

filename <- paste(qopt, ".png", sep = "")
y_title <- paste("Admixture proportions (k=", k, ")", sep = "")
png(filename = filename, width = 900, height = 300)
par(mar = c(5,4,1,1))

barplot(t(q)[,ord],
        col = c(4,5,8,1,3),
        names = pop[ord],
        las = 2,
        ylab = y_title,
        main = qopt,
        cex.axis = 1, 
        cex.names = 1)

dev.off()

# order according to population only
filename <- paste(qopt, "_b.png", sep = "")
png(filename = filename, width = 500, height = 300)
par(mar = c(5,4,1,1))

barplot(t(q)[,ord],
        col = c(4,5,8,1,3),
        space = 0,
        border = NA,
        xlab = "Individuals from each population",
        ylab = y_title,
        main = qopt,
        cex.axis = 1)
text(tapply(1:length(pop), pop[ord], mean), -0.05, unique(pop[ord]), xpd=T)
abline(v = cumsum(sapply(unique(pop[ord]), function(x){sum(pop[ord]==x)})), col=1, lwd=1.2)
dev.off()
