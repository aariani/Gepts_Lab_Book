## script for 2d plot of PCA (1-2)

a=read.table('allinfo.csv', sep=',', header=T)

## The allinfo.csv file contains the ID of the genotype, the PC1 and PC2 coordinates (with adegenet) and the clustering with sNMF
## divided with A, M, and H as symbols

colC=NULL
x=1
for (i in unique(a$ID)){
	s=subset(a, a$ID==i)
	if (s$sNMF_K2=='M') {
		colC[x]=rgb(red=1, green=0, blue=0, alpha=0.6)
		}
	else if (s$sNMF_K2=='A'){
		colC[x]=rgb(red=0, green=0, blue=1, alpha=0.6)
		}
	else {colC[x]=rgb(red=1, green=1, blue=1, alpha=0.4)}
	x=x+1
	}

plot(a$PC1, a$PC2, pch=1, cex=2.2, xlab='PC1', ylab='PC2', main='PCA analysis of Brazilian Landraces')
abline(h=0, v=0)
points(a$PC1, a$PC2, pch=19, col=colC, cex=2)

legend('topright', c('Mesoamerican', 'Andean', 'Hybrid'), col='black', pch=1, text.col='white', bty='n')

legCol=c(rgb(red=1, green=0, blue=0, alpha=0.6), rgb(red=0, green=0, blue=1, alpha=0.6), rgb(red=1, green=1, blue=1, alpha=0.4))
legend('topright', c('Mesoamerican', 'Andean', 'Hybrid'), text.col='white', col='black', pch=1, bty='n')
legend('topright', c('Mesoamerican', 'Andean', 'Hybrid'), col=legCol, pch=19, bty='n')

