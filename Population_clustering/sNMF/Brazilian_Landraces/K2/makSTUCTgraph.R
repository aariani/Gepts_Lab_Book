## Script for creating barplot of K2
## The genotypes are sorted based on PC1 in the file

a=read.table('K2_matrix.csv', sep=',', header=T)
## The K2 matrix is the Ancestri coefficienty (col 2 and 3)
## And the PC1 coordinates (for ordering the data)
mat=a[,2:3]
colP=c('blue', 'red')
pdf('K2STUCTplot.pdf', width=20, height=6)
par(mar=c(0.5,4,1,0.5))
barplot(t(mat), col=colP, ylim=c(-0.2, 1.2), xlim=c(0,350), yaxt='n', ylab='Q') #main='Brazilian collection K2')
axis(2, at=seq(0, 1, 0.2))
legend('right', c('Andean', 'Mesoamerican'), text.col='white', col=colP, pch=15, bty='n', cex=1.4)
legend('right', c('Andean', 'Mesoamerican'), col='black', pch=0, bty='n', cex=1.4)



