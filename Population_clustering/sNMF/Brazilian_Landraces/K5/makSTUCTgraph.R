## Script for creating barplot of K5
## The genotypes are sorted based on PC1 for mesoamerican and Andean,
## While within Mesoamerican are sorted by PC3
## The genotypes belonging to K1 are called Meso 4 cause they cluster with Meso in K2, even though in our analysis they are admized

a=read.table('allInfo.csv', sep=',', header=T)
mat=a[,2:6]
colP=c('#a50f15', '#fcae91', 'blue', '#fb6a4a', '#de2d26')

pdf('K5STUCTplot.pdf', width=20, height=6)
par(mar=c(0.5,4,1,0.5))
barplot(t(mat), col=colP, ylim=c(-0.2, 1.2), xlim=c(0,355), yaxt='n', ylab='Q') #main='Brazilian collection K2')
axis(2, at=seq(0, 1, 0.2))

colL=c('blue' ,'#fcae91','#fb6a4a', '#de2d26','#a50f15')
legend('right', c('Andean', 'Mesoamerican 1', 'Mesoamerican 2', 'Mesoamerican 3', 'Mesoamerican 4'), text.col='white', col=colL, pch=15, bty='n', cex=1.4)
legend('right', c('Andean', 'Mesoamerican 1', 'Mesoamerican 2', 'Mesoamerican 3', 'Mesoamerican 4'), col='black', pch=0, bty='n', cex=1.4)



