## script for 3d plot of PCA (1-3)

library(rgl)
a=read.table('allInfo.csv', sep=',', header=T,  as.is=T, comment.char='')

### The allInfo.csv table contains the 3 PC coordinates (based on adegenet), and the colors based on K5 clustering and Q >=0.7

plot3d(a$PC1, a$PC3, a$PC2, type='s', size=2, xlab='PC1', ylab='PC3', zlab='PC2', box=F, col=a$colP)

## For saving this kind of plot you have to do it manually

snapshot3d(filename = '3dplot.png', fmt = 'png')
