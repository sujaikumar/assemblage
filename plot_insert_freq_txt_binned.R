#!/usr/bin/env Rscript

# 2012-07-09
# works with R 15.0 and ggplot 0.9. 
# Check ggplot2 help forums or contact sujai.kumar@ed.ac.uk if something doesn't run
# because of updated programs/packages

args <- commandArgs(trailingOnly = TRUE)
arg_input_file <- args[1]

d <- read.table(arg_input_file,header=FALSE,col.names=c("lib","orientation","dist","freq"))
library(ggplot2)
library(grid)
theme_set(theme_bw())
png( paste(arg_input_file,"png",sep="."),1200,300 * length(levels(d$lib)) )
pairlin1 <- qplot(dist,geom="histogram",weight=freq,data=d,facets = lib ~ orientation, xlab='Insert size',ylab='Count',binwidth=10) + xlim(0,800)
pairlin2 <- qplot(dist,geom="histogram",weight=freq,data=d,facets = lib ~ orientation, xlab='Insert size',ylab='Count',binwidth=100) + xlim(0,6000)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
vplayout <- function(x, y)
viewport(layout.pos.row = x, layout.pos.col = y)
print(pairlin1, vp = vplayout(1, 1))
print(pairlin2, vp = vplayout(1, 2))
dev.off()
