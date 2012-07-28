#!/usr/bin/env Rscript

# 2012-04-18
# works with R 15.0 and ggplot 0.9. 
# Check ggplot2 help forums or contact sujai.kumar@ed.ac.uk if something doesn't run
# because of updated programs/packages

blob.plot<-function(dfilt, nc, nr) {
    theme_set(theme_bw())

    colorbrewer=c("#DDDDDD", "#777777", "#4575B4", "#ABD9E9", "#313695", "#E0F3F8", "#FEE090", "#74ADD1", "#FFFFBF", "#F46D43", "#D73027", "#FDAE61", "#A50026")
    paultol    =c("#DDDDDD", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#44AA99", "#999933", "#AA4499", "#882255", "#777777")

    g <- ggplot() + scale_colour_manual(values=paultol, name="Taxonomic\nClassification", limits=levels(dfilt$taxon) )
    for (t in levels(dfilt$taxon)) {
      g <- g + geom_point(data=dfilt[dfilt$taxon==t,],aes(gc, cov, colour=taxon), size=2, alpha=I(1/3))
    }
    axis_breaks = c(1,2,5,10,20,50,100,200,500,1000);
    g +
      facet_wrap(~read_set, nrow=nr, ncol=nc) + 
      xlim(0,1) + ylim(0,max(dfilt$cov)) + scale_y_log10(breaks = axis_breaks, labels = axis_breaks) +
      labs(x="GC content", y="Read coverage") + 
      guides(colour = guide_legend(override.aes = list(alpha = 1,size=10))) + 
      opts (
        strip.text.x = theme_text(colour = "black", size = 25, vjust = 0.5),
        axis.text.x  = theme_text(colour = "black", size = 25, vjust = 1),
        axis.text.y  = theme_text(colour = "black", size = 25, vjust = 0.5),
        axis.title.x = theme_text(colour = "black", size = 25, vjust = 0),
        axis.title.y = theme_text(colour = "black", size = 25, hjust = 0.5, vjust = 0.5, angle=90),
        legend.text  = theme_text(colour = "black", size = 25, vjust = 0),
        legend.title = theme_text(colour = "black", size = 25, vjust = 0, hjust = 0, lineheight=1),
        legend.key.height = unit(2.2,"line"),
        legend.justification=c(1,1), legend.position=c(1,1)
      )
}

#Function to ignore low frequency annotations:

clean.blobs<-function(d,threshold) {
    annotated<-d[d$taxon!="Not annotated",]
    total<-dim(annotated)[1]
    levels(d$taxon)[which(table(d$taxon)<threshold*total)]<-"Not annotated"
    return(d)
}

#Load data from file and generate plot:

library(ggplot2)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
arg_input_file <- args[1]
arg_ignore_below_prop=as.numeric(args[2])
arg_num_cols=as.numeric(args[3])
arg_num_rows=as.numeric(args[4])
orig <- read.delim(arg_input_file,header=FALSE,sep="\t", col.names=c("read_set","contigid","len","cov","gc","taxon"))
dfilt <- clean.blobs(orig,arg_ignore_below_prop)
taxa  <- rownames(rev(sort(table(dfilt$taxon))))
taxa  <- taxa[taxa != "Not annotated"]
dfilt$taxon<-ordered(dfilt$taxon,levels=c("Not annotated",taxa)) #put the taxon category with the most hits first so it is plotted first (background)
png(paste(arg_input_file,"png",sep="."), (arg_num_cols * 1000), (arg_num_rows * 1000),units="px",res=100)
blob.plot(dfilt,arg_num_cols,arg_num_rows)
dev.off()
