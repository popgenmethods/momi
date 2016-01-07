### R plotting script called by benchmark.py

library(ggplot2)
library(data.table)
library(reshape2)
library(plyr)

main <- function(){
    args <- commandArgs(TRUE)
    if(args[1] == 'timing'){
        timing()
    } else if (args[1] == 'accuracy'){
        accuracy()
    } else {
        print("Unrecognized argument")
    }
}

timing <- function(){
    dt <- data.table(read.table(file('stdin'),head=T,stringsAsFactors=F))
    dt$D <- as.factor(dt$D)

    dt <- dt[, list(time=mean(time)), by=list(method, n, D, component)]

    p <- ggplot(dt, aes(x=n, y=time, color=D)) + geom_line() + geom_point() + ylab("Seconds") + theme(legend.position="top") + facet_grid(component ~ method, scales='free_y')
    ggsave('figures/timing.pdf',p)
#    ggsave('figures/timing_bw.pdf',p+scale_color_grey())

    p <- ggplot(dt, aes(x=n, y=time, color=D, linetype=method)) + geom_line() + geom_point(aes(shape=method)) + ylab("Seconds") + theme(legend.position="top") + facet_grid(component ~ .) + scale_y_log10(breaks=c(1e-4,1e-2,1,1e2,1e4)) + scale_x_log10(breaks=c(2,4,8,16,32,64,128,256))
    ggsave('figures/timing_log.pdf',p)
#    ggsave('figures/timing_log_bw.pdf',p+scale_color_grey())    

    write.table(format(dcast(subset(dt, component=='Per SNP'), formula=n+D~method, value.var='time'),digits=4,scientific=T), file='figures/per_snp.txt', quote=F, row.names=F, sep=' & ', eol=' \\\\ \\hline \n')

    write.table(format(dcast(subset(dt, component=='Precomputation'), formula=n+D~method, value.var='time'),digits=4,scientific=T), file='figures/precomputation.txt', quote=F, row.names=F, sep=' & ', eol=' \\\\ \\hline \n')
}

accuracy <- function(){
    dt <- read.table(file('stdin'), head=T)

    dt <- dt[sample(nrow(dt)),]
    dt$n <- factor(dt$n)

    library(scales)
    sgn_log <- function(x) sign(x) * log10(1 + abs(x))
    inv_sgn_log <- function(x) sign(x) * (10^(abs(x)) - 1)
    sgn_log_trans <- function() trans_new("sgn_log", sgn_log, inv_sgn_log, trans_breaks(sgn_log, inv_sgn_log))
    ##p <- ggplot(dt, aes(x=momi,y=Chen,color=n)) + geom_abline(color='red', linetype='dashed') + geom_point() + scale_y_continuous(trans = sgn_log_trans()) + scale_x_continuous(trans = sgn_log_trans(), labels=function(x) as.character(round(x,2))) + scale_colour_brewer(palette="YlOrRd")
    p <- ggplot(dt, aes(x=momi,y=Chen,color=n,shape=n)) + geom_abline(color='red', linetype='dashed') + geom_point(aes(shape=n)) + scale_y_continuous(trans = sgn_log_trans()) + scale_x_continuous(trans = sgn_log_trans(), labels=function(x) as.character(round(x,2))) + scale_shape_manual(values=seq(0,8))
    ggsave('figures/accuracy.png', p)
#    ggsave('figures/accuracy_bw.png', p+scale_color_grey())
}


main()
