# plotting script, meant to be called by compare_dadi.py

library(ggplot2)

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
    table <- read.table(file('stdin'),head=T)
    table$G <- factor(table$G)
    levels(table$G)[1] <- 'momi'
    p <- ggplot(table, aes(x=n,y=seconds,color=G,linetype=G)) + scale_y_log10() + geom_point(aes(shape=G)) + geom_line() + scale_shape_manual(values=seq(0,8))
    ggsave('dadi_timing.pdf',p)
}

accuracy <- function(){
    table <- read.table(file('stdin'),head=T)
    table$abs.ratio <- abs(table$val / table$momi)
    table$sign <- '>=0'
    table$sign[table$val < 0] <- '<0'
    table$G <- table$n.pts

    table <- table[sample(nrow(table)),]

    leg = guide_legend(title="sign (dadi)")
    p <- ggplot(table, aes(x=momi,y=abs.ratio, color=sign, shape=sign)) + geom_point() + geom_abline(intercept=0,slope=0) + scale_x_log10() + scale_y_log10() + ylab('|dadi/momi|') + facet_grid(G~n, labeller=label_both) + guides(color=leg, shape=leg) + scale_shape_manual(values=c(8,1))

    ggsave('dadi_accuracy.png',p,width=10,height=10)

    f <- file('ratio_summary.txt')
    cat(summary(table$val / table$momi), file=f)
    close(f)
}

main()
