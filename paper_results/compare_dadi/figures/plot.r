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
    p <- ggplot(table, aes(x=n,y=seconds,color=G,linetype=G)) + scale_y_log10() + geom_point() + geom_line()
    ggsave('timing.pdf',p)
}

accuracy <- function(){
    table <- read.table(file('stdin'),head=T)
    table$abs.ratio <- abs(table$val / table$momi)
    table$sign <- '>=0'
    table$sign[table$val < 0] <- '<0'
    table$G <- table$n.pts

    table <- table[sample(nrow(table)),]

    p <- ggplot(table, aes(x=momi,y=abs.ratio, color=sign)) + geom_point() + geom_abline(intercept=0,slope=0) + scale_x_log10() + scale_y_log10() + ylab('|dadi/momi|') + facet_grid(G~n, labeller=label_both) + guides(color=guide_legend(title="sign (dadi)"))

    #p <- ggplot(table, aes(x=momi,y=val,color=factor(n))) + geom_point() + geom_abline(intercept=0,slope=1)  + scale_x_log10() + scale_y_log10() + facet_wrap(~ n.pts)

    ggsave('accuracy.png',p,width=10,height=10)

    f <- file('ratio_summary.txt')
    cat(summary(table$val / table$momi), file=f)
    close(f)
}

main()
