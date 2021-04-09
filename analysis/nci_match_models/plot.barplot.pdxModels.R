
f_table <- './drugTargetGenes.mergedTable.pdxmodel.out'

d_merge_all <- fread(f_table, data.table = F)

library(ggplot2)
func_barPlot <- function(df){
    
    #set_x_axis = c("BLCA", "BRCA", "CESC", "COAD", "CSCC", "GBM", "GIAD", "GIST", "HNSC", "KIRC", "LUAD", "LUAS", "LUSC", "MCC", "MESO", "OV", "PAAD", "PRAD", "READ", "SARC", "SCLC", "SKCM", "STAD", "UCEC", "UCS", "Other")
    set_x_axis = c("BLCA", "BRCA", "CESC", "COAD", "CSCC", "GBM", "GIAD", "HNSC", "KIRC", "LUAD", "LUSC", "MCC", "MESO", "OV", "PAAD", "PRAD", "READ", "SARC", "SCLC", "SKCM", "STAD", "UCEC", "UCS", "Other")
    
    N_sum <- sum(df$n)
    
    p <- ggplot(df, aes(x=CancerType, y=n)) +
            geom_bar(stat="identity", width=0.6, position=position_dodge(), fill="#A8427D") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
            geom_text(aes(label = n), position = position_stack(vjust = 1.2), size=2.5) +
            theme(axis.title.x = element_blank(),
                  axis.title = element_text(size=6),
                  #axis.text.y = element_text(size=6),
                  axis.text.y = element_blank()
                 ) +
            theme(
                  axis.line.y = element_blank(),
                  axis.line.x = element_line(colour = 'black', size = 0.2),
                  axis.ticks = element_blank()
                 ) +
            ylab("# of unique \nPDX models") +
            scale_x_discrete(limit = set_x_axis) +
            annotate("text", x=23, y=60, label=paste0("Total N=",N_sum), size=2.5)

    
    return(p)
}



pdf(paste0("drugGene.pdxmodel.barplot.pdxmodel.", Sys.Date(), ".pdf"), width = 4.4, height = 1)
func_barPlot(d_merge_all)
dev.off()

