# 4/10/2019
# Hua Sun
# Rscript estimates.R -input exp.input -name example -outdir /path/result

# exp.input
# GeneSymbol	C3L-00004	C3L-00010	C3L-00011
# A2ML1	0.040932437	0.053669561	0.048588662
# A2MP1	0.857802521	0.686595892	0.143925057
# A3GALT2	0.026564435	0.020165087	0.211973116


#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
#library(estimate)


library(GetoptLong)

GetoptLong(
	"input=s", "matrix",
	"prefix=s", "prefix name",
	"outdir=s", "out dir folder"
)


library(estimate)

# output_file: filtered genes
filtered_f = paste0(outdir, "/", prefix, "_genes.gct")
# output_file: tumor purity score
score_f = paste0(outdir, "/", prefix, "_estimate_score.gct")
# output_final
output = paste0(outdir, "/", prefix, "_estimate_score.tsv")

# only get <=10412 genes
filterCommonGenes(input.f=input, output.f=filtered_f, id="GeneSymbol")

estimateScore(input.ds = filtered_f, output.ds=score_f, platform="affymetrix")
scores=read.table(score_f, skip = 2, header = T)

rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])  # important. Don't change
# scores - (sample name as row)
#		StromalScore	ImmuneScore	ESTIMATEScore	TumorPurity

# if the sample as X8888 then format as 8888 
rownames(scores) <- sub('X', '', rownames(scores))


# output the tumor purity information
write.table(scores, output, sep = "\t", quote = FALSE)






# draw the purity plot
library(ggplot2)


if (nrow(scores)>30){

	p <- ggplot(as.data.frame(scores), aes(x=rownames(scores), y=TumorPurity)) +
				geom_bar(fill = "#45B8AC", stat = "identity") +
				#   geom_text(aes(label = Freq), vjust = 0, hjust = .5, size = 3) + 
				#   theme_pubclean() + 
				theme_classic() + 
				geom_hline(yintercept=.3, linetype="dashed", color = "red") +
				geom_hline(yintercept=.5, linetype="dashed", color = "blue") +
				aes(x=reorder(rownames(scores), -TumorPurity)) +
				labs(x="Samples", y="Tumor purity") + 
				ylim(0, 1) + 
				theme(axis.text.x=element_blank(),
        			axis.ticks.x=element_blank())
				
} else {

	p <- ggplot(as.data.frame(scores), aes(x=rownames(scores), y=TumorPurity)) +
				geom_bar(fill = "#45B8AC", stat = "identity") +
				#   geom_text(aes(label = Freq), vjust = 0, hjust = .5, size = 3) + 
				#   theme_pubclean() + 
				theme_classic() + 
				geom_hline(yintercept=.3, linetype="dashed", color = "red") +
				geom_hline(yintercept=.5, linetype="dashed", color = "blue") +
				aes(x=reorder(rownames(scores), -TumorPurity)) +
				labs(x="Samples", y="Tumor purity") + 
				ylim(0, 1) + 
				theme(axis.text.x=element_text(angle=90, hjust=1), axis.text=element_text(size=6))
				
}


pdf(paste0(output, ".pdf"), width = 5, height = 2)
print(p)
dev.off()

