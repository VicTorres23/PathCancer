library(maftools)
library(ggplot2)

laml <- read.maf("uploads/tcga_laml.maf")

output_dir <- "C:/Users/vemma/Documents/Master in Bioinformatics/Fall_2024_Semester/Post-Genomic Analysis/Project/pathcancer/static/images/"

pw <- pathways(laml)
pw_df <- as.data.frame(pw)

ggplot(pw_df, aes(x = reorder(Pathway, -Mutated_samples), y = Mutated_samples)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() + coord_flip() + labs(title = "Number of Mutated Samples per Pathway", x = "Pathway", y = "Number of Mutated Samples") +
  ggsave("static/images/pathway_barplot.png", width = 10, height = 6)
