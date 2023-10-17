library(ggplot2)
library(viridis)
library(ggsci)
library(ggrepel)
cols <- c(
  "#BC3C29FF", "#E18727FF", "#0072B5FF", "#20854EFF", "#7876B1FF",
  "#6F99ADFF", "#FFDC91FF"
)
freqBarChart <- function(plotdb) {
  p <- ggplot(plotdb, aes(fill = mutType, y = freq, x = study)) +
    # geom_bar(position = "stack", stat = "identity", width = 0.85) +
    geom_bar(position = "stack", stat = "identity", width = 0.8) +
    theme_classic() +
    # coord_flip() +
    xlab("") +
    scale_y_continuous(
      limits = c(0, 10),
      expand = c(0, 0),
      breaks = seq(0, 10, 2),
      name = "Mutation Frequency (%)"
    ) +
    theme(
      legend.position = c(0.8, 0.5),
      # # legend.direction = "vertical",
      text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1)
    )
  return(p)
}
mutdb <- readRDS("./results/mutdb.rds")
mutdb$type <- mutdb$Variant_Classification

mutdb$type <- gsub("Frame_Shift_Del", "Frameshift Deletion", mutdb$type)
mutdb$type <- gsub("In_Frame_Ins", "InFrame Insertion", mutdb$type)
mutdb$type <- gsub("_", " ", mutdb$type)



clindb <- readRDS("./results/clindb.rds")
clindb$cumFreq <- 0
mutTypes <- unique(mutdb$type)

mutTypeFreqs <- list() # <- setNames(rep(0, length(mutTypes)), mutTypes)

plotdb <- data.frame()
for (i in seq_len(nrow(clindb))) {
  curStudy <- clindb$study[i]
  for (curMutType in mutTypes) {
    submutdb <- mutdb[mutdb$study == curStudy & mutdb$type == curMutType, ]
    freq <- nrow(submutdb) / clindb$size[i] * 100
    if (freq == 0) next()
    t <- data.frame(study = curStudy, mutType = curMutType, freq = freq)
    plotdb <- rbind(t, plotdb)
    clindb$cumFreq[i] <- clindb$cumFreq[i] + freq
    mutTypeFreqs[[curMutType]] <- append(freq, mutTypeFreqs[[curMutType]])
  }
}

meanMutTypeFreqs <- sapply(mutTypeFreqs, mean)
meanMutTypeFreqs <- sort(meanMutTypeFreqs, decreasing = TRUE)
mutTypeOrder <- names(meanMutTypeFreqs)

studyOrder <- clindb$study[order(clindb$cumFreq, decreasing = TRUE)]
plotdb$study <- factor(plotdb$study, levels = studyOrder)
plotdb$mutType <- factor(plotdb$mutType, levels = mutTypeOrder)

plotdb <- plotdb[plotdb$freq > 0, ]


pdf("./figures/bar.pdf", height = 4, width = 8)
plot(freqBarChart(plotdb) + scale_fill_manual(values = cols, name = "Mutation types"))
# scale_fill_nejm("default", name = "Mutation types"))

dev.off()
print("done")
