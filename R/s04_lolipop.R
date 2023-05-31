library(ggplot2)
library(viridis)
library(ggsci)
library(ggrepel)
cols <- c(
    "#BC3C29FF", "#E18727FF", "#0072B5FF", "#20854EFF", "#7876B1FF",
    "#6F99ADFF", "#FFDC91FF"
)
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

mutdb$IMPACT <- factor(mutdb$IMPACT, levels = c("LOW", "MODERATE", "HIGH"))
mutdb$pos <- as.numeric(gsub("/656", "", mutdb$Protein_position))
mutdb$type <- factor(mutdb$type, levels = mutTypeOrder)

newMutDb <- as.data.frame(table(mutdb$pos, mutdb$type, mutdb$IMPACT, mutdb$Amino_acids))
colnames(newMutDb) <- c("pos", "type", "impact", "AA", "num")
newMutDb <- newMutDb[newMutDb$num > 0, ]
print(head(newMutDb))
newMutDb$pos <- as.numeric(as.character(newMutDb$pos))
newMutDb$num <- as.numeric(newMutDb$num)
# newMutDb$impact <- factor(newMutDb$impact, levels = c("LOW", "MODERATE", "HIGH"))
newMutDb$impact <- factor(newMutDb$impact, levels = c("HIGH", "MODERATE", "LOW"))
newMutDb$type <- factor(newMutDb$type, levels = mutTypeOrder)
length(levels(newMutDb$type))
phdDomRows <- newMutDb$pos >= 28 & newMutDb$pos <= 76
#sdb <- newMutDb[newMutDb$impact == "HIGH", ]
sdb <- newMutDb[phdDomRows | newMutDb$impact == "HIGH", ]
sdb$label <- apply(sdb, 1, 
                   function(x) {
                   v <- paste0(unlist(strsplit(x["AA"], "/"))[1], x["pos"], unlist(strsplit(x["AA"], "/"))[2])
                   v <- gsub(" ", "", v)
                   if(grepl("NA", v)) {
                       v <- gsub("NA", "", v)
                       v <- paste0(v, substring(v,1, 1))
                   }
                   return(v)
}
)
#newMutDb$num <- newMutDb$num + 0.5
# newMutDb$num <- newMutDb$num * ifelse(phdDomRows, 1, -1)

# ymin <- -.5
# ymax <- .5
ymin <- -1
ymax <- 0

set.seed(123)
p <- ggplot(newMutDb, aes(x = pos,y = num, fill = type, color = type, group = type)) + # size = impact, ,  , 
    annotate("rect",
        #xmin = 60, xmax = 110, 
        xmin = 28, xmax = 76, 
        ymin = 0, ymax = 4,
        alpha = 1, fill = "#F9F1E4"
    ) +
    geom_segment(
        #mapping = aes(xend = pos, yend = 0), 
        mapping = aes(xend = pos, yend = -.5), 
        color = "grey50",
        size = .5, alpha = .85
    ) +
          # better colors
          # put all points up
    geom_hline(yintercept = (ymax + ymin) / 2, linewidth = 3, color = "grey60") +
    geom_hline(yintercept = 1, linewidth = 0.5, linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = 2, linewidth = 0.5, linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = 3, linewidth = 0.5, linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = 4, linewidth = 0.5, linetype = "dashed", color = "grey80") +
    geom_rect(xmin = 28, xmax = 76,
              ymin = ymin, ymax = ymax,
              color = NA, fill = "green") + # PHD
    geom_text(x = (28 + 76) / 2,
              y = (ymin + ymax) / 2,
              color = "white",
              label = "PHD") +
    geom_rect(xmin = 163, xmax = 208,
              ymin = ymin, ymax = ymax,
              color = NA, fill = "blue") + # ZF , zinc finger domain
    geom_text(x = (163 + 208) / 2,
              y = (ymin + ymax) / 2,
              color = "white",
              label = "ZF") +
    geom_rect(xmin = 400, xmax = 636,
              ymin = ymin, ymax = ymax,
              color = NA, fill = "purple") + #zf-CpG_bind_C: CpG binding protein zinc finger C terminal domain
    geom_text(x = (400 + 636) / 2,
              y = (ymin + ymax) / 2,
              color = "white",
              label = "ZF-CpG binding") +
    geom_point(alpha = 0.85, size = 2) +
    # https://ggrepel.slowkow.com/articles/examples.html
    geom_text_repel(
        data = sdb, 
        mapping = aes(label = label),
        #ylim = c(3, 4.5),
        ylim = c(3.5, 4.5),
        color = "red", 
        size = 2.5, # 2
        angle = 0,
        nudge_y = 2,
        segment.curvature = -1e-20,
        fontface = "bold",
    ) +
    scale_fill_manual(values = cols, name = "Mutation types") +
    scale_color_manual(values = cols, name = "Mutation types") +
    theme_classic() +
    scale_x_continuous(
        expand = c(0.01, 0.01),
        limits = c(0, 656),
        breaks = c(seq(0, 650, 50), 28, 76),
        ) +
    scale_y_continuous(
        limits = c(-1, 5),
        breaks = c(seq(1, 4, 1)),
        labels = c(1:4),
        expand = c(0, 0),
        name = "Number of Mutations"
        ) +
    #coord_cartesian(ylim=c(0,1)) +
    theme(
        legend.position = "top",
        legend.direction = "horizontal",
        #axis.line.x = element_blank()
    )
# # table(mutdb$IMPACT)
pdf("./lolipop.pdf", width = 20, height = 3)
plot(p)
dev.off()
