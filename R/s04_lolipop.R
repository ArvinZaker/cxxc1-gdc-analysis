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
# sdb <- newMutDb[newMutDb$impact == "HIGH", ]
sdb <- newMutDb[phdDomRows | newMutDb$impact == "HIGH", ]
sdb$label <- apply(
    sdb, 1,
    function(x) {
        v <- paste0(unlist(strsplit(x["AA"], "/"))[1], x["pos"], unlist(strsplit(x["AA"], "/"))[2])
        v <- gsub(" ", "", v)
        if (grepl("NA", v)) {
            v <- gsub("NA", "", v)
            v <- paste0(v, substring(v, 1, 1))
        }
        return(v)
    }
)
# newMutDb$num <- newMutDb$num + 0.5
# newMutDb$num <- newMutDb$num * ifelse(phdDomRows, 1, -1)

# ymin <- -.5
# ymax <- .5
ymin <- -1
ymax <- 0
legsize <- 8

set.seed(123)
p <- ggplot(newMutDb, aes(x = pos, y = num, fill = type, color = type, group = type)) + # size = impact, ,  ,
    annotate("rect",
        xmin = 28, xmax = 76,
        ymin = 0, ymax = 4,
        alpha = 1, fill = "#F9F1E4"
    ) +
    # 1-3 bg dashed lines
    geom_hline(yintercept = 1, linewidth = 0.5, linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = 2, linewidth = 0.5, linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = 3, linewidth = 0.5, linetype = "dashed", color = "grey80") +
    # https://ggrepel.slowkow.com/articles/examples.html
    # annotations
    geom_text_repel(
        data = sdb,
        mapping = aes(label = label),
        size = legsize / 2,
        ylim = c(3.5, 4),
        alpha = .75,
        color = "#373640",
        size = 2.5, # 2
        angle = 0,
        nudge_y = 2,
        segment.curvature = -1e-20,
        fontface = "bold",
    ) +
    # the lines that go to dots
    geom_segment(
        mapping = aes(xend = pos, yend = -.5),
        color = "black",
        size = .25, alpha = .85
    ) +
    geom_point(alpha = 1, size = 2) +
    # baseline of the prot
    geom_hline(yintercept = (ymax + ymin) / 2, linewidth = 3, color = "grey60") +
    geom_rect(
        xmin = 28, xmax = 76,
        ymin = ymin, ymax = ymax,
        color = NA,
        fill = "#d95f02"
    ) + # PHD 3"green" "#C4A484"
    geom_text(
        x = (28 + 76) / 2,
        y = (ymin + ymax) / 2,
        color = "white",
        size = legsize,
        label = "PHD"
    ) +
    geom_rect(
        xmin = 163, xmax = 208,
        ymin = ymin, ymax = ymax,
        color = NA,
        fill = "#1b9e77"
    ) + # ZF , zinc finger domain
    geom_text(
        x = (163 + 208) / 2,
        y = (ymin + ymax) / 2,
        color = "white",
        size = legsize,
        label = "ZF"
    ) +
    geom_rect(
        xmin = 400, xmax = 636,
        ymin = ymin, ymax = ymax,
        color = NA,
        fill = "#7570b3"
    ) + # zf-CpG_bind_C: CpG binding protein zinc finger C terminal domain
    geom_text(
        x = (400 + 636) / 2,
        y = (ymin + ymax) / 2,
        color = "white",
        size = legsize,
        label = "ZF-CpG binding"
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
        # limits = c(-1, 4),
        limits = c(-1, 5),
        breaks = 1:3,
        labels = 1:3,
        expand = c(0, 0),
        name = "Number of Mutations"
    ) +
    guides(fill = "none") +
    xlab("") +
    theme(
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(size = legsize * 2),
        axis.text.y = element_text(size = legsize * 2),
        axis.title.x = element_text(size = legsize * 2.5),
        axis.title.y = element_text(size = legsize * 2.5),
    )
pdf("./figures/lolipop.pdf", width = 20, height = 20 / 4)
plot(p)
dev.off()
