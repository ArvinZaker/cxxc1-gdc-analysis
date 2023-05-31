# purpose: download RNASeq of TCGA PRAD from TCGA biolinker
library(TCGAbiolinks)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(viridis)
library(rstatix)

projects <- getGDCprojects()$id

mutdb <- data.frame()
clindb <- data.frame()
for (project in projects) {
    fmut <- sprintf("./data/%s-mut.rds", project)
    fclin <- sprintf("./data/%s-clinical.rds", project)
    if (file.exists(fmut)) {
        print(project)
        mut <- readRDS(fmut)
        mut <- mut[mut$Hugo_Symbol %in% "CXXC1", ]
        if (nrow(mut) == 0) next()
        mut <- mut[
            ,
            c(
                "IMPACT", "VARIANT_CLASS", "CLIN_SIG", # "SOMATIC", "EXON",
                "INTRON", "Protein_position", "Amino_acids",
                "t_depth", "t_ref_count", "t_alt_count",
                "Exon_Number", "Mutation_Status",
                "Hugo_Symbol", "Variant_Classification", "Variant_Type",
                "Tumor_Sample_Barcode"
            )
        ]
        mut$study <- project
        mutdb <- rbind(mut, mutdb)
    }
    # clinical data
    if (file.exists(fclin)) {
        clin <- readRDS(fclin)
        size <- max(sapply(clin, nrow))
        t <- data.frame(study = project, size = size)
        clindb <- rbind(t, clindb)
    }
}

saveRDS(mutdb, "./results/mutdb.rds")
saveRDS(clindb, "./results/clindb.rds")
print("done")
