# purpose: download mutation info of GDC
library(TCGAbiolinks)
projects <- getGDCprojects()$id


wd <- getwd()
for (project in projects) {
  setwd(wd)
  tryCatch(
    {
      f <- sprintf("./data/%s-mut.rds", project)
      if (file.exists(f)) next()
      query <- GDCquery(
        project = project,
        data.category = "Simple Nucleotide Variation",
        access = "open",
        # legacy = FALSE,
        data.type = "Masked Somatic Mutation"
        # workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
      )
      GDCdownload(
        query = query
      )
      dataPrep <- GDCprepare(
        query = query,
        save = TRUE
      )
      setwd(wd)
      saveRDS(dataPrep, f)
      # clinical
      f <- sprintf("./data/%s-clinical.rds", project)
      if (file.exists(f)) next()
      print(project)
      query <- GDCquery(
        project = project,
        data.category = "Clinical",
        data.type = "Clinical Supplement",
        data.format = "BCR Biotab"
      )
      setwd("~/.cache/TCGAbiolinks/")
      GDCdownload(query)
      dataPrep <- GDCprepare(query)
      setwd(wd)
      saveRDS(dataPrep, f)
    },
    error = function(cond) {
      message(cond)
    }
  )
  setwd(wd)
}
print("done")
