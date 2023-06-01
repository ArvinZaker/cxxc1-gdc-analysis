#!/bin/sh
Rscript ./R/s01_gdc_download.R
Rscript ./R/s02_process.R
Rscript ./R/s03_bar.R
Rscript ./R/s04_lolipop.R
