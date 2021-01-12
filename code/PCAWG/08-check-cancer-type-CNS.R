library(sigminer)

cancer_type_files <- list.files("data/cancer_types/", pattern = "PCAWG_CN176", full.names = TRUE)
pcawg <- lapply(cancer_type_files, readRDS)

bp_show_survey2(pcawg[[1]]) # 11
bp_show_survey2(pcawg[[2]]) # 14
bp_show_survey2(pcawg[[3]]) # 13
bp_show_survey2(pcawg[[4]]) # 8
bp_show_survey2(pcawg[[5]]) # 14
bp_show_survey2(pcawg[[6]]) # 11
bp_show_survey2(pcawg[[7]]) # 12
bp_show_survey2(pcawg[[8]]) # 8
bp_show_survey2(pcawg[[9]]) # 7
bp_show_survey2(pcawg[[10]]) # 5
bp_show_survey2(pcawg[[11]]) # 13
bp_show_survey2(pcawg[[12]]) # 9
bp_show_survey2(pcawg[[13]]) # 15
bp_show_survey2(pcawg[[14]]) # 7
bp_show_survey2(pcawg[[15]]) # 11
bp_show_survey2(pcawg[[16]]) # 12
bp_show_survey2(pcawg[[17]]) # 12
bp_show_survey2(pcawg[[18]]) # 17
bp_show_survey2(pcawg[[19]]) # 12
bp_show_survey2(pcawg[[20]]) # 6
bp_show_survey2(pcawg[[21]]) # 3
bp_show_survey2(pcawg[[22]]) # 7
bp_show_survey2(pcawg[[23]]) # 11
bp_show_survey2(pcawg[[24]]) # 19
bp_show_survey2(pcawg[[25]]) # 8
bp_show_survey2(pcawg[[26]]) # 10
bp_show_survey2(pcawg[[27]]) # 11
bp_show_survey2(pcawg[[28]]) # 9
bp_show_survey2(pcawg[[29]]) # 9
bp_show_survey2(pcawg[[30]]) # 15
bp_show_survey2(pcawg[[31]]) # 6
bp_show_survey2(pcawg[[32]]) # 17

data.frame(
  cancer_type = sub(".*_CN176_([^_]+).rds", "\\1", basename(cancer_type_files)),
  signum = c(11, 14, 13, 8, 14, 11, 12, 8,
             7, 5, 13, 9, 15, 7, 11, 12,
             12, 17, 12, 6, 3, 7, 11, 19,
             8, 10, 11, 9, 9, 15, 6, 17)
)


