folder <- "C:/Thunderstorm test/Nik/GG1"
#folder <- "D:/Maarten/GG1"
condition <- "GG1 IR"

source("0 Functions.R")
#import data, import only once, next time just run load
SMOLR_IMPORT(folder = folder ,profile = "thunderstorm_merged",basename = "STORM_",
             condition=condition, channel = c("BRCA2","RAD51"),sep_chfiles = T, extension = ".csv",use_data.table = FALSE)

SMOLR_LOAD(folder = folder,statistics=F)
#head(localizations[[2]])

#load and apply ROI files to have only localization inside the nucleus
roi_files <- file.path(folder,paste0(seq(1:length(localizations)),".roi"))


for(i in 1:length(localizations)){
  localizations[[i]] <- IJROI_subset(localizations[[i]],file = roi_files[[i]],pxsize = 5)
}

save(localizations,file = file.path(folder,"localizations_ROI.Rdata"))
