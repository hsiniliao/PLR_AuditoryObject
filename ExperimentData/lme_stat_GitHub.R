library("lme4")

path <- c("/Volumes/Install macOS Sierra/iMac_20200409/covert_auditory_attention/data/tables/")
filename <- c("E4Norm")
data_p <- read.csv(paste(path,filename,c("_p.csv"),sep=''), header=TRUE)
data_gLum <- read.csv(paste(path,filename,c("_gLum.csv"),sep=''), header=TRUE)

tArray_gLum <- vector()
tArray_lum <- vector()

for (i in 1:(length(data_p)-5))
{
  p <- data_p[[i]]
  gLum <- data_gLum[[i]]
  lum <- data_p[["lum"]]
  subj <- data_p[["subj"]]
  
  fitTmp <- lmer(p ~ gLum+lum+(1|subj))
  tArray_gLum <- append(tArray_gLum, coef(summary(fitTmp))[, "t value"][2])
  tArray_lum <- append(tArray_lum, coef(summary(fitTmp))[, "t value"][3])
}

write.csv(tArray_gLum, file = paste(path,filename,c("_p_tArray_gLum_GitHub.csv"),sep=''))
write.csv(tArray_lum, file = paste(path,filename,c("_p_tArray_lum_GitHub.csv"),sep=''))

