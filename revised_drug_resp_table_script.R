cells <- c("MFM-223","CAL-51","CAL-120","HDQ-P1","DU-4475","Hs-578-T","HCC38","CAL-85-1","MDA-MB-468","MDA-MB-436","MDA-MB-231","HCC1187","HCC1143","HCC1395","HCC2157","BT-20","HCC1937")
# gdsc1 <- read.csv("D:/work/drug_resp/GDSC1_fitted_dose_response_25Feb20.txt",stringsAsFactors = F,sep = "\t")
# gdsc2 <-read.csv("D:/work/drug_resp/GDSC2_fitted_dose_response_25Feb20.txt",stringsAsFactors = F,sep = "\t")
gdsc2 <-read.csv("D:/work/drug_resp/GDSC2.txt",stringsAsFactors = F,sep = "\t",header = T)
# overall_gdsc <- rbind(gdsc1,gdsc2)
overall_gdsc <- gdsc2
tnbc_gdsc <- overall_gdsc[overall_gdsc$CELL_LINE_NAME %in% cells,]
# cells <- as.character(unique(gdsc$CELL_LINE_NAME))
tnbc_gdsc$tmp <- exp(tnbc_gdsc$LN_IC50)

tnbc_gdsc$resp <- ifelse(tnbc_gdsc$tmp > tnbc_gdsc$MAX_CONC, "R", "S")

# for i in nrow(tnbc_gdsc){
#     ifelse(tnbc_gdsc[tnbc_gdsc$CELL_LINE_NAME[[i,]] == tnbc_gdsc$MAX_CONC, "R", "S"
# ,}
library(reshape)
tmp <- tnbc_gdsc[,c("CELL_LINE_NAME","DRUG_NAME","resp")]
tmp <- melt(tmp,id = c("CELL_LINE_NAME","DRUG_NAME"))
temp <- cast(tmp, CELL_LINE_NAME+DRUG_NAME~value)
# remove the ambigousity and assign the final response then remove the NA drugs from the table
fnl_resp <- ifelse(temp$R == temp$S, NA, ifelse(temp$R > temp$S, "R", "S"))
temp$resp <- fnl_resp
temp<- na.omit(temp)
resp_tnbc <- temp[,c("CELL_LINE_NAME","DRUG_NAME","resp")]



library(tidyr) ##### Long to wide data format   
drug_resp <- spread(resp_tnbc,DRUG_NAME, resp)
temp <- drug_resp[,complete.cases(t(drug_resp))]

# resp_tnbc <- as.data.frame(t(trans))
rownames(temp) <- temp[,1]
temp <- temp[,-1]
resp_tnbc <- temp
drugs <- as.character(names(resp_tnbc))
# tmp <- as.data.frame(t(resp_tnbc))
# rm(drug_resp)
# rm(tmp)
# rm(gdsc2)
# rm(trans)
# rm(fnl_resp)

########
# change "R" to "0" in resp
resp_bin <- as.data.frame(sapply(resp_tnbc, function(x) {
  ifelse(x == "R", 0, 1)}))
rownames(resp_bin) <- rownames(resp_tnbc)
colnames(resp_bin) <- colnames(resp_tnbc)
temp <- resp_bin[,colSums(resp_bin) >= 3 & colSums(resp_bin) <= 14]
final_resp <- as.data.frame(sapply(temp, function(x) {
  ifelse(x == 0, "R", "S")}))
rownames(final_resp) <- rownames(temp)
filtered_drugs <- names(final_resp)
write.table(final_resp,row.names = T,col.names = T,quote = F,sep = "\t", "revised_drug_response.txt")
targets <- tnbc_gdsc[tnbc_gdsc$DRUG_NAME %in% names(final_resp),]
targets  <- spread(targets,DRUG_NAME, resp)