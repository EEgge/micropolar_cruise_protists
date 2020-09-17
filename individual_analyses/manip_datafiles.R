library(plyr); library(dplyr)
names(mulpat_tab_list[[1]])
names(mulpat_tab_list[[1]])[9] <- "newOTUid_wgen"

syndimmulitpatt <- left_join(syndi, mulpat_tab_list[[1]], by = "newOTUid_wgen")
# Warning message:
#   Column `newOTUid_wgen` joining factor and character vector, coercing into character vector 
head(syndimmulitpatt)
dim(syndimmulitpatt)



indextab <- data.frame(ses = c("Jan", "Mar", "May", "Aug", "Nov", "Jan+Mar", 
                               "Jan+May", "Jan+Aug", "Jan+Nov", "Mar+May", "Mar+Aug",
                               "Mar+Nov", "May+Aug", "May+Nov", "Aug+Nov", 
                               "Jan+Mar+May", "Jan+Mar+Aug", "Jan+Mar+Nov",
                               "Jan+May+Aug", "Jan+May+Nov", "Jan+Aug+Nov",
                               "Mar+May+Aug", "Mar+May+Nov", "Mar+Aug+Nov",
                               "May+Aug+Nov", "Jan+Mar+May+Aug", "Jan+May+Aug+Nov",
                                "Jan+Mar+Aug+Nov", "Jan+Mar+May+Nov", "Mar+May+Aug+Nov"),
                       index = c(1:30))
indextab

syndimultipat_epi <- left_join(syndimmulitpatt, indextab, by = "index")

write.table(syndimultipat_epi[,c(142:162)], "syndimultipatt_epi.txt")

write.table(syndimmulitpatt[,c(144:161)], "syndimultpatt_043_meso.txt")