library(here)
pico <- readLines(here("data","pico_descnames.txt"))
three180 <- readLines(here("data","three180_descnames.txt"))

otutab_sfsep0 <- read.table(here("data", "OTUtab_nozooembr_minsize5_prop_wtax_wnewseqid_20191022_point_frep3.txt"), header = T, sep = "\t")[,-c(141:155)]
otutab_sfsep <- otutab_sfsep0[,c(5:8,1:4,9:122,127:130,123:126,131:151)]

otutab_sfsep_num <- otutab_sfsep %>%  select_if(is.numeric)
otutab_sfsep_tax <- otutab_sfsep %>%  select_if(negate(is.numeric))
otutab_sfsep_prop <- sweep(otutab_sfsep_num, 2 , colSums(otutab_sfsep_num), FUN = "/")

colSums(otutab_sfsep_prop)


#Make OTUtable no. 5

MP3_P1_1_3.200<- apply(otutab_sfsep_prop[, c("MP3_P1_1_3.10", "MP3_P1_1_10.50", "MP3_P1_1_50.200")], 1, mean)
MP3_P1_20_3.200<- apply(otutab_sfsep_prop[, c("MP3_P1_20_3.10", "MP3_P1_20_10.50", "MP3_P1_20_50.200")], 1, mean)
MP3_P1_417_3.200<- apply(otutab_sfsep_prop[, c("MP3_P1_417_3.10", "MP3_P1_417_10.50", "MP3_P1_417_50.200")], 1, mean)

MP3_P3_1_3.200<- apply(otutab_sfsep_prop[, c("MP3_P3_1_3.10", "MP3_P3_1_10.50", "MP3_P3_1_50.200")], 1, mean)
MP3_P3_15_3.200<- apply(otutab_sfsep_prop[, c("MP3_P3_15_3.10", "MP3_P3_15_10.50", "MP3_P3_15_50.200")], 1, mean)
MP3_P3_447_3.200<- apply(otutab_sfsep_prop[, c("MP3_P3_447_3.10", "MP3_P3_447_10.50", "MP3_P3_447_50.200")], 1, mean)


MP3_P4_1_3.200<- apply(otutab_sfsep_prop[, c("MP3_P4_1_3.10", "MP3_P4_1_10.50", "MP3_P4_1_50.200")], 1, mean)
MP3_P4_15_3.200<- apply(otutab_sfsep_prop[, c("MP3_P4_15_3.10", "MP3_P4_15_10.50", "MP3_P4_15_50.200")], 1, mean)
MP3_P4_500_3.200<- apply(otutab_sfsep_prop[, c("MP3_P4_500_3.10", "MP3_P4_500_10.50", "MP3_P4_500_50.200")], 1, mean)
MP3_P4_1000_3.200<- apply(otutab_sfsep_prop[, c("MP3_P4_1000_3.10", "MP3_P4_1000_10.50", "MP3_P4_1000_50.200")], 1, mean)

MP4_P5_1_3.200<- apply(otutab_sfsep_prop[, c("MP4_P5_1_3.10", "MP4_P5_1_10.50", "MP4_P5_1_50.200")], 1, mean)
MP4_P5_20_3.200<- apply(otutab_sfsep_prop[, c("MP4_P5_20_3.10", "MP4_P5_20_10.50", "MP4_P5_20_50.200")], 1, mean)
MP4_P5_213_3.200<- apply(otutab_sfsep_prop[, c("MP4_P5_213_3.10", "MP4_P5_213_10.50", "MP4_P5_213_50.200")], 1, mean)

MP4_P6_1_3.200<- apply(otutab_sfsep_prop[, c("MP4_P6_1_3.10", "MP4_P6_1_10.50", "MP4_P6_1_50.200")], 1, mean)
MP4_P6_24_3.200<- apply(otutab_sfsep_prop[, c("MP4_P6_24_3.10", "MP4_P6_24_10.50", "MP4_P6_24_50.200")], 1, mean)
MP4_P6_500_3.200<- apply(otutab_sfsep_prop[, c("MP4_P6_500_3.10", "MP4_P6_500_10.50", "MP4_P6_500_50.200")], 1, mean)
MP4_P6_1000_3.200<- apply(otutab_sfsep_prop[, c("MP4_P6_1000_3.10", "MP4_P6_1000_10.50", "MP4_P6_1000_50.200")], 1, mean)

MP4_P7_1_3.200<- apply(otutab_sfsep_prop[, c("MP4_P7_1_3.10", "MP4_P7_1_10.50", "MP4_P7_1_50.200")], 1, mean)
MP4_P7_25_3.200<- apply(otutab_sfsep_prop[, c("MP4_P7_25_3.10", "MP4_P7_25_10.50", "MP4_P7_25_50.200")], 1, mean)
MP4_P7_500_3.200<- apply(otutab_sfsep_prop[, c("MP4_P7_500_3.10", "MP4_P7_500_10.50", "MP4_P7_500_50.200")], 1, mean)
MP4_P7_1000_3.200<- apply(otutab_sfsep_prop[, c("MP4_P7_1000_3.10", "MP4_P7_1000_10.50", "MP4_P7_1000_50.200")], 1, mean)


MP5_MP5.2_20_3.200<- apply(otutab_sfsep_prop[, c("MP5_MP5.2_20_3.10", "MP5_MP5.2_20_10.50", "MP5_MP5.2_20_50.200")], 1, mean)


MP5_MP5.3_20_3.200<- apply(otutab_sfsep_prop[, c("MP5_MP5.3_20_3.10", "MP5_MP5.3_20_10.50", "MP5_MP5.3_20_50.200")], 1, mean)
MP5_MP5.3_300_3.200<- apply(otutab_sfsep_prop[, c("MP5_MP5.3_300_3.10", "MP5_MP5.3_300_10.50", "MP5_MP5.3_300_50.200")], 1, mean)


MP5_MP5.4_20_3.200<- apply(otutab_sfsep_prop[, c("MP5_MP5.4_20_3.10", "MP5_MP5.4_20_10.50", "MP5_MP5.4_20_50.200")], 1, mean)
MP5_MP5.4_1000_3.200<- apply(otutab_sfsep_prop[, c("MP5_MP5.4_1000_3.10", "MP5_MP5.4_1000_10.50", "MP5_MP5.4_1000_50.200")], 1, mean)



otutab_sfsep_3.200 <- as.data.frame(cbind(MP3_P1_1_3.200, MP3_P1_20_3.200, MP3_P1_417_3.200, MP3_P3_1_3.200, MP3_P3_15_3.200, MP3_P3_447_3.200, MP3_P4_1_3.200, MP3_P4_15_3.200, 
                                                MP3_P4_500_3.200, MP3_P4_1000_3.200, MP4_P5_1_3.200, MP4_P5_20_3.200, MP4_P5_213_3.200,
                                                MP4_P6_1_3.200, MP4_P6_24_3.200, MP4_P6_500_3.200, MP4_P6_1000_3.200,  MP4_P7_1_3.200, MP4_P7_25_3.200,
                                                MP4_P7_500_3.200, MP4_P7_1000_3.200,  MP5_MP5.2_20_3.200,  MP5_MP5.3_20_3.200, MP5_MP5.3_300_3.200,
                                                MP5_MP5.4_20_3.200, MP5_MP5.4_1000_3.200))

otutab_sfsep_3.200all <- as.data.frame(cbind(otutab_sfsep[,-c(141:151)], otutab_sfsep_3.200))
otutab_3.200all_mtax <- as.data.frame(cbind(otutab_sfsep_3.200all, otutab_sfsep_tax))

write.table(otutab_3.200all_mtax, file = "otutab_3.200all_mtax.txt", sep = "\t")

otutab_sf_tax_num <- otutab_sfsep_syndi_3.200all
rownames(otutab_sf_tax_num) <- otutab_sfsep_syndi$newOTUid_wgen
