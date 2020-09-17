library(shiny)
library(here)
library(ggplot2)
library(tidyverse)
library(vegan)
library(cowplot)
library(ggpubr)
library(mvabund)
library(RColorBrewer)
library(reshape2)
library(purrr)
library(dplyr)
library(plotly)
library(indicspecies)
library(adespatial)
library(fossil)

bpps <- 0.6
otutab_sfsep0 <- read.table(here("data", "OTUtab_nozooembr_minsize5_prop_wtax_wnewseqid_20191022_point_frep2.txt"), header = T, sep = "\t")[,-c(141:155)]
#otutab_sfsep0 <- read.table(here("data", "OTUtab_nozooembr_minsize5_prop_wtax_wnewseqid_14032019_point2.txt"), header = T, sep = "\t")[,-c(141:155)]
#syndi <- otutab_sfsep %>% filter(Divisionlong %in% "Dinoflagellata_Syndiniales")
#syndi$total = syndi %>% select_if(is.numeric) %>% rowSums()
# 
#syndi <- syndi %>% arrange(desc(total))
# view(head(syndi, n= 50))
# write.table(syndi[,c(143:152)], "syndi.txt")

dim(otutab_sfsep0)
otutab_sfsep <- otutab_sfsep0[,c(5:8,1:4,9:122,127:130,123:126,131:151)]
# otutab_sfsep_strings1 <- otutab_sfsep %>% select_if(negate(is.numeric))
# otutab_sfsep_strings2 <- otutab_sfsep_strings1 %>% mutate_all(list(~str_replace_all(., "Dino-Group", "DG")))
# otutab_sfsep_strings3 <- otutab_sfsep_strings2 %>% mutate_all(list(~str_replace_all(., "unclassified", "uncl")))
# otutab_sfsep_num <- otutab_sfsep %>% select_if(is.numeric)
# 
# otutab_sfsep <- cbind(otutab_sfsep_num, otutab_sfsep_strings3)
# 
# names(otutab_sfsep) <- names(otutab_sfsep1)

#names(otutab_sfsep)
rm(otutab_sfsep0)
sample_ordersfsep <- tibble(sample = names(otutab_sfsep)[c(1:140)])
pico <- readLines(here("data","pico_descnames.txt"))
three10 <- readLines(here("data","three10_descnames.txt"))
three180 <- readLines(here("data","three180_descnames.txt"))
ten50 <- readLines(here("data","ten50_descnames.txt"))
fifty200 <- readLines(here("data","fifty200_descnames.txt"))
net_all <- readLines(here("data","net_all_descnames.txt"))
sfnames <- list("sf0.4.3" = pico, "sf3.10" = three10, "sf10.50" = ten50, "sf50.200" = fifty200, "sf3.180" = three180)
sfs <- names(sfnames)

npico <- length(pico)
nthree10 <- length(three10)
nten50 <- length(ten50)
nfifty200 <- length(fifty200)
nthree180 <- length(three180)

nmp1 <- 8
nmp2 <- 10
nmp3 <- 10
nmp4 <- 11
nmp5 <- 5

mp1names <- c(pico[c(1:nmp1)], three180[c(1:nmp1)])
mp2names <- c(pico[c((nmp1+1):(nmp1+nmp2))], three180[c((nmp1+1):(nmp1+nmp2))])
mp3names <- c(pico[c((nmp1+nmp2+1):(nmp1+nmp2+nmp3))], three10[c(1:nmp3)], ten50[c(1:nmp3)], fifty200[c(1:nmp3)])
mp4names <- c(pico[c((nmp1+nmp2+nmp3+1):(nmp1+nmp2+nmp3+nmp4))], three10[c((nmp3+1):(nmp3+nmp4))], ten50[c((nmp3+1):(nmp3+nmp4))], fifty200[c((nmp3+1):(nmp3+nmp4))])
mp5names <- c(pico[c((nmp1+nmp2+nmp3+nmp4+1):(nmp1+nmp2+nmp3+nmp4+nmp5))], three10[c((nmp3+nmp4+1):(nmp3+nmp4+nmp5))], ten50[c((nmp3+nmp4+1):(nmp3+nmp4+nmp5))], fifty200[c((nmp3+nmp4+1):(nmp3+nmp4+nmp5))])


meta_table <- read.table(here("data", "MP_CB_profiles_20170802_wsampnames.txt"), header = T, sep = "\t", fill = NA)

depthbin0.4.3 <- meta_table %>% filter(sf0.4.3 %in% sfnames$sf0.4.3) %>% pull(depthbin) %>% as.character()
depthbin3.10 <- meta_table %>% filter(sf3.10 %in% sfnames$sf3.10) %>% pull(depthbin) %>% as.character()
depthbin10.50 <- meta_table %>% filter(sf10.50 %in% sfnames$sf10.50) %>% pull(depthbin) %>% as.character()
depthbin50.200 <- meta_table %>% filter(sf50.200 %in% sfnames$sf50.200) %>% pull(depthbin) %>% as.character()
depthbin3.180 <- meta_table %>% filter(sf3.180 %in% sfnames$sf3.180) %>% pull(depthbin) %>% as.character()

depthbinsfsep <- c(depthbin0.4.3, depthbin3.10, depthbin10.50, depthbin50.200, depthbin3.180)

samples0.4.3 <- meta_table %>% filter(sf0.4.3 %in% sfnames$sf0.4.3) %>% pull(sf0.4.3) %>% as.character()
samples3.10 <- meta_table %>% filter(sf3.10 %in% sfnames$sf3.10) %>% pull(sf3.10) %>% as.character()
samples10.50 <- meta_table %>% filter(sf10.50 %in% sfnames$sf10.50) %>% pull(sf10.50) %>% as.character()
samples50.200 <- meta_table %>% filter(sf50.200 %in% sfnames$sf50.200) %>% pull(sf50.200) %>% as.character()
samples3.180 <- meta_table %>% filter(sf3.180 %in% sfnames$sf3.180) %>% pull(sf3.180) %>% as.character()

samplessfsep <- c(samples0.4.3, samples3.10, samples10.50, samples50.200, samples3.180)

depthbinsfsep_tb0 <- tibble(depthbin = depthbinsfsep, sample = samplessfsep)

depthbinsfsep_tb <- left_join(sample_ordersfsep, depthbinsfsep_tb0, by = "sample")

depthlab0.4.3 <- meta_table %>% filter(sf0.4.3 %in% sfnames$sf0.4.3) %>% pull(DEPTH_M) %>% as.character()
depthlab3.10 <- meta_table %>% filter(sf3.10 %in% sfnames$sf3.10) %>% pull(DEPTH_M) %>% as.character()
depthlab10.50 <- meta_table %>% filter(sf10.50 %in% sfnames$sf10.50) %>% pull(DEPTH_M) %>% as.character()
depthlab50.200 <- meta_table %>% filter(sf50.200 %in% sfnames$sf50.200) %>% pull(DEPTH_M) %>% as.character()
depthlab3.180 <- meta_table %>% filter(sf3.180 %in% sfnames$sf3.180) %>% pull(DEPTH_M) %>% as.character()

depthlabsfsep <-  c(depthlab0.4.3, depthlab3.10, depthlab10.50, depthlab50.200, depthlab3.180)                    

depthlabsfsep_tb0 <- tibble(depthlab = depthlabsfsep, sample = samplessfsep)

depthlabsfsep_tb <- left_join(sample_ordersfsep, depthlabsfsep_tb0, by = "sample")

cruise0.4.3 <- meta_table %>% filter(sf0.4.3 %in% sfnames$sf0.4.3) %>% pull(Cruise) %>% as.character()
cruise3.10 <- meta_table %>% filter(sf3.10 %in% sfnames$sf3.10) %>% pull(Cruise) %>% as.character()
cruise10.50 <- meta_table %>% filter(sf10.50 %in% sfnames$sf10.50) %>% pull(Cruise) %>% as.character()
cruise50.200 <- meta_table %>% filter(sf50.200 %in% sfnames$sf50.200) %>% pull(Cruise) %>% as.character()
cruise3.180 <- meta_table %>% filter(sf3.180 %in% sfnames$sf3.180) %>% pull(Cruise) %>% as.character()

cruisesfsep <-  c(cruise0.4.3, cruise3.10, cruise10.50, cruise50.200, cruise3.180)                    

cruisesfsep_tb0 <- tibble(cruise = cruisesfsep, sample = samplessfsep)

cruisesfsep_tb <- left_join(sample_ordersfsep, cruisesfsep_tb0, by = "sample")

Station0.4.3 <- meta_table %>% filter(sf0.4.3 %in% sfnames$sf0.4.3) %>% pull(Station) %>% as.character()
Station3.10 <- meta_table %>% filter(sf3.10 %in% sfnames$sf3.10) %>% pull(Station) %>% as.character()
Station10.50 <- meta_table %>% filter(sf10.50 %in% sfnames$sf10.50) %>% pull(Station) %>% as.character()
Station50.200 <- meta_table %>% filter(sf50.200 %in% sfnames$sf50.200) %>% pull(Station) %>% as.character()
Station3.180 <- meta_table %>% filter(sf3.180 %in% sfnames$sf3.180) %>% pull(Station) %>% as.character()

Stationsfsep <-  c(Station0.4.3, Station3.10, Station10.50, Station50.200, Station3.180)                    

Stationsfsep_tb0 <- tibble(Station = Stationsfsep, sample = samplessfsep)

Stationsfsep_tb <- left_join(sample_ordersfsep, Stationsfsep_tb0, by = "sample")

StationDepthsfsep_tb <- paste(Stationsfsep_tb$Station, depthlabsfsep_tb$depthlab, sep = "_")
StationDepthsfsep_tb <- factor(StationDepthsfsep_tb, ordered = T, levels = unique(StationDepthsfsep_tb))

####
#Chloroplast
chltab <- read.table(here("data", "chltab_prop_mtax.txt"), header = T, sep = "\t", row.names = 1)
chlmeta <- read.table(here("data", "chl_samplesmeta_dnavsrna.txt"), header = T, sep = "\t", row.names = 1)
chlmeta$sample <- rownames(chlmeta)

#Chl meta tables

StationDepths_tb_chl <- paste(chlmeta$station, chlmeta$depth, sep = "_")
StationDepths_tb_chl <- factor(StationDepths_tb_chl, ordered = T, levels = unique(StationDepths_tb_chl))




source(here("src", "nmdsplot.R"))
source(here("src", "mp_nmds.R"))
source(here("src", "taxa_plot.R"))




sf <- sfnames[["sf3.10"]]
otutab_sf <- otutab_sfsep[c(as.character(sf), "Divisionlong", "newOTUid_wgen")]
sampnames_meta_sf <- meta_table["sf3.10"]
meta_sf_all <- meta_table[match(colnames(otutab_sf %>% select(-Divisionlong, -newOTUid_wgen)),as.character(sampnames_meta_sf[,1])),]

meta_sf <- droplevels(meta_sf_all)

#if (input$taxo_group != "All") {
  
  otutab_sf_tax <- otutab_sf %>% filter(Divisionlong %in% "Dinoflagellata_Syndiniales")
#} else {
 # otutab_sf_tax <- otutab_sf
#}



otutab_sf_tax_num <- otutab_sf_tax %>% select(-Divisionlong, -newOTUid_wgen)
#hmm <- names(otutab_sf_tax) #"Divisionlong is last column
numotus <- length(which(rowSums(otutab_sf_tax_num)>0))
rownames(otutab_sf_tax) <- otutab_sf_tax$newOTUid_wgen




# if (input$taxo_group == "Stramenopiles_X" && input$sizefract1 == "sf3.10") {
#   nmds_tax <- mp_nmds(otutab_sf_tax[,-26] %>% select(-Divisionlong, -newOTUid_wgen))
#   nmds_table <- data.frame(cbind(nmds_tax$points[,1], nmds_tax$points[,2], meta_sf[-26,]))
# } else {
  nmds_tax <- mp_nmds(otutab_sf_tax %>% select(-Divisionlong, -newOTUid_wgen)) #function mp_nmds transposes otu table
  nmds_table <- data.frame(cbind(nmds_tax$points[,1], nmds_tax$points[,2], meta_sf))
#}




names(nmds_table)[c(1,2)] <- c("nmds_axis1", "nmds_axis2")

#Hierarchical clustering of samples
braydist <- vegdist(t(otutab_sf_tax_num), distance = "bray")
otuclust <- hclust(braydist, "average")
plotclust <- plot(otuclust)
plot(otuclust)
huh <- rect.hclust(otuclust, k = 4)
grp <- cutree(otuclust, k =4)
#plot(otuclust)
#rect.hclust(otuclust, k = 10)
#grp <- cutree(otuclust, k = 10)
#plot(ord, display = "sites")
ok <- plot(nmds_tax, display = "sites")
ordih <- ordihull(nmds_tax, grp, col = "red")

numclust <- 4
otustand <- decostand(t(otutab_sf_tax_num), method = "hellinger")
otudist <- vegdist(otustand, distance = "jaccard")
otuclust <- hclust(otudist, "average")
plot(otuclust)
rect.hclust(otuclust, k = 4)
otuprop <- sweep(otutab_sf_tax_num, 2 , colSums(otutab_sf_tax_num), FUN = "/")
otupropt <- data.frame(t(otuprop))
names(otupropt) <- otutab_sf_tax$newOTUid_wgen
head(otupropt)
#otuprop <- decostand(t(otutab_sf_tax_num), method = "total")
grp <- cutree(otuclust, k =4)


iva.clust <- multipatt(otupropt, grp, func = "r.g", max.order = 4, control = how(nperm = 999))
res <- data.frame(otu = rownames(iva.clust$str), pval = iva.clust$sign$p.value) %>% filter(pval<= 0.05)

mulpat_tab <- iva.clust$sign %>% mutate(newOTUid_wgen = rownames(.)) %>% filter(p.value <= 0.05)
mulpat_tab$Family <- otutab_sfsep %>% filter(newOTUid_wgen %in% mulpat_tab$newOTUid_wgen) %>% pull(Family)
mulpat_forbp <- mulpat_tab %>% group_by(Family) %>% count(index)
mulpat_tab %>% count(index)
mulpat_tab %>% count(index) %>% select(n) %>% sum
ggplotly(ggplot(mulpat_forbp, aes(x=index, y=n, fill = Family, text = sprintf("Taxon: %s<br>n OTUs: %s", Family, n)))+
           geom_bar(stat = "identity", position = "stack")  )


otuclust0 <- 
hclust(braydist, "average")
plot(otuclust0)
otuclust1 <- as.dendrogram(otuclust0)

oturect <- rect.dendrogram(otuclust1, k = numclus)

# ordihl <- list()
#for (i in c(1:length(ordih))) {
# ordihl[i] <- as.data.frame(cbind(ordih$i[,1], ordih[i][,2]))

#}
rownames(ordih$`1`) <- c()
ordih1 <- as.data.frame(cbind(ordih$'1'[,1], ordih$'1'[,2]))

if (input$sizefract1 == "sf0.4.3") {
  Month <-c("January", "March", "May", "August", "November")
} else { if (input$sizefract1 == "sf3.180") {
  Month <- c("January", "March")
} else {
  Month <- c("May", "August", "November")}}


#mvaotutab <- mvabund(t)




difftab = as.data.frame(cbind(Month, admp, anglmmp))
names(difftab) <- c("Month", "p-value adon", "p-value glm")

###Plotting NMDS
plot_nmds <- nmdsplot(nmds_table)+
  xlim(min(nmds_table$nmds_axis1)-.1,max(nmds_table$nmds_axis1)+.1)+
  ylim(min(nmds_table$nmds_axis2)-.1,max(nmds_table$nmds_axis2)+.1)#+
#geom_path(aes(x = V1, y = V2), data = ordih1)

#coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
nmdsplotly <- ggplotly(plot_nmds, tooltip = "text")





# Optimal number of clusters according as per indicator species
# analysis (IndVal, Dufrene-Legendre; package: labdsv)
otutabnum <-  otutab_sf_tax %>% select(-Divisionlong, -newOTUid_wgen) * 10000
otutabnumr <- round(otutabnum,0)
otutabnumr <- otutabnumr[rowSums(otutabnumr !=0) > 0,]
dim(otutabnumr)
otutabnumt <- t(otutabnumr)

IndVal <- numeric(nrow(otutabnumt))
ng <- numeric(nrow(otutabnumt))
for (k in 2:10)
{
  iva <- indval(otutabnumt, cutree(otuclust, k = k), func = "r.g",
                   control = how(nperm = 999))
  gr <- factor(iva$maxcls[iva$pval <= 0.05])
  ng[k] <- length(levels(gr)) / k
  iv <- iva$indcls[iva$pval <= 0.05]
  IndVal[k] <- sum(iv)
}
k.best <- which.max(IndVal[ng == 1]) + 1
col3 <- rep(1, nrow(otutabnumt))
col3[ng == 1] <- 3

dev.new(
  title = "IndVal-based search for optimal number of clusters",
  width = 12,
  height = 6,
  noRStudioGD = F
)
par(mfrow = c(1, 2))
plot(
  1:nrow(otutabnumt),
  IndVal,
  type = "h",
  main = "IndVal-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "IndVal sum",
  col = col3
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(
  which.max(IndVal),
  max(IndVal),
  pch = 16,
  col = "red",
  cex = 1.5
)
text(28, 15.7, "a", cex = 1.8)

plot(
  1:nrow(otutabnumt),
  ng,
  type = "h",
  xlab = "k (number of clusters)",
  ylab = "Ratio",
  main = "Proportion of clusters with significant indicator species",
  col = col3
)
axis(1,
     k.best,
     paste("optimum", k.best, sep = "\n"),
     col = "red",
     font = 2,
     col.axis = "red")
points(k.best,
       max(ng),
       pch = 16,
       col = "red",
       cex = 1.5)
text(28, 0.98, "b", cex = 1.8)

