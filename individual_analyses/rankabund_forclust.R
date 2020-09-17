#if (input$taxo_group_clust != "All") {
  
  otutab_tax <- otutab_sfsep %>% filter(Divisionlong %in% "Ochrophyta_Bacillariophyta")
# } else {
#   otutab_tax <- otutab_sfsep
# }

otutab_tax_sf0 <- otutab_tax[,c(as.character(sfnames$sf10.50))]



#if (input$propchoice == "selectgroup") {
  otutab_tax_sf1 <- sweep(otutab_tax_sf0, 2 , colSums(otutab_tax_sf0), FUN = "/")
#} else {
#   otutab_tax_sf1 <- otutab_tax_sf0
# }


otutab_tax_sf <- cbind.data.frame(otutab_tax_sf1, otutab_tax %>% select_if(negate(is.numeric)))
otutab_tax_sf[is.na(otutab_tax_sf)] <- 0
sjekk <- names(otutab_tax_sf)

otutab_tax_sf$total <- otutab_tax_sf %>% select_if(is.numeric) %>% rowSums()
#taxlevtab <- otutab_tax_sf %>% select(input$taxlevel_clustplot) %>% droplevels.data.frame()

#####
#Presence-absence for OTU count per taxonomic group
otutab_sfsep_pa <- otutab_sfsep
otutab_sfsep_pa[otutab_sfsep_pa>0] <- 1
otutab_tax_sf_pa <- otutab_tax_sf
otutab_tax_sf_pa[otutab_tax_sf_pa>0] <- 1

taxgroupspre <-  select(otutab_tax_sf, -total) %>%   group_by_(input$taxlevel_clustplot) %>% summarise_if(is.numeric, sum)

taxgroupspre_pa <- select(otutab_tax_sf_pa, -total) %>%   group_by_(input$taxlevel_clustplot) %>% summarise_if(is.numeric, sum)




#taxgroups <- melt(taxgroupspre)

#
limfun <- function(x) {
  ifelse(x>=0.05,1,0)
}

taxgroupspre_mat <- as.matrix(taxgroupspre[,-1, drop = FALSE])
taxgroupspre_bin <- apply(taxgroupspre_mat, 2, limfun)

#sjekk <- head(taxgroupspre_bin)

if (is.vector(taxgroupspre_bin)) {
  taxgroupspre_bin2 <- as.data.frame(as.list(taxgroupspre_bin))
} else {
  taxgroupspre_bin2 <- as.data.frame(taxgroupspre_bin)}


taxgroupspre_bin_sums <- taxgroupspre_bin2 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
taxgroupspre_bin_sums$taxgroups <- as.vector(taxgroupspre[[taxlevel]])
taxgroupspre_bin_sums_yes <- filter(taxgroupspre_bin_sums, rowsumm >0) %>% select(taxgroups)


taxgroups_select <- taxgroupspre %>% filter(get(taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)

taxgroups_other <- taxgroupspre %>% filter(!get(taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)

taxgroups_pa_select <- taxgroupspre_pa %>% filter(get(input$taxlevel_clustplot) %in% taxgroupspre_bin_sums_yes$taxgroups)

taxgroups_pa_other <- taxgroupspre_pa %>% filter(!get(input$taxlevel_clustplot) %in% taxgroupspre_bin_sums_yes$taxgroups)


if (dim(taxgroupspre_bin_sums_yes)[1] < dim(taxgroupspre_bin_sums)[1]) {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups,"Other")
  othersum <- colSums(taxgroups_other[,-1])
} else {
  Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups)
  othersum <- NULL
}



taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)


#
#  #
#  #
taxgroups_select31 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)
#taxgroups_select3_forcluster_sf <- taxgroups_select3 %>% select(sf)
taxgroups_select30 <- taxgroups_select31[,-1][,otuclust0$order]
taxgroups_select3 <- cbind(taxgroups_select31[,1], taxgroups_select30)
names(taxgroups_select3)[1] <- taxlevel
sjekk <- names(taxgroups_select3)

taxgroups_select3t <- t(taxgroups_select3)
taxgroups_select3t2 <- taxgroups_select3t[-1,]
taxgroups_select3t3 <- apply(taxgroups_select3t2, 2, as.numeric)
taxgroups_select3tdf0 <- as.data.frame(taxgroups_select3t3)
names(taxgroups_select3tdf0) <- Taxonomic_group
taxgroups_select3tdf0$sample <- factor(names(taxgroups_select3[,-1]), levels = names(taxgroups_select3[,-1]), order =T)
#sjekk <- "sjekkja"

taxgroups_select3tdfmelt <- melt(taxgroups_select3tdf0)

barplot_clust <- ggplot(taxgroups_select3tdfmelt, aes(x=sample, y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, sample, round(value*100,3))))+
  geom_bar(stat = "identity", position = "stack")+
  #facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(legend.text=element_text(size=8))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(legend.key.size = unit(0.1, "cm"))+
  coord_flip()
barplot_clust



barplot_clustly <- ggplotly(barplot_clust, tooltip = "text") %>%  layout(legend = list(x = 100, y = 0.5))

plotclust <- barplot_clustly
#####
#####################barplot for cluster end

######################rankotus for cluster start
#####
#####
# if (input$taxo_group_ra != "All") {
#
#   otutab_ra <- otutab_sfsep %>% filter(Divisionlong %in% input$taxo_group_ra)
# } else {
#   otutab_ra <- otutab_sfsep
# }


otutab_raclust_num0 <- otutab_tax_sf %>% select_if(is.numeric)
otutab_raclust_tax0 <- otutab_tax_sf %>% select_if(negate(is.numeric))
otutab_raclust_tax <- apply(otutab_raclust_tax0, 2, as.character)

#### must order otutab_raclust_num column according to cluster

otustand <- decostand(t(otutab_raclust_num0), method = "hellinger")
otudist <- vegdist(otustand, distance = "jaccard")
otuclust0 <- hclust(otudist, "average")


otutab_raclust_num <- otutab_raclust_num0[,otuclust0$order]




rownames(otutab_raclust_num) <- otutab_tax_sf$newOTUid_wgen
rownames(otutab_raclust_tax) <- otutab_tax_sf$newOTUid_wgen



rankmat <- matrix(, nrow = dim(otutab_raclust_num)[1], ncol = dim(otutab_raclust_num)[2])
for (i in c(1:dim(otutab_raclust_num)[2])) {
  rankmat[,i] <- sort(otutab_raclust_num[,i], decreasing = F)
}



#Transform into data frame
rankmatdf <- as.data.frame(rankmat)
names(rankmatdf) <- names(otutab_raclust_num)

rankotus <- matrix(, dim(otutab_raclust_num)[1], ncol = dim(otutab_raclust_num)[2])
for (i in c(1:dim(otutab_raclust_num)[2])) {
  rankotus[,i] <- rownames(otutab_raclust_num[order(otutab_raclust_num[,i]),, drop = FALSE])
}

notus <- 10 #number of OTUs to include in the "rank-abundance barplot"

#Select only the notu highest abundance values, and vectorize for ggplot.
#Remember that they are sorted in increasing order, so the highest values are at the bottom of the matrix
rankvals100 <- c(rankmat[c(((dim(rankmat)[1]-notus)+1):dim(rankmat)[1]),])

rankphylum <- matrix(, nrow = notus, ncol = dim(otutab_raclust_num)[2])
for (i in c(1:notus)) {
  for (j in c(1:dim(otutab_raclust_num)[2])){
    rankphylum[i,j] <- as.character(otutab_raclust_tax[rankotus[dim(otutab_tax_sf)[1]-(notus-i),j],"newOTUid_wgen"])
  }
}

#Vectorise this table for ggplot
rankphyls100 <- c(rankphylum)

#Group the "Divisionlongs" that do not have a representative among the 5 most abundant OTUs in any sample in "Other" category
rankphyls100.2 <- character()
for (i in 1:length(rankphyls100)){
  ifelse(rankphyls100[i] %in% unique(c(rankphylum[c((notus-4):notus),])), rankphyls100.2[i] <- rankphyls100[i], rankphyls100.2[i] <- "Other")
}

#Dummy matrix with sample names for ggplot
samples <- matrix(, nrow = notus, ncol = dim(otutab_raclust_num)[2])
for (i in c(1:dim(otutab_raclust_num)[2])) {
  samples[,i] <- rep(names(otutab_raclust_num)[i],notus)
}

#Vectorise sample names for ggplot
samplesvec <- c(samples)
samplesvecord <- factor(samplesvec, ordered = TRUE, levels = unique(samplesvec))

rep100 <- function(x) {
  rep(x, notus)
}

rankotus_notus <- as.matrix(rankotus[c((dim(otutab_raclust_num)[1]-(notus-1)):dim(otutab_raclust_num)[1]),]) #18841:18940
rankotus_notus2 <- rankotus_notus
#names(rankotus_notus2) <- names(OTUtab_merged_mean2)

rankotus_notus_vec <- c(rankotus_notus)


#Combine abundance values, taxonomic level and dummy variables into data frame
#if ()
seqra <- c(rep(c(1:notus),(nmp3+nmp4+nmp5)))

dfra <- cbind.data.frame(rankvals100, rankphyls100.2, samplesvecord, rankotus_notus_vec, seqra)


barplot_clust_ra <- ggplot(dfra, aes(x=samplesvecord, y = rankvals100, fill = rankphyls100.2, group=seqra, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", rankphyls100.2, samplesvec, round(rankvals100*100,3))))+
  geom_bar(stat = "identity", position = "stack")+
  #facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
  #theme(axis.text.x = element_blank())+
  #theme(axis.text.y = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(legend.text=element_text(size=8))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(legend.key.size = unit(0.1, "cm"))+
  coord_flip()
barplot_clust_ra
ggplotly(barplot_clust_ra)
