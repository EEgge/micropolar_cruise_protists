library(pvclust)
library(SpadeR)
library(grid)
library(pals)

dinocol <- read.table(here("data","colors_dinocil.txt"), header = T, sep = "\t")
dinocol_ra <- read.table(here("data","col_dinocil_ra.txt"), header = T, sep = "\t")
dinocolvec <- as.character(dinocol$color)
names(dinocolvec) <- dinocol$Order
dinocolvec_ra <- as.character(dinocol_ra$color)
names(dinocolvec_ra) <- dinocol_ra$OTU

leg_plot <- ggplot(dinocol, aes(x = Order, fill = Order))+
  geom_bar()+
  scale_fill_manual(values = dinocolvec)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 0.7)),
         color = guide_legend(override.aes = list(size = 0.7)))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"))
leg_plot
dino_legnd_10200 <- get_legend(leg_plot)

leg_plot_ra <- ggplot(dinocol_ra, aes(x = OTU, fill = OTU))+
  geom_bar()+
  scale_fill_manual(values = dinocolvec_ra)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 0.7)),
         color = guide_legend(override.aes = list(size = 0.7)))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "lines"))
leg_plot_ra
dino_legnd_10200_ra <- get_legend(leg_plot_ra)

###Bacillariophyta legend
bacillcol <- read.table(here("data","col_bacill.txt"), header = T, sep = "\t", comment.char = "")
bacillcol_ra <- read.table(here("data","col_bacill_ra.txt"), header = T, sep = "\t", comment.char = "")
bacillcolvec <- as.character(bacillcol$color)
names(bacillcolvec) <- bacillcol$Genus
bacillcolvec_ra <- as.character(bacillcol_ra$color)
names(bacillcolvec_ra) <- bacillcol_ra$OTU

bacill_leg_plot <- ggplot(bacillcol, aes(x = Genus, fill = Genus))+
  geom_bar()+
  scale_fill_manual(values = bacillcolvec)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 0.7)),
         color = guide_legend(override.aes = list(size = 0.7)),
         fill = guide_legend(nrow = 4, byrow = F))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"))
bacill_leg_plot
bacill_legnd_10200 <- get_legend(bacill_leg_plot)


bacill_leg_plot_ra <- ggplot(bacillcol_ra, aes(x = OTU, fill = OTU))+
  geom_bar()+
  scale_fill_manual(values = bacillcolvec_ra)+
  theme(legend.position = "bottom")+
  guides(shape = guide_legend(override.aes = list(size = 0.7)),
         color = guide_legend(override.aes = list(size = 0.7)),
         fill = guide_legend(nrow = 4, byrow = F))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "lines"))
bacill_leg_plot_ra
bacill_legnd_10200_ra <- get_legend(bacill_leg_plot_ra)

#####
  

#Make figures for paper¨
  currentvec <- bacillcolvec
  currentvec_ra <- bacillcolvec_ra
  #currentylab <- "% of Dinophyceae+Ciliophora reads"
  currentylab <- "% of Bacillariophyta reads"
  sf <- "sf50.200"
  samplessf <- sfnames[[sf]]
  sjekk <- sf
  taxlevel <- "Genus"
  taxo_group <- c("Ochrophyta_Bacillariophyta")
  numclus <- 2
  otutab_sf <- otutab_sfsep[c(as.character(samplessf), "Divisionlong", "newOTUid_wgen")]
  sampnames_meta_sf <- meta_table[sf]
  meta_sf_all <- meta_table[match(colnames(otutab_sf %>% select(-Divisionlong, -newOTUid_wgen)),as.character(sampnames_meta_sf[,1])),]
  
  meta_sf <- droplevels(meta_sf_all)
  
  
  otutab_sf_tax <- otutab_sf %>% filter(Divisionlong %in% taxo_group)

  # 
  # 
  otutab_sf_tax_num <- otutab_sf_tax %>% select(-Divisionlong, -newOTUid_wgen)
  #hmm <- names(otutab_sf_tax) #"Divisionlong is last column
  numotus <- length(which(rowSums(otutab_sf_tax_num)>0))
  rownames(otutab_sf_tax) <- otutab_sf_tax$newOTUid_wgen
  
   
  #Hierarchical clustering of samples
  
  otustand <- decostand(t(otutab_sf_tax_num), method = "hellinger")
  otudist <- vegdist(otustand, distance = "jaccard")
  otuclust0 <- hclust(otudist, "average")
  
  library(cluster)
  #testpvclust <- pvclust(t(otustand), method.dist = "euc", method.hclust = "average", parallel = TRUE)
  plot(testpvclust)
  
  #k     <- 4
  clust <- cutree(otuclust0,k=numclus)  # k clusters
  
  
  dendr    <- dendro_data(otuclust0, type="rectangle") # convert for ggplot
  clust.df <- data.frame(label=rownames(t(otutab_sf_tax_num)), cluster=factor(clust))
  dendr[["labels"]]   <- merge(dendr[["labels"]],clust.df, by="label")
  rect <- aggregate(x~cluster,label(dendr),range)
  rect <- data.frame(rect$cluster,rect$x)
  ymax <- mean(otuclust0$height[length(otuclust0$height)-((numclus-2):(numclus-1))])
  
  dendroplot <- ggplot() +
    #theme(plot.margin = unit(c(1,1,1,1), "lines"))+
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster),
              size=2.5) +
    geom_rect(data=rect, aes(xmin=X1-.3, xmax=X2+.3, ymin=0, ymax=ymax),
              color="red", fill=NA)+
    #geom_hline(yintercept=0.33, color="blue")+
    coord_flip() + scale_y_reverse(expand=c(0.4, 0)) + #adjust first number in "expand" to fit sample names!
    theme(legend.position = "none") +
    theme_dendro()+
    theme(plot.margin = margin(5.5, 0, 5.5, -15, "pt"))
  dendroplot
  
  # gt <- ggplot_gtable(ggplot_build(dendroplot))
  # gt$layout$clip[gt$layout$name == "panel"] <- "off"
  # grid.draw(gt) not working
  
  #https://stackoverflow.com/questions/24140339/tree-cut-and-rectangles-around-clusters-for-a-horizontal-dendrogram-in-r
  
  
  
  
  grp <- cutree(otuclust0, k =numclus)
  
  
  #####################barplot for cluster
    
    otutab_tax <- otutab_sfsep %>% filter(Divisionlong %in% taxo_group)
  
  
  otutab_tax_sf0 <- otutab_tax[,c(as.character(samplessf))]
  
  
  
  
    otutab_tax_sf1 <- sweep(otutab_tax_sf0, 2 , colSums(otutab_tax_sf0), FUN = "/")
    otutab_tax_sf1[is.na(otutab_tax_sf1)] <- 0
  
  
  otutab_tax_sf <- cbind.data.frame(otutab_tax_sf1, otutab_tax %>% select_if(negate(is.numeric)))
  otutab_tax_sf[is.na(otutab_tax_sf)] <- 0
  sjekk <- names(otutab_tax_sf)
  
  otutab_tax_sf$total <- otutab_tax_sf %>% select_if(is.numeric) %>% rowSums()
  taxlevtab <- otutab_tax_sf %>% select(taxlevel) %>% droplevels.data.frame()
  
  #####
  #Presence-absence for OTU count per taxonomic group
  otutab_sfsep_pa <- otutab_sfsep
  otutab_sfsep_pa[otutab_sfsep_pa>0] <- 1
  otutab_tax_sf_pa <- otutab_tax_sf
  otutab_tax_sf_pa[otutab_tax_sf_pa>0] <- 1
  
  taxgroupspre <-  select(otutab_tax_sf, -total) %>%   group_by_(taxlevel) %>% summarise_if(is.numeric, sum)
  
  taxgroupspre_pa <- select(otutab_tax_sf_pa, -total) %>%   group_by_(taxlevel) %>% summarise_if(is.numeric, sum)
  
  
  
  
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
  
  #taxgroups_pa_select <- taxgroupspre_pa %>% filter(get(input$taxlevel_clustplot) %in% taxgroupspre_bin_sums_yes$taxgroups)
  
#  taxgroups_pa_other <- taxgroupspre_pa %>% filter(!get(input$taxlevel_clustplot) %in% taxgroupspre_bin_sums_yes$taxgroups)
  
  
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
  
  
  if (sf == "sf10.50") {
    taxgroups_select3tdf0 <- taxgroups_select3tdf0 %>% 
    mutate(Thoracosphaeraceae = Thoracosphaeraceae_uncl+Thoracosphaeraceae_X) %>% 
    mutate(Chytriodiniaceae = Chytriodiniaceae_uncl + Chytriodinium) %>% 
      select(-c(Thoracosphaeraceae_uncl, Thoracosphaeraceae_X, Chytriodiniaceae_uncl, Chytriodinium)) #%>% 
      #select(c(1:10, 12,13,11))
  } else {
    if (sf == "sf50.200") {
      taxgroups_select3tdf0 <- taxgroups_select3tdf0 %>% 
        rename(Thoracosphaeraceae = Thoracosphaeraceae_X)
    } else {
    taxgroups_select3tdf0 <- taxgroups_select3tdf0
  }
  }
  
  taxgroups_select3tdfmelt <- melt(taxgroups_select3tdf0)
  head(taxgroups_select3tdfmelt)
  levels(taxgroups_select3tdfmelt$variable)
  
  ncolrs <- length(levels(taxgroups_select3tdfmelt$variable))-1
  #newtol <- c(tol(),"blue", "cyan")
  #barpcol <- c(alphabet2()[c(1:ncolrs)], "grey")
  #names(barpcol) <- levels(taxgroups_select3tdfmelt$variable)
  
  barplot_clust <- ggplot(taxgroups_select3tdfmelt, aes(x=sample, y = value, fill = variable))+
    geom_bar(stat = "identity", position = "stack")+
    scale_fill_manual(values = currentvec)+
    theme(axis.text.x = element_text(size =10 ), 
          axis.title.x = element_text(size =10))+
    theme(axis.text.y = element_blank())+
    theme(axis.title.y = element_blank())+
    scale_y_continuous(labels = scales::percent)+
    labs(y = currentylab, fill = "Clade")+
    theme(legend.position = "bottom")+
    guides(shape = guide_legend(override.aes = list(size = 0.7)),
           color = guide_legend(override.aes = list(size = 0.7)))+
    theme(legend.title = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.3, "lines"))+
    theme(plot.margin = margin(5.5, 1, 5.5, 1, "pt"))+
    coord_flip()
  barplot_clust
  
  
  
  #####################barplot for cluster end
  
  ######################rankotus for cluster start
  #####
  #####
  
  
  otutab_raclust_num0 <- otutab_tax_sf %>% select(-total) %>% select_if(is.numeric)
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
  
  notus <- 100 #number of OTUs to include in the "rank-abundance barplot"
  
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
    ifelse(rankphyls100[i] %in% unique(c(rankphylum[c((notus-1):notus),])), rankphyls100.2[i] <- rankphyls100[i], rankphyls100.2[i] <- "Other")
  }
  
  unique(rankphyls100.2)
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
  if (sf == "sf0.4.3") {
    seqra <- c(rep(c(1:notus),(nmp1+nmp2+nmp3+nmp4+nmp5)))
  } else { if (sf == "sf3.180") {
    seqra <- c(rep(c(1:notus),(nmp1+nmp2)))
  } else {
    seqra <- c(rep(c(1:notus),(nmp3+nmp4+nmp5)))
  }}
  
  
  dfra <- cbind.data.frame(rankvals100, rankphyls100, rankphyls100.2, samplesvecord, rankotus_notus_vec, seqra)
  
  rankfyls <- unique(rankphyls100.2) %>% str_replace("Other", replacement = NA_character_) %>% na.omit()
  rankfyls2 <- factor(c(rankfyls, "Other"), ordered = T, levels = c(rankfyls, "Other"))
  notucols <- length(unique(rankphyls100.2))-1
  otus <- tibble(rankfyls = rankfyls2) %>% mutate(color = c(alphabet(notucols),"grey"))
  otusvec <- as.character(otus$color)
  names(otusvec) <- otus$rankfyls
  
  barplot_clust_ra <- ggplot(dfra, aes(x=samplesvecord, y = rankvals100, fill = rankphyls100.2, group=seqra))+
    geom_bar(stat = "identity", position = "stack", color = "white", size = 0.1)+
    scale_fill_manual(values = currentvec_ra)+
    theme(axis.text.y = element_blank())+
    theme(axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10))+
    theme(axis.title.y = element_blank())+
    scale_y_continuous(labels = scales::percent)+
    theme(legend.text=element_text(size=8))+
    theme(legend.title = element_blank())+
    theme(legend.position = "bottom")+
    theme(legend.key.size = unit(0.5, "cm"))+
    labs(y = currentylab, fill = "Clade")+
    theme(plot.margin = margin(5.5, 1, 5.5, 1, "pt"))+
    coord_flip()
  barplot_clust_ra
  
  ####Richness####
  otutab_sf_tax_num_int <- round(otutab_tax_sf1*40000)
  rich <- specnumber(t(otutab_sf_tax_num_int))
  richdf <- tibble(sample = as.factor(names(rich)), numotus = unname(rich))
  rich_ord <- tibble(sample = unique(samplesvecord))
  richdf <- richdf %>% right_join(rich_ord, richdf, by = "sample") 
  
  richdf <- richdf %>% mutate(sample = factor(sample, ordered = T, levels = unique(sample)))
  
  richplot_clust_ra <- ggplot(data = richdf, aes(x = sample, y = numotus))+
    geom_bar(stat = "identity")+
    theme(axis.text.y = element_blank())+
    theme(axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10))+
    theme(axis.title.y = element_blank())+
    labs(y = "# OTUs")+
    theme(plot.margin = margin(5.5, 1, 5.5, 1, "pt"))+
    coord_flip()
  richplot_clust_ra
  
  ####slope raref curve###
  rarslopes <- c(NULL)
  for (i in 1:dim(otutab_sf_tax_num_int)[2]) {
    rarslopes[i] <- rareslope(t(otutab_sf_tax_num_int[,i]), sample=sum(otutab_sf_tax_num_int[,i])-1)  
  }
  names(rarslopes) <- names(otutab_sf_tax_num)
  
  slopedf <- tibble(sample = as.factor(names(rarslopes)), slope = unname(rarslopes))
  slope_ord <- tibble(sample = unique(samplesvecord))
  slopedf <- slopedf %>% right_join(slope_ord, slopedf, by = "sample") 
  
  slopedf <- slopedf %>% mutate(sample = factor(sample, ordered = T, levels = unique(sample))) %>% mutate(dummy = 1)
  
  slopeplot_clust_ra <- ggplot(data = slopedf, aes(x = sample, y = dummy))+
    geom_tile(aes(fill = slope), colour="white") + 
    geom_text(aes(label = round(slope, 4)), colour = "white", size = 4)+
    scale_fill_gradient(low="lightgrey", high= "black", na.value="white")+
    theme(axis.text.y = element_blank())+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())+
    theme(axis.title.y = element_blank())+
    theme(legend.position = "none")+
    theme(plot.margin = margin(5.5, 1, 5.5, 1, "pt"))+
    coord_flip()
  slopeplot_clust_ra
  
  
  # dendroplot_syndi_043 <- dendroplot
  # barplot_clust_syndi_043 <- barplot_clust
  # barplot_clust_ra_syndi_043 <- barplot_clust_ra
  # legend_barplot_clust_syndi_043 <- get_legend(barplot_clust_syndi_043)

  # dendroplot_syndi_310 <- dendroplot
  # barplot_clust_syndi_310 <- barplot_clust
  # barplot_clust_ra_syndi_310 <- barplot_clust_ra
  # legend_barplot_clust_syndi_310 <- get_legend(barplot_clust_syndi_043)
  # #

  # dendroplot_diatom_1050 <- dendroplot
  # barplot_clust_diatom_1050 <- barplot_clust
  # barplot_clust_ra_diatom_1050 <- barplot_clust_ra
  # # # legend_barplot_clust_diatom_1050 <- get_legend(barplot_clust_diatom_1050)
  # richplot_clust_ra_diatom_1050 <- richplot_clust_ra
  # slopeplot_clust_ra_diatom_1050 <- slopeplot_clust_ra

  dendroplot_diatom_50200 <- dendroplot
  barplot_clust_diatom_50200 <- barplot_clust
  barplot_clust_ra_diatom_50200 <- barplot_clust_ra
  #legend_barplot_clust_diatom_50200 <- get_legend(barplot_clust_diatom_1050)
  richplot_clust_ra_diatom_50200 <- richplot_clust_ra
  slopeplot_clust_ra_diatom_50200 <- slopeplot_clust_ra

  # dendroplot_dino_310 <- dendroplot
  # barplot_clust_dino_310 <- barplot_clust
  # barplot_clust_ra_dino_310 <- barplot_clust_ra
  # legend_barplot_clust_dino_310 <- get_legend(barplot_clust_syndi_043)
  # #

  # dendroplot_dinocil_1050 <- dendroplot
  # barplot_clust_dinocil_1050 <- barplot_clust
  # barplot_clust_ra_dinocil_1050 <- barplot_clust_ra
  # #legend_barplot_clust_dinocil_1050 <- get_legend(barplot_clust_dino_1050)
  # richplot_clust_ra_dinocil_1050 <- richplot_clust_ra
  # slopeplot_clust_ra_dinocil_1050 <- slopeplot_clust_ra

  # dendroplot_dinocil_50200 <- dendroplot
  # barplot_clust_dinocil_50200 <- barplot_clust
  # barplot_clust_ra_dinocil_50200 <- barplot_clust_ra
  # #legend_barplot_clust_dinocil_50200 <- get_legend(barplot_clust_dinocil_50200)
  # richplot_clust_ra_dinocil_50200 <- richplot_clust_ra
  # slopeplot_clust_ra_dinocil_50200 <- slopeplot_clust_ra
  # #
  
  testplot <- ggdraw()+
    draw_plot(dendroplot, x = 0, y = 0.5, width = 0.32, height = 0.5)+
    draw_plot(barplot_clust+theme(legend.position = "none"), x = 0.35, y = 0.5, width = 0.25, height = 0.485)+
    draw_plot(barplot_clust_ra+theme(legend.position = "none"), x = 0.6, y = 0.5, width = 0.25, height = 0.485)+
    draw_plot(ggarrange(richplot_clust_ra, slopeplot_clust_ra, align = "h", ncol = 2, widths = c(3, 1)),
              x = 0.85, y = 0.5, width = 0.15, height = 0.485)+
    draw_grob(dino_legnd, x = 0.5, y = -0.1)
  
  testplot
    
  ggdraw() +
    draw_plot(dendroplot_syndi_043, x = 0, y = 0.5, width = 0.5, height = 0.5)+
    draw_plot(barplot_clust_syndi_043+theme(legend.position = "none"), x = 0.5, y = 0.5, width = 0.25, height = 0.49)+
    draw_plot(barplot_clust_ra_syndi_043, x = 0.75, y = 0.5, width = 0.25, height = 0.49)+
    draw_plot(dendroplot_syndi_310, x = 0, y = 0, width = 0.5, height = 0.5)+
    draw_plot(barplot_clust_syndi_310+theme(legend.position = "none"), x = 0.5, y = 0, width = 0.25, height = 0.49)+
    draw_plot(barplot_clust_ra_syndi_310, x = 0.75, y = 0, width = 0.25, height = 0.49)

  ggdraw() +
    draw_plot(dendroplot_diatom_1050, x = 0, y = 0.5, width = 0.5, height = 0.5)+
    draw_plot(barplot_clust_diatom_1050+theme(legend.position = "none"), x = 0.5, y = 0.5, width = 0.25, height = 0.49)+
    draw_plot(barplot_clust_ra_diatom_1050, x = 0.75, y = 0.5, width = 0.25, height = 0.49)+
    draw_plot(dendroplot_diatom_50200, x = 0, y = 0, width = 0.5, height = 0.5)+
    draw_plot(barplot_clust_diatom_50200+theme(legend.position = "none"), x = 0.5, y = 0, width = 0.25, height = 0.49)+
    draw_plot(barplot_clust_ra_diatom_50200, x = 0.75, y = 0, width = 0.25, height = 0.49)
  
  dinocilplot_10 <- ggdraw()+
    draw_plot(dendroplot_dinocil_1050, x = 0, y = 0.5, width = 0.32, height = 0.5)+
    draw_plot(barplot_clust_dinocil_1050+theme(legend.position = "none")+ theme(axis.title.x = element_blank()), x = 0.35, y = 0.53, width = 0.25, height = 0.455)+
    draw_plot(barplot_clust_ra_dinocil_1050+theme(legend.position = "none")+ theme(axis.title.x = element_blank()), x = 0.6, y = 0.53, width = 0.25, height = 0.455)+
    draw_plot(ggarrange(richplot_clust_ra_dinocil_1050+ theme(axis.title.x = element_blank()), slopeplot_clust_ra_dinocil_1050, align = "h", ncol = 2, widths = c(3, 1)),
              x = 0.85, y = 0.53, width = 0.15, height = 0.455)+
    draw_plot(dendroplot_dinocil_50200, x = 0, y = 0.07, width = 0.32, height = 0.5)+
    draw_plot(barplot_clust_dinocil_50200+theme(legend.position = "none"), x = 0.35, y = 0.07, width = 0.25, height = 0.485)+
    draw_plot(barplot_clust_ra_dinocil_50200+theme(legend.position = "none"), x = 0.6, y = 0.07, width = 0.25, height = 0.485)+
    draw_plot(ggarrange(richplot_clust_ra_dinocil_50200, slopeplot_clust_ra_dinocil_50200, align = "h", ncol = 2, widths = c(3, 1)),
              x = 0.85, y = 0.07, width = 0.15, height = 0.485)+
    draw_grob(dino_legnd_10200, x = 0.01, y = -0.45)+
    draw_grob(dino_legnd_10200_ra, x = 0.7, y = .02, width = 0.01, height = 0.05, scale = 0.05)+
    draw_grob(text_grob("Order"), x = -0.12, y = 0.48)+
    draw_grob(text_grob("Rank-abundance OTUs"), x = 0.2, y = 0.48)+
    draw_grob(text_grob("Slope"), x = 0.48, y = 0.48)+
    draw_grob(text_grob("10-50", rot = 90), x = -0.48, y = 0.28)+
    draw_grob(text_grob("50-200", rot = 90), x = -0.48, y = -0.20)
  dinocilplot_10
  
####diatoms 10-200
  diatomplot_10 <- ggdraw()+
    draw_plot(dendroplot_diatom_1050, x = 0, y = 0.5, width = 0.32, height = 0.5)+
    draw_plot(barplot_clust_diatom_1050+theme(legend.position = "none")+ theme(axis.title.x = element_blank()), x = 0.35, y = 0.53, width = 0.25, height = 0.455)+
    draw_plot(barplot_clust_ra_diatom_1050+theme(legend.position = "none")+ theme(axis.title.x = element_blank()), x = 0.6, y = 0.53, width = 0.25, height = 0.455)+
    draw_plot(ggarrange(richplot_clust_ra_diatom_1050+ theme(axis.title.x = element_blank()), slopeplot_clust_ra_diatom_1050, align = "h", ncol = 2, widths = c(3, 1)),
              x = 0.85, y = 0.53, width = 0.15, height = 0.455)+
    draw_plot(dendroplot_diatom_50200, x = 0, y = 0.07, width = 0.32, height = 0.5)+
    draw_plot(barplot_clust_diatom_50200+theme(legend.position = "none"), x = 0.35, y = 0.07, width = 0.25, height = 0.485)+
    draw_plot(barplot_clust_ra_diatom_50200+theme(legend.position = "none"), x = 0.6, y = 0.07, width = 0.25, height = 0.485)+
    draw_plot(ggarrange(richplot_clust_ra_diatom_50200, slopeplot_clust_ra_diatom_50200, align = "h", ncol = 2, widths = c(3, 1)),
              x = 0.85, y = 0.07, width = 0.15, height = 0.485)+
    draw_grob(bacill_legnd_10200, x = 0.01, y = -0.45)+
    draw_grob(bacill_legnd_10200_ra, x = 0.72, y = .02, width = 0.01, height = 0.05, scale = 0.05)+
    draw_grob(text_grob("Genus"), x = -0.12, y = 0.48)+
    draw_grob(text_grob("Rank-abundance OTUs"), x = 0.2, y = 0.48)+
    draw_grob(text_grob("Slope"), x = 0.48, y = 0.48)+
    draw_grob(text_grob("10-50", rot = 90), x = -0.48, y = 0.28)+
    draw_grob(text_grob("50-200", rot = 90), x = -0.48, y = -0.20)
  diatomplot_10
  ###
  
  
  ###Test varpart
  otutab_sf_tax_num_hel <- decostand(t(otutab_sf_tax_num), "hellinger")
  vp.dino.310 <- varpart(otutab_sf_tax_num_hel, meta_sf$Cruise, as.numeric(as.character(meta_sf$DEPTH_M)))
  plot(vp.dino.310)
  vp.dino.310
  str(vp.dino.310)
  vp.dino.310$part$indfract$Adj.R.squared[c(1,3)]
  
  #Test of fraction [a+b]
  anova1.dino.310 <- anova(rda(otutab_sf_tax_num_hel, meta_sf$Cruise))
  anova1.dino.310
  #Test of fraction [b+c]
  anova2.dino.310 <- anova(rda(otutab_sf_tax_num_hel, as.numeric(as.character(meta_sf$DEPTH_M))))
  anova2.dino.310
  #Test of fraction [a+b+c]
  c_d <- cbind.data.frame(cruise=as.factor(meta_sf$Cruise), depth=as.numeric(as.character(meta_sf$DEPTH_M)))
  anova3.dino.310 <- anova(rda(otutab_sf_tax_num_hel, c_d))
  anova3.dino.310
  #Test of fraction [a]
  anova4.dino.310 <- anova(rda(otutab_sf_tax_num_hel, meta_sf$Cruise, as.numeric(as.character(meta_sf$DEPTH_M))))
  anova4.dino.310
  #Test of fraction [c]
  anova5.dino.310 <- anova(rda(otutab_sf_tax_num_hel, as.numeric(as.character(meta_sf$DEPTH_M)), meta_sf$Cruise))
  anova5.dino.310
  
  
  ####
  otutab_sf_tax_num_hel <- decostand(t(otutab_sf_tax_num), "hellinger")
  vp.syndi.310 <- varpart(otutab_sf_tax_num_hel, meta_sf$Cruise, as.numeric(as.character(meta_sf$DEPTH_M)))
  plot(vp.syndi.310)
  vp.syndi.310
  
  #Test of fraction [a+b]
  anova1.syndi.310 <- anova(rda(otutab_sf_tax_num_hel, meta_sf$Cruise))
  anova1.syndi.310
  #Test of fraction [b+c]
  anova2.syndi.310 <- anova(rda(otutab_sf_tax_num_hel, as.numeric(as.character(meta_sf$DEPTH_M))))
  anova2.syndi.310
  #Test of fraction [a+b+c]
  c_d <- cbind.data.frame(cruise=as.factor(meta_sf$Cruise), depth=as.numeric(as.character(meta_sf$DEPTH_M)))
  anova3.syndi.310 <- anova(rda(otutab_sf_tax_num_hel, c_d))
  anova3.syndi.310
  #Test of fraction [a]
  anova4.syndi.310 <- anova(rda(otutab_sf_tax_num_hel, meta_sf$Cruise, as.numeric(as.character(meta_sf$DEPTH_M))))
  anova4.syndi.310
  #Test of fraction [c]
  anova5.syndi.310 <- anova(rda(otutab_sf_tax_num_hel, as.numeric(as.character(meta_sf$DEPTH_M)), meta_sf$Cruise))
  anova5.syndi.310
  
  
  ####Rarefaction
  1/min(otutab_tax_sf1[otutab_tax_sf1>0])
  otutab_sf_tax_num_int <- round(otutab_tax_sf1*40000)
  min(otutab_sf_tax_num_int[otutab_sf_tax_num_int>0])
 rarcurv <- rarecurve(t(round(otutab_sf_tax_num_int)), sample = 40000, step = 1000)

 plot(attr(rarcurv[[1]], which = "Subsample"), rarcurv[[1]], xlim=c(0,30000), ylim=c(0,700), type="l")
 for (i in 1:dim(otutab_sf_tax_num)[2]) {
   lines(attr(rarcurv[[i]], which = "Subsample"), rarcurv[[i]])
 }
 
 rarslopes <- c(NULL)
 for (i in 1:dim(otutab_sf_tax_num)[2]) {
   rarslopes[i] <- rareslope(t(otutab_sf_tax_num_int[,i]), sample=sum(otutab_sf_tax_num_int[,i])-1)  
 }
 names(rarslopes) <- names(otutab_sf_tax_num)
 
 ChaoSpecies(t(otutab_sf_tax_num_int), datatype = "abundance")
 
 Diversity(t(otutab_sf_tax_num_int)[19,], datatype = "abundance")
 diversity(t(otutab_sf_tax_num_int))
 specnumber(t(otutab_sf_tax_num_int))
 plot(diversity(t(otutab_sf_tax_num_int)), specnumber(t(otutab_sf_tax_num_int)))
 