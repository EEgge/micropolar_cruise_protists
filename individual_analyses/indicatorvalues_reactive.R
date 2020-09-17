###Indicator values   
indval <- reactive({
  sfsamp <- sfnames[[input$sizefract_iva]]
  #zone <- input$zone_iva
  tax <- input$taxo_group_iva
  numclust <- input$numclust_iva
  
  sjekk1 <- "test"
  sjekk2 <- input$taxo_group_iva
  # sjekk2 <- get(input$zone_iva)
  otutab_sf <- otutab_sfsep[c(as.character(sfsamp), "Divisionlong", "newOTUid_wgen")]
  
  otutab_tax <- otutab_sf %>% filter(Divisionlong %in% tax) %>% droplevels()
  otutab_tax_num <- otutab_tax %>% select_if(is.numeric)
  otutab_tax_notnum <- otutab_tax %>% select_if(Negate(is.numeric))
  otutab_tax_prop0 <- sweep(otutab_tax_num, 2, colSums(otutab_tax_num), FUN = "/")
  otutab_tax_prop <- cbind(otutab_tax_prop0, otutab_tax_notnum)
  otutab_tax_prop[is.na(otutab_tax_prop)] <- 0 # ok to replace NAs with 0, because when the taxo group is 0, the OTU is also 0
  otutab_tax <- otutab_tax_prop
  
  otustand <- decostand(t(otutab_tax_num), method = "hellinger")
  braydist <- vegdist(otustand, distance = "jaccard")
  
  otuclust0 <- hclust(braydist, "average")
  plot(otuclust0)
  grp <- cutree(otuclust0, k =input$clusters)
  # 
  # 
  otutab_sf_list <- list()
  for (i in 1:length(sfnames)) {
    otutab_sf_list[[i]] <- otutab_tax[c(as.character(sfnames[[i]]), "Divisionlong", "newOTUid_wgen")]
  }
  # 
  # dim(otutab_sf_list[[1]])
  # 
  sampnames_meta_sf_list <- list()
  for (i in 1:length(sfnames)) {
    sampnames_meta_sf_list[[i]] <- meta_table[names(sfnames)[i]]
  }
  # 
  meta_sf_list <- list()
  for (i in 1:length(sfnames)) {
    meta_sf_list[[i]] <- meta_table[match(colnames(otutab_sf_list[[i]] %>% select(-Divisionlong, -newOTUid_wgen)),as.character(sampnames_meta_sf_list[[i]][,1])),] %>% droplevels()
  }
  # 
  # 
  otutab_tax_num_list <- list()
  for (i in 1:length(sfnames)) {
    otutab_tax_num_list[[i]] <- otutab_sf_list[[i]] %>% select_if(is.numeric) %>% t() %>% as.data.frame()
    
  }
  
  for (i in 1:length(sfnames)) {
    names(otutab_tax_num_list[[i]]) <- otutab_sf_list[[i]]$newOTUid_wgen
  }
  
  for (i in 1:length(sfnames)) {
    otutab_tax_num_list[[i]]$sample <- rownames(otutab_tax_num_list[[i]])
  }
  
  zone_samp_list <- list()
  for (i in 1:length(sfnames)) {
    zone_samp_list[[i]] <- meta_sf_list[[i]] %>% filter(depthbin %in% zone) %>% select(sfs[i])
  }
  
  
  # 
  # 
  # #get cruise names
  zone_cruise_list <- list()
  for (i in 1:length(sfnames)) {
    zone_cruise_list[[i]] <- meta_sf_list[[i]] %>% filter(depthbin %in% zone) %>% select(Cruise)
  }
  
  
  # 
  otutab_tax_num_zone_list <- list()
  for (i in 1:length(sfnames)) {
    otutab_tax_num_zone_list[[i]] <- otutab_tax_num_list[[i]] %>% filter(sample %in% as.character(zone_samp_list[[i]][,1])) %>% select_if(is.numeric)
  }
  # 
  otutab_tax_num_zone_pa_list <- list()
  for (i in 1:length(sfnames)) {
    otutab_tax_num_zone_pa_list[[i]] <- ifelse(otutab_tax_num_zone_list[[i]]>0,1,0)
  }
  
  iva.clust <- multipatt(otutab_tax_prop0, grp, func = "r.g", max.order = 4, control = how(nperm = 999))
  
  iva.tax.zone.list <- list()
  for (i in 1:length(sfnames)) {
    iva.tax.zone.list[[i]] <- multipatt(otutab_tax_num_zone_list[[i]], zone_cruise_list[[i]]$Cruise,
                                        func = "r.g",
                                        max.order = 5,
                                        control = how(nperm = 999))
  }
  # 
  ptab.tax.zone.list <- list()
  for (i in 1:length(sfnames)) {
    ptab.tax.zone.list[[i]] <- data.frame(otu = rownames(iva.tax.zone.list[[i]]$str), pval = iva.tax.zone.list[[i]]$sign$p.value) %>% filter(pval<= 0.05)
  }
  
  otutab_sfsep_tax <- otutab_sfsep %>% filter(Divisionlong == input$taxo_group_iva) %>% droplevels.data.frame()
  
  otusign.mtax.list <- list()
  for (i in 1:length(sfnames)) {
    otusign.mtax.list[[i]] <- otutab_sfsep_tax %>% filter(newOTUid_wgen %in% ptab.tax.zone.list[[i]]$otu) %>% select_if(Negate(is.numeric))
  }
  
  
  taxlev_sign_count_list <- list()
  for (i in 1:length(sfnames)) {
    taxlev_sign_count_list[[i]] <- as.data.frame(table(otusign.mtax.list[[i]]$Family))
  }
  
  
  otutab_tax_zone_list <- list()
  for (i in 1:length(sfnames)) {
    otutab_tax_zone_list[[i]] <- otutab_tax %>% select(c(as.character(zone_samp_list[[i]][,1]), "Divisionlong", "Class", "Order", "Family", "Genus", "newOTUid_wgen"))
  }
  
  sjekk1 <- dim(otutab_tax_zone_list[[1]])[2]
  # 
  otutab_tax_zone_pres_list <- list()
  for (i in 1:length(sfnames)) {
    otutab_tax_zone_pres_list[[i]] <- otutab_tax_zone_list[[i]] %>% mutate(total = rowSums(select_if(otutab_tax_zone_list[[i]], is.numeric))) %>% filter(total > 0)
  }
  
  taxlev_total_count_list <- list()
  for (i in 1:length(sfnames)) {
    taxlev_total_count_list[[i]] <- as.data.frame(table(otutab_tax_zone_pres_list[[i]]$Family))
  }
  
  indicspec_table_list <- list()
  for (i in 1:length(sfnames)) {
    indicspec_table_list[[i]] <- data.frame("fam" = taxlev_total_count_list[[i]]$Var1, "total" = taxlev_total_count_list[[i]]$Freq, "sign" = taxlev_sign_count_list[[i]]$Freq) %>%
      mutate("perc" = (sign/total)*100)
  }
  
  # 
  indicplot_tax_list <- list()
  for (i in 1:length(sfnames)) {
    indicplot_tax_list[[i]] <- ggplotly(ggplot(indicspec_table_list[[i]], aes(x = fam, y = perc, text = sprintf("Family: %s<br>Total: %s<br>Sign.: %s ", fam, total, sign)))+
                                          geom_bar(stat = "identity")+
                                          geom_text(aes(x=fam, y = perc, label =round(perc, 1)))+
                                          theme_minimal()+
                                          theme(axis.text.x = element_text(angle = 90))+
                                          ylab("Percent seasonal OTUs"), tooltip = "text")
  }
  
  mulpat_tab_list <- list()
  for (i in 1:length(sfnames)) {
    mulpat_tab_list[[i]] <- iva.tax.zone.list[[i]]$sign %>% mutate(otuid = rownames(.)) %>% filter(p.value <= 0.05)
  }
  
  for (i in 1:length(sfnames)) {
    mulpat_tab_list[[i]]$Family <- otutab_sfsep %>% filter(newOTUid_wgen %in% mulpat_tab_list[[i]]$otuid) %>% pull(Family)
  }
  
  mulpat_tab_all_list <- list()
  for (i in 1:length(sfnames)) {
    mulpat_tab_all_list[[i]] <- iva.tax.zone.list[[i]]$sign %>% mutate(otuid = rownames(.))
  }
  
  for (i in 1:length(sfnames)) {
    mulpat_tab_all_list[[i]]$Family <- otutab_sfsep %>% filter(newOTUid_wgen %in% mulpat_tab_all_list[[i]]$otuid) %>% pull(Family)
  }
  
  mulpat_forbp_list <- list()
  for (i in 1:length(sfnames)) {
    mulpat_forbp_list[[i]] <- mulpat_tab_list[[i]] %>% group_by(Family) %>% count(index)
  }
  
  season_mulpatplot_list <- list()
  for (i in 1:length(sfnames)) {
    season_mulpatplot_list[[i]] <- ggplotly(ggplot(mulpat_forbp_list[[i]], aes(x=index, y=n, fill = Family, text = sprintf("Taxon: %s<br>n OTUs: %s", Family, n)))+
                                              geom_bar(stat = "identity", position = "stack")  )
  }
  
  
  
  for (i in 1:length(sfnames)) {
    k <- which(names(mulpat_tab_all_list[[i]]) == "otuid")
    names(mulpat_tab_all_list[[i]])[k] <- "newOTUid_wgen"
  }
  
  
  
  otutab_tax <- otutab_sfsep %>% filter(Divisionlong %in% tax) %>% droplevels()
  
  tax_for_multipatt <- otutab_tax
  tax_for_multipatt$total = tax_for_multipatt %>% select_if(is.numeric) %>% rowSums()
  # 
  tax_for_multipatt <- tax_for_multipatt %>% arrange(desc(total))
  
  names(mulpat_tab_list[[1]])
  
  tax_multipatt_list <- list()
  for (i in 1:length(sfnames)) {
    tax_multipatt_list[[i]] <- left_join(tax_for_multipatt, mulpat_tab_all_list[[i]], by = "newOTUid_wgen")
  }
  
  indextab043 <- data.frame(ses = c("Jan", "Mar", "May", "Aug", "Nov", "Jan+Mar", 
                                    "Jan+May", "Jan+Aug", "Jan+Nov", "Mar+May", "Mar+Aug",
                                    "Mar+Nov", "May+Aug", "May+Nov", "Aug+Nov", 
                                    "Jan+Mar+May", "Jan+Mar+Aug", "Jan+Mar+Nov",
                                    "Jan+May+Aug", "Jan+May+Nov", "Jan+Aug+Nov",
                                    "Mar+May+Aug", "Mar+May+Nov", "Mar+Aug+Nov",
                                    "May+Aug+Nov", "Jan+Mar+May+Aug", "Jan+May+Aug+Nov",
                                    "Jan+Mar+Aug+Nov", "Jan+Mar+May+Nov", "Mar+May+Aug+Nov"),
                            index = c(1:30))
  
  indextab3200 <- data.frame(ses = c("May", "Aug", "Nov", "May+Aug", "May+Nov", "Aug+Nov", "May+Aug+Nov"),
                             index = c(1:7))
  
  indextab3180 <- data.frame(ses = c("Jan", "Mar", "Jan+Mar"),
                             index = c(1:3))
  
  tax_multipatt_list_final <- list()
  
  tax_multipatt_list_final[[1]] <- left_join(tax_multipatt_list[[1]][,c(142:161)], indextab043, by = "index")
  
  for (i in 2:4) {
    tax_multipatt_list_final[[i]] <- left_join(tax_multipatt_list[[i]][,c(140:159)], indextab3200, by = "index")
  }
  
  tax_multipatt_list_final[[5]] <- left_join(tax_multipatt_list[[5]][,c(139:158)], indextab3180, by = "index")
  
  list(iva.tax.zone.list = iva.tax.zone.list, indicplot_tax_list = indicplot_tax_list, season_mulpatplot_list = season_mulpatplot_list, sjekk1 = sjekk1, sjekk2 = sjekk2,
       tax_multipatt_list_final =tax_multipatt_list_final)
})
