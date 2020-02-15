#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#detach("package:here", unload=TRUE)
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

bpps <- 0.6
otutab_sfsep0 <- read.table(here("data", "OTUtab_nozooembr_minsize5_prop_wtax_wnewseqid_20191022_point.txt"), header = T, sep = "\t")[,-c(141:155)]
#otutab_sfsep0 <- read.table(here("data", "OTUtab_nozooembr_minsize5_prop_wtax_wnewseqid_14032019_point2.txt"), header = T, sep = "\t")[,-c(141:155)]
# syndi <- otutab_sfsep %>% filter(Divisionlong %in% "Dinoflagellata_Syndiniales")
# syndi$total = syndi %>% select_if(is.numeric) %>% rowSums()
# 
# syndi <- syndi %>% arrange(desc(total))
# view(head(syndi, n= 50))
# write.table(syndi[,c(143:152)], "syndi.txt")

dim(otutab_sfsep0)
otutab_sfsep <- otutab_sfsep0[,c(5:8,1:4,9:122,127:130,123:126,131:151)]
names(otutab_sfsep)
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



# Define UI for application that analyses taxonomic composition
#options(error= recover)
ui <- navbarPage(
   
   # Application title
  "Taxonomic composition of MicroPolar protist samples",
  
  #Tab showing barplot of proportional read abundance and proportional OTU richness
  tabPanel("Proportional abundance",
           fluidPage(sidebarLayout(
             sidebarPanel(
               checkboxGroupInput("taxo_group2", label = "Division", 
                                  choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                              "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "Pseudofungi",  "Radiolaria",  "All"), 
                                  selected = "All"),
               radioButtons("taxlevel", label = "Taxonomic level",
                            choices = c("Class", "Order", "Family", "Genus", "Species")),
               conditionalPanel(condition = "input.taxo_group2 == 'Dinoflagellata_Syndiniales'",
                                radioButtons("family", label = "Family:", choices = c("Clade 1", "Clade 2")))
               
             ),
          
             mainPanel(plotlyOutput(outputId = "propPlot1", height = 400, width = 850),
                       plotlyOutput(outputId = "propPlot2", height = 600, width = 1500),
                      plotlyOutput(outputId = "propPlotTot",height = 200, width = 1500),
                       plotlyOutput(outputId = "propPlot_pa1", height = 400, width = 850),
                       plotlyOutput(outputId = "propPlot_pa2", height = 600, width = 1500),
                      plotlyOutput(outputId = "propPlotTot_pa",height = 200, width = 1500)
                       )))),
  
  # Tab with NMDS plot, hclust dendrogram and table with p-values for test of differential OTU composition between depths
  tabPanel("Beta diversity - NMDS",
   
   # Sidebar with radio buttons to select size fraction, checkbox to select one or more taxonomic groups
   sidebarLayout(
     sidebarPanel(
       radioButtons("sizefract1", label = "Size fraction", 
                    choices = c("sf0.4.3", "sf3.180",  "sf3.10", "sf10.50", "sf50.200"), selected = NULL),
       checkboxGroupInput("taxo_group", label = "Division", 
                    choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "All"), 
                    selected = "All"),
       radioButtons("manyglm", label = "manyglm test of different OTU composition between depth zones (may take some time to compute)", choices = c("yes", "no"), selected = "no"),
       #Buttons to download results
       downloadButton("downloadData_glmp3", "Download p-values manyglm May"),
       downloadButton("downloadData_glmp4", "Download p-values manyglm Aug"),
       downloadButton("downloadData_glmp5", "Download p-values manyglm Nov"),
       downloadButton("downloadOTUtab", "Download OTU table"),
       downloadButton("downloadnmdstab", "Download nmds table")
     ),
      
       
     mainPanel(
       plotlyOutput(outputId = "distPlot"),
       textOutput(outputId = "stress"),
       plotOutput(outputId = "plotclust"),
       tableOutput(outputId = "difftab"),
       textOutput(outputId = "numotus"),
       textOutput(outputId = "anovaja"),
       tableOutput(outputId = "sigmp1"),
       tableOutput(outputId = "sigmp2"),
       tableOutput(outputId = "sigmp3"),
       tableOutput(outputId = "sigmp4"),
       tableOutput(outputId = "sigmp5")
     ))),
  
  #Multivariate test of different OTU communities between cruises
  tabPanel("manyglm test - season",
           
           # Sidebar with radio buttons to select size fraction, checkbox to select one or more taxonomic groups
           sidebarLayout(
             sidebarPanel(
               radioButtons("sizefract3", label = "Size fraction", 
                            choices = c("sf0.4.3", "sf3.180",  "sf3.10", "sf10.50", "sf50.200"), selected = NULL),
               checkboxGroupInput("taxo_group4", label = "Division", 
                                  choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                              "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "All"), 
                                  selected = "All"),
               radioButtons("zone", label = "Depth zone", choices = c("epi", "meso"), selected = NULL),
               radioButtons("manyglm2", label = "manyglm test of different OTU composition by cruise (may take some time to compute)", choices = c("yes", "no"), selected = "no")),
               
             
             
             mainPanel(
               tableOutput(outputId = "tax_month_sf_sign_1_zone")
               ))),
  
  
  #Tab that shows rank-abundance of the OTUs within the selected Division(s). Color code shows the selected taxonomic level.
  tabPanel("Rank-abundance",
           fluidPage(sidebarLayout(
             sidebarPanel(
               checkboxGroupInput("taxo_group_ra", label = "Division", 
                                  choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                              "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "Radiolaria",  "All"), 
                                  selected = "All"),
               radioButtons("taxlevel_ra", label = "Taxonomic level",
                            choices = c("Class"="Class", "Order"="Order", "Family"="Family", "Genus"="Genus", "OTU"="newOTUid_wgen"))
               
             ),
             mainPanel(
               plotlyOutput(outputId = "raPlot1", height = 400, width = 1500),
               plotlyOutput(outputId = "raPlot2", height = 600, width = 1500)
             )))),
  
  #Tab with plot of alpha diversity (OTU richness) by cruise and size fraction for the selected division(s)
  tabPanel("Alpha diversity - richness & evenness",
           fluidPage(sidebarLayout(
             sidebarPanel(
               checkboxGroupInput("taxo_group_adiv", label = "Division", 
                                  choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                              "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "Radiolaria",  "All"), 
                                  selected = "All")),
             mainPanel(
               fluidRow(column(5, plotlyOutput(outputId = "richPlot1", height = 350)),
               column(7, plotlyOutput(outputId = "richPlot2", height = 600))),
               fluidRow(column(5, plotlyOutput(outputId = "evenPlot1", height = 350)),
               column(7, plotlyOutput(outputId = "evenPlot2", height = 600)))
             )))),
  
  #Tab with plot of so-called "Local contribution to beta diversity"
  tabPanel("Beta diversity",
           fluidPage(sidebarLayout(
             sidebarPanel(
               radioButtons("sizefract2", label = "Size fraction", 
                            choices = c("sf0.4.3", "sf3.180",  "sf3.10", "sf10.50", "sf50.200"), selected = NULL),
               
               checkboxGroupInput("taxo_group_bdiv", label = "Division", 
                                  choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                              "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "Radiolaria",  "All"), 
                                  selected = "All")
             ),
             mainPanel("Interpretation of Local contribution to beta diversity (lcbd): see Legendre & De Cáceres (2013)",
               fluidRow(column(12, plotlyOutput(outputId = "lcbdPlott", height = 350))),
               textOutput(outputId = "BDtotal" ))
             
                        
             
           ))),
  #Tab with proportional read abundance of Chloroplast 16S data (DNA and RNA)
  tabPanel("Chloroplast 16S data",
           sidebarLayout(
             sidebarPanel(
               checkboxGroupInput("taxo_group3", label = "Division", 
                                  choices = c("Chlorophyta", 
                                              "Cryptophyta", "Dinophyta", "Discoba",
                                              "Haptophyta", "Ochrophyta", "Rappemonad",  "All"), 
                                  selected = "All"),
               radioButtons("taxlevel_chl", label = "Taxonomic level",
                            choices = c("Class", "SubClass", "Order", "SubOrder",  "Family", "Genus", "Species"))
               
             ),
             mainPanel(
               plotlyOutput(outputId = "propPlotchl", height = 1000, width = 1000)
             )
           )
  )
  ,
  tabPanel("Characteristic OTUs",
            sidebarLayout(
              sidebarPanel(
                # checkboxGroupInput("month", label = "Month",
                #                    choices = c("January", "March", "May", "August", "November", "All"),
                #                    selected = "All"),
                # radioButtons("sizefract_iva", label = "Size fraction",
                #                    choices = c("sf0.4.3","sf3.180", "sf3.10", "sf10.50", "sf50.200"),
                #                    selected = NULL),
                radioButtons("taxo_group_iva", label = "Division", 
                                   choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                               "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "Radiolaria"), 
                                   selected = "Dinoflagellata_Syndiniales"),
                radioButtons("zone_iva", label = "Depth zone",
                                   choices = c("Epipelagic"="epi", "Mesopelagic"="meso"),
                                   selected = NULL)
                ),
                mainPanel(verbatimTextOutput(outputId = "iva"),
                          textOutput(outputId = "sjekk1"),
                          textOutput(outputId = "sjekk2"),
                          plotlyOutput(outputId = "iva0.4.3", height = 600, width = 1500),
                          plotlyOutput(outputId = "iva3.180", height = 600, width = 1500),
                          plotlyOutput(outputId = "iva3.10", height = 600, width = 1500),
                          plotlyOutput(outputId = "iva10.50", height = 600, width = 1500),
                          plotlyOutput(outputId = "iva50.200", height = 600, width = 1500)
                          )
                
              )
            )
  # ,
  # tabPanel("OTU counts",
  #          sidebarLayout(
  #            sidebarPanel(
  #              # checkboxGroupInput("month", label = "Month",
  #              #                    choices = c("January", "March", "May", "August", "November", "All"),
  #              #                    selected = "All"),
  #              radioButtons("sizefract_iva", label = "Size fraction",
  #                           choices = c("sf0.4.3","sf3.180", "sf3.10", "sf10.50", "sf50.200"),
  #                           selected = NULL),
  #              checkboxGroupInput("taxo_group_iva", label = "Division", 
  #                                 choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
  #                                             "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "Picozoa", "Stramenopiles_X", "Radiolaria",  "All"), 
  #                                 selected = "All"),
  #              radioButtons("zone_iva", label = "Depth zone",
  #                           choices = c("Epipelagic", "Mesopelagic"),
  #                           selected = NULL)
  #            ),
  #            mainPanel(verbatimTextOutput(outputId = "iva"),
  #                      textOutput(outputId = "sjekk1"))
  #            
  #          )
  #          )
  # 
  )
   



#Create barplot with proportional abundance for the selected Division(s), at the selected taxonomic level  
server <- function(input, output, session) {
  prop_tax <- reactive({
    taxxx <- input$taxlevel
    
    if (input$taxo_group2 != "All") {
      
      otutab_tax <- otutab_sfsep %>% filter(Divisionlong %in% input$taxo_group2)
    } else {
      otutab_tax <- otutab_sfsep
    }
    
    # otutab_tax_num <- otutab_tax %>% select_if(is.numeric) %>% select(-total)
    # otutab_tax_numt <- as.data.frame(t(otutab_tax_num))
    # names(otutab_tax_numt) <- otutab_tax$newOTUid_wgen
    # otutab_tax_numt$sample <- rownames(otutab_tax_numt)
    
    otutab_tax$total <- otutab_tax %>% select_if(is.numeric) %>% rowSums()
    taxlevtab <- otutab_tax %>% select(input$taxlevel) %>% droplevels.data.frame()
    
    #Presence-absence for OTU count per taxonomic group
    otutab_sfsep_pa <- otutab_sfsep
    otutab_sfsep_pa[otutab_sfsep_pa>0] <- 1
    otutab_tax_pa <- otutab_tax
    otutab_tax_pa[otutab_tax_pa>0] <- 1
    
    taxgroupspre <-  select(otutab_tax, -total) %>%   group_by_(input$taxlevel) %>% summarise_if(is.numeric, sum)
    
    taxgroupspre_pa <- select(otutab_tax_pa, -total) %>%   group_by_(input$taxlevel) %>% summarise_if(is.numeric, sum)
    
    
    taxgroups <- melt(taxgroupspre)
    
    
    limfun <- function(x) {
      ifelse(x>=0.05,1,0)
    }
    
    taxgroupspre_mat <- as.matrix(taxgroupspre[,-1, drop = FALSE])
    taxgroupspre_bin <- apply(taxgroupspre_mat, 2, limfun)
    
    if (is.vector(taxgroupspre_bin)) {
      taxgroupspre_bin2 <- as.data.frame(as.list(taxgroupspre_bin))
    } else {
      taxgroupspre_bin2 <- as.data.frame(taxgroupspre_bin)}
    
  
    taxgroupspre_bin_sums <- taxgroupspre_bin2 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
    taxgroupspre_bin_sums$taxgroups <- as.vector(taxgroupspre[[taxxx]])
    taxgroupspre_bin_sums_yes <- filter(taxgroupspre_bin_sums, rowsumm >0) %>% select(taxgroups)
    
    
    taxgroups_select <- taxgroupspre %>% filter(get(input$taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)
   
    taxgroups_other <- taxgroupspre %>% filter(!get(input$taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)
    
    taxgroups_pa_select <- taxgroupspre_pa %>% filter(get(input$taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)
    
    taxgroups_pa_other <- taxgroupspre_pa %>% filter(!get(input$taxlevel) %in% taxgroupspre_bin_sums_yes$taxgroups)
  
    
    if (dim(taxgroupspre_bin_sums_yes)[1] < dim(taxgroupspre_bin_sums)[1]) {
      Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups,"Other")
      othersum <- colSums(taxgroups_other[,-1])
    } else {
      Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups)
      othersum <- NULL
    }
    
    
    
    taxgroups_select2 <- rbind(taxgroups_select[,-1],othersum)
    
    
    
    taxgroups_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_select2)
    
    taxgroups_select3t <- t(taxgroups_select3)
    taxgroups_select3t2 <- taxgroups_select3t[-1,]
    taxgroups_select3t3 <- apply(taxgroups_select3t2, 2, as.numeric)
    taxgroups_select3tdf0 <- as.data.frame(taxgroups_select3t3)
    names(taxgroups_select3tdf0) <- Taxonomic_group
    
    total04.3 <- taxgroups_select3 %>% select(samples0.4.3) %>% mutate(rowsumm = rowSums(.))
    total3.10 <- taxgroups_select3 %>% select(samples3.10) %>% mutate(rowsumm = rowSums(.))
    total10.50 <- taxgroups_select3 %>% select(samples10.50) %>% mutate(rowsumm = rowSums(.))
    total50.200<- taxgroups_select3 %>% select(samples50.200) %>% mutate(rowsumm = rowSums(.))
    total3.180 <- taxgroups_select3 %>% select(samples3.180) %>% mutate(rowsumm = rowSums(.))
    
    
    totaltotal04.3 <- otutab_sfsep %>% select(samples0.4.3) %>% mutate(rowsumm = rowSums(.))
    totaltotal3.10 <- otutab_sfsep %>% select(samples3.10) %>% mutate(rowsumm = rowSums(.))
    totaltotal10.50 <- otutab_sfsep %>% select(samples10.50) %>% mutate(rowsumm = rowSums(.))
    totaltotal50.200<- otutab_sfsep %>% select(samples50.200) %>% mutate(rowsumm = rowSums(.))
    totaltotal3.180 <- otutab_sfsep %>% select(samples3.180) %>% mutate(rowsumm = rowSums(.))
    
    taxgroups_select3tdf <- rbind(taxgroups_select3tdf0, total04.3$rowsumm/sum(totaltotal04.3$rowsumm), total3.10$rowsumm/sum(totaltotal3.10$rowsumm), total10.50$rowsumm/sum(totaltotal10.50$rowsumm), total50.200$rowsumm/sum(totaltotal50.200$rowsumm), total3.180$rowsumm/sum(totaltotal3.180$rowsumm))
    
    taxgroups_select3tdf$sf <- factor(c(rep("sf0.4.3", npico), rep("sf3.10", nthree10), 
                                        rep("sf10.50", nten50), rep("sf50.200", nfifty200), rep("sf3.180", nthree180), c("sf0.4.3", "sf3.10", "sf10.50", "sf50.200", "sf3.180")), ordered = T, 
                                      levels = c("sf0.4.3", "sf3.10", "sf3.180", "sf10.50", "sf50.200"))
    
    taxgroups_select3tdf$stdep <-  factor(c(as.character(StationDepthsfsep_tb), as.character(rep("total",5))), levels = c(levels(StationDepthsfsep_tb), "total"), ordered = T)
    taxgroups_select3tdf$cruise <-  factor(c(cruisesfsep_tb$cruise, rep("total",5)), ordered = T)
    
    
    # taxgroups_select3tdf$sf <- factor(c(rep("sf0.4.3", npico), rep("sf3.10", nthree10), 
    #                      rep("sf10.50", nten50), rep("sf50.200", nfifty200), rep("sf3.180", nthree180)), ordered = T, 
    #                      levels = c("sf0.4.3", "sf3.10", "sf3.180", "sf10.50", "sf50.200"))
    # 
    # taxgroups_select3tdf$stdep <-  StationDepthsfsep_tb
    # taxgroups_select3tdf$cruise <-  factor(cruisesfsep_tb$cruise, ordered = T)
    # 
   
    
    taxgroups_select3tdf_mp12 <- taxgroups_select3tdf %>% filter(cruise %in% c("MP1", "MP2"))
    
    taxgroups_select3tdf_mp345 <- taxgroups_select3tdf %>% filter(cruise %in% c("MP3", "MP4", "MP5"))
    
    taxgroups_select3tdf_total <- taxgroups_select3tdf %>% filter(cruise %in% c("total"))
    
    
   taxgroups_select3tdf_mp12melt <- melt(taxgroups_select3tdf_mp12)
   taxgroups_select3tdf_mp345melt <- melt(taxgroups_select3tdf_mp345)
   taxgroups_select3tdf_totalmelt <- melt(taxgroups_select3tdf_total)
    
   p12 <- ggplot(taxgroups_select3tdf_mp12melt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
     labs(title = "Proportional read abundance")+ 
     geom_bar(stat = "identity", position = "stack")+
      facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
     theme(axis.title.y = element_blank(),
           axis.title.x = element_blank())+
     ylim(0,1.01)+
      coord_flip()
    p12
    
    MP12plotly <- ggplotly(p12, tooltip = "text")
    
    p345 <- ggplot(taxgroups_select3tdf_mp345melt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
      geom_bar(stat = "identity", position = "stack")+
      facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      ylim(0,1.01)+
      coord_flip()
    p345
    
    MP345plotly <- ggplotly(p345, tooltip ="text")
    
    ptotal <- ggplot(taxgroups_select3tdf_totalmelt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
      geom_bar(stat = "identity", position = "stack")+
      facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank())+
      ylim(0,1.01)+
      coord_flip()
    ptotal
    
    totalplotly <- ggplotly(ptotal, tooltip = "text")
    
    totalplotly
    
    #Proportional OTU richness
    if (dim(taxgroupspre_bin_sums_yes)[1] < dim(taxgroupspre_bin_sums)[1]) {
      Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups,"Other")
      othersum_pa <- colSums(taxgroups_pa_other[,-1])
    } else {
      Taxonomic_group <- c(taxgroupspre_bin_sums_yes$taxgroups)
      othersum_pa <- NULL
    }
    
    taxgroups_pa_select2 <- rbind(taxgroups_pa_select[,-1],othersum_pa)
    
    otutab_sfsep_pa_num <- otutab_sfsep_pa %>% select_if(is.numeric)
    otutab_sfsep_pa_num_prop <- sweep(otutab_sfsep_pa_num, 2 , colSums(otutab_sfsep_pa_num), FUN = "/")
    
    
    taxgroups_pa_select2_prop <- sweep(taxgroups_pa_select2, 2, colSums(otutab_sfsep_pa_num), FUN = "/")
    
    
   # taxgroups_pa_select2_prop <- sweep(taxgroups_pa_select2, 2, colSums(taxgroups_pa_select2), FUN = "/")
    
    
    
    taxgroups_pa_select3 <- cbind.data.frame(Taxonomic_group,taxgroups_pa_select2_prop)
    
    taxgroups_pa_select3t <- t(taxgroups_pa_select3)
    taxgroups_pa_select3t2 <- taxgroups_pa_select3t[-1,]
    taxgroups_pa_select3t3 <- apply(taxgroups_pa_select3t2, 2, as.numeric)
    taxgroups_pa_select3tdf0 <- as.data.frame(taxgroups_pa_select3t3)
    names(taxgroups_pa_select3tdf0) <- Taxonomic_group
    
    patotal04.3 <- taxgroups_pa_select3 %>% select(samples0.4.3) %>% mutate(rowsumm = rowSums(.))
    patotal3.10 <- taxgroups_pa_select3 %>% select(samples3.10) %>% mutate(rowsumm = rowSums(.))
    patotal10.50 <- taxgroups_pa_select3 %>% select(samples10.50) %>% mutate(rowsumm = rowSums(.))
    patotal50.200<- taxgroups_pa_select3 %>% select(samples50.200) %>% mutate(rowsumm = rowSums(.))
    patotal3.180 <- taxgroups_pa_select3 %>% select(samples3.180) %>% mutate(rowsumm = rowSums(.))
    
    
    patotaltotal04.3 <- otutab_sfsep_pa_num_prop %>% select(samples0.4.3) %>% mutate(rowsumm = rowSums(.))
    patotaltotal3.10 <- otutab_sfsep_pa_num_prop %>% select(samples3.10) %>% mutate(rowsumm = rowSums(.))
    patotaltotal10.50 <- otutab_sfsep_pa_num_prop %>% select(samples10.50) %>% mutate(rowsumm = rowSums(.))
    patotaltotal50.200<- otutab_sfsep_pa_num_prop %>% select(samples50.200) %>% mutate(rowsumm = rowSums(.))
    patotaltotal3.180 <- otutab_sfsep_pa_num_prop %>% select(samples3.180) %>% mutate(rowsumm = rowSums(.))
    
    taxgroups_pa_select3tdf <- rbind(taxgroups_pa_select3tdf0, patotal04.3$rowsumm/sum(patotaltotal04.3$rowsumm), patotal3.10$rowsumm/sum(patotaltotal3.10$rowsumm), patotal10.50$rowsumm/sum(patotaltotal10.50$rowsumm), patotal50.200$rowsumm/sum(patotaltotal50.200$rowsumm), patotal3.180$rowsumm/sum(patotaltotal3.180$rowsumm))
    
    taxgroups_pa_select3tdf$sf <- factor(c(rep("sf0.4.3", npico), rep("sf3.10", nthree10), 
                                           rep("sf10.50", nten50), rep("sf50.200", nfifty200), rep("sf3.180", nthree180), c("sf0.4.3", "sf3.10", "sf10.50", "sf50.200", "sf3.180")), ordered = T, 
                                         levels = c("sf0.4.3", "sf3.10", "sf3.180", "sf10.50", "sf50.200"))
    
    taxgroups_pa_select3tdf$stdep <-  factor(c(as.character(StationDepthsfsep_tb), as.character(rep("total",5))), levels = c(levels(StationDepthsfsep_tb), "total"), ordered = T)
    taxgroups_pa_select3tdf$cruise <-  factor(c(cruisesfsep_tb$cruise, rep("total",5)), ordered = T)
    
    
    taxgroups_pa_select3tdf_mp12 <- taxgroups_pa_select3tdf %>% filter(cruise %in% c("MP1", "MP2"))
    
    taxgroups_pa_select3tdf_mp345 <- taxgroups_pa_select3tdf %>% filter(cruise %in% c("MP3", "MP4", "MP5"))
    
    taxgroups_pa_select3tdf_mp12melt <- melt(taxgroups_pa_select3tdf_mp12)
    taxgroups_pa_select3tdf_mp345melt <- melt(taxgroups_pa_select3tdf_mp345)
    
    taxgroups_pa_select3tdf_total <- taxgroups_pa_select3tdf %>% filter(cruise %in% c("total"))
    taxgroups_pa_select3tdf_totalmelt <- melt(taxgroups_pa_select3tdf_total)
    
    
    p12_pa <- ggplot(taxgroups_pa_select3tdf_mp12melt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
      labs(title = "Proportional OTU richness")+
      geom_bar(stat = "identity", position = "stack")+
      facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      ylim(0,1.01)+
      coord_flip()
    p12_pa
    
    MP12plotly_pa <- ggplotly(p12_pa, tooltip = "text")
    
    p345_pa <- ggplot(taxgroups_pa_select3tdf_mp345melt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
      geom_bar(stat = "identity", position = "stack")+
      facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
      ylim(0,1.01)+
      coord_flip()
    p345_pa
    
    MP345plotly_pa <- ggplotly(p345_pa, tooltip ="text")
    
    ptotal_pa <- ggplot(taxgroups_pa_select3tdf_totalmelt, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
      labs(title = "Proportional OTU richness")+
      geom_bar(stat = "identity", position = "stack")+
      facet_grid(rows = vars(cruise), cols = vars(sf), scale = "free_y", space = "free_y")+
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())+
     # ylim(0,1)+
      coord_flip()
    ptotal_pa
    
    totalplotly_pa <- ggplotly(ptotal_pa, tooltip = "text")
    
    totalplotly_pa
    
    list(taxgroups = taxgroups, MP12plotly = MP12plotly, MP345plotly = MP345plotly, totalplotly = totalplotly, MP12plotly_pa = MP12plotly_pa, MP345plotly_pa = MP345plotly_pa, totalplotly_pa =totalplotly_pa)
    
  })
    
    ###Alpha diversity 
  alpha_div <- reactive({
    if (input$taxo_group_adiv != "All") {
      
      otutab_tax <- otutab_sfsep %>% filter(Divisionlong %in% input$taxo_group_adiv)
    } else {
      otutab_tax <- otutab_sfsep
    }
    
    richtab <- otutab_tax %>% select_if(is.numeric)
    
    
    
    richtab_t <- t(richtab)
    
    richtab_tdf <- as.data.frame(richtab_t)
    
    richtab_tdf$richness <- specnumber(richtab_t)
    richtab_tdf$H <- diversity(richtab_t)
    
    richtab_tdf$N1 <- exp(richtab_tdf$H)
    richtab_tdf$N2 <- diversity(richtab_t, "inv")
    richtab_tdf$E10 <- richtab_tdf$N1/richtab_tdf$richness
    richtab_tdf$E20 <- richtab_tdf$N2/richtab_tdf$richness
    ######
  
    percofmax <- c(NULL)
    percofmax[c(1:npico)] <- 100 * richtab_tdf$richness[c(1:npico)]/max(richtab_tdf$richness[c(1:npico)])
    percofmax[c(45:70)] <- 100 * richtab_tdf$richness[c(45:70)]/max(richtab_tdf$richness[c(45:70)])
    percofmax[c(71:96)] <- 100 * richtab_tdf$richness[c(71:96)]/max(richtab_tdf$richness[c(71:96)])
    percofmax[c(97:122)] <- 100 * richtab_tdf$richness[c(97:122)]/max(richtab_tdf$richness[c(97:122)])
    percofmax[c(123:140)] <- 100 * richtab_tdf$richness[c(123:140)]/max(richtab_tdf$richness[c(123:140)])
    
    richtab_tdf$percofmax <- round(percofmax,0)
    
    
    
    
    #######
    
    richtab_tdf$sample <- names(richtab)
    
    richtab_tdf$sf <- c(rep("sf0.4.3", npico), rep("sf3.10", nthree10), 
                        rep("sf10.50", nten50), rep("sf50.200", nfifty200), rep("sf3.180", nthree180))
    richtab_tdf <- full_join(richtab_tdf, depthbinsfsep_tb, by = "sample")
    
    richtab_tdf$sf <- factor(richtab_tdf$sf, levels = c("sf0.4.3", "sf3.10", "sf10.50", "sf50.200", "sf3.180"), ordered =T)
    
    richtab_tdf$stdep <- StationDepthsfsep_tb
    richtab_tdf$cruise <- factor(cruisesfsep_tb$cruise, ordered = T)
    
    richtab_tdfforplot <- richtab_tdf %>% select(richness, H, N1, N2, E10, E20, sample, sf, depthbin, stdep, cruise, percofmax)
    
    
    #whichvar <- "N1"
    
    
    richtab_tdf_12 <- richtab_tdfforplot %>% filter(cruise %in% c("MP1", "MP2"))
    richtab_tdf_345 <- richtab_tdfforplot %>% filter(cruise %in% c("MP3", "MP4", "MP5"))
    rich12 <- ggplot(richtab_tdf_12, aes(x=cruise, y = richness, colour = depthbin, text = sprintf("Richness: %s<br>Sample: %s<br>Percent of max: %s ", richness, stdep, percofmax)))+
      geom_point( size = bpps)+
      facet_grid(rows = vars(sf))+
      ylim(0,2200)
    rich12
    
    rich12ly <- ggplotly(rich12, tooltip = "text")
    rich12ly
    
    rich345 <- ggplot(richtab_tdf_345, aes(x=cruise, y = richness, colour = depthbin, text = sprintf("Richness: %s<br>Sample: %s<br>Percent of max: %s ", richness, stdep, percofmax)))+
      geom_point(size = bpps)+
      facet_grid(rows = vars(sf))+
      ylim(0,2200)
    rich345
    
    rich345ly <- ggplotly(rich345, tooltip = "text")
    rich345ly
    
    even12 <- ggplot(richtab_tdf_12, aes(x=cruise, y = E10, colour = depthbin, text = sprintf("Evenness: %s<br>Sample: %s ", round(E10,2), stdep)))+
      geom_point(size = bpps)+
      facet_grid(rows = vars(sf))
    
    even12ly <- ggplotly(even12, tooltip = "text")
    even12ly
    
    even345 <- ggplot(richtab_tdf_345, aes(x=cruise, y = E10, colour = depthbin, text = sprintf("Evenness: %s<br>Sample: %s ", E10, stdep)))+
      geom_point(size = bpps)+
      facet_grid(rows = vars(sf))
    
    even345ly <- ggplotly(even345, tooltip = "text")
    even345ly
    
    list(rich12ly = rich12ly, rich345ly = rich345ly, even12ly = even12ly, even345ly = even345ly)
    
  })
  
  ###########
  #Chloroplast data
  prop_tax_chl <- reactive({
    taxlev_chl <- input$taxlevel_chl
    
    if (input$taxo_group3 != "All") {
      
      otutab0_tax_chl <- chltab %>% filter(Division %in% input$taxo_group3)
    } else {
      otutab0_tax_chl <- chltab
    }
    
    
    
    if (is.vector(otutab0_tax_chl)) {
      otutab_tax_chl <- as.data.frame(as.list(otutab0_tax_chl))
    } else {
      otutab_tax_chl <- otutab0_tax_chl}
    
    tabtest <- head(otutab_tax_chl)
    
    otutab_tax_chl$total <- otutab_tax_chl %>% select_if(is.numeric) %>% rowSums()
    
    
    taxlevtab_chl <- otutab_tax_chl %>% select(input$taxlevel_chl) %>% droplevels.data.frame()
    
   
    taxgroupspre_chl <-  select(otutab_tax_chl, -total) %>%   group_by_(input$taxlevel_chl) %>% summarise_if(is.numeric, sum)
     
    
    limfun_chl <- function(x) {
       ifelse(x>=0.0005,1,0)
      }
    # 
  
    taxgroupspre_mat_chl <- as.matrix(taxgroupspre_chl[,-1, drop = FALSE])
    taxgroupspre_bin_chl <- apply(taxgroupspre_mat_chl, 2, limfun_chl)
    
    
    
    # 
    if (is.vector(taxgroupspre_bin_chl)) {
       taxgroupspre_bin2_chl <- as.data.frame(as.list(taxgroupspre_bin_chl))
     } else {
       taxgroupspre_bin2_chl <- as.data.frame(taxgroupspre_bin_chl)}
    # 
    # 
    # 
    taxgroupspre_bin_sums_chl <- taxgroupspre_bin2_chl %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
    taxgroupspre_bin_sums_chl$taxgroups <- as.vector(taxgroupspre_chl[[taxlev_chl]])
    taxgroupspre_bin_sums_yes_chl <- filter(taxgroupspre_bin_sums_chl, rowsumm >0) %>% select(taxgroups)
    
    
    # 
    # 
    taxgroups_select_chl <- taxgroupspre_chl %>% filter(get(input$taxlevel_chl) %in% taxgroupspre_bin_sums_yes_chl$taxgroups)
    
       # 
    # 
    taxgroups_other_chl <- taxgroupspre_chl %>% filter(!get(input$taxlevel_chl) %in% taxgroupspre_bin_sums_yes_chl$taxgroups)
    # 
    # 
    # 
     if (dim(taxgroupspre_bin_sums_yes_chl)[1] < dim(taxgroupspre_bin_sums_chl)[1]) {
       Taxonomic_group_chl <- c(taxgroupspre_bin_sums_yes_chl$taxgroups,"Other")
       othersum_chl <- colSums(taxgroups_other_chl[,-1])
     } else {
       Taxonomic_group_chl <- c(taxgroupspre_bin_sums_yes_chl$taxgroups)
       othersum_chl <- NULL
     }
    # 
    # 
    # 
     taxgroups_select2_chl <- rbind(taxgroups_select_chl[,-1],othersum_chl)
    # 
    # 
    # 
     taxgroups_select3_chl <- cbind.data.frame(Taxonomic_group_chl,taxgroups_select2_chl)
    # 
    # 
     taxgroups_select3t_chl <- t(taxgroups_select3_chl)
     taxgroups_select3t2_chl <- taxgroups_select3t_chl[-1,]
     taxgroups_select3t3_chl <- apply(taxgroups_select3t2_chl, 2, as.numeric)
     taxgroups_select3tdf_chl <- as.data.frame(taxgroups_select3t3_chl)
     names(taxgroups_select3tdf_chl) <- Taxonomic_group_chl
    
     
     
     # 
    # 
    taxgroups_select3tdf_chl$stdep <-  StationDepths_tb_chl
    taxgroups_select3tdf_chl$cruise <-  factor(chlmeta$cruise, ordered = T, levels = c("Jan", "Mar", "May", "Aug", "Nov"))
    taxgroups_select3tdf_chl$nuclac <- factor(chlmeta$nuclac)
    
    # 
    taxgroups_select3tdf_melt_chl <- melt(taxgroups_select3tdf_chl)
    # 
    
    # 
    chloroplot <- ggplot(taxgroups_select3tdf_melt_chl, aes(x=reorder(stdep, desc(stdep)), y = value, fill = variable, text = sprintf("Taxon: %s<br>Sample: %s<br>Percent: %s ", variable, stdep, round(value*100,3))))+
       geom_bar(stat = "identity", position = "stack")+
       facet_grid(rows = vars(cruise), cols = vars(nuclac), scale = "free_y", space = "free_y")+
       theme(axis.title.y = element_blank())+
       coord_flip()
    # 
    # 
     chloroplotly <- ggplotly(chloroplot, tooltip = "text")
    # 
    # 
    list(tabtest = tabtest, chloroplotly = chloroplotly)
  })
  
  
  ###########
  
  sftest <- reactive({
    sf <- sfnames[[input$sizefract1]]
    otutab_sf <- otutab_sfsep[c(as.character(sf), "Divisionlong", "newOTUid_wgen")]
   sampnames_meta_sf <- meta_table[input$sizefract1]
   meta_sf_all <- meta_table[match(colnames(otutab_sf %>% select(-Divisionlong, -newOTUid_wgen)),as.character(sampnames_meta_sf[,1])),]
   
   meta_sf <- droplevels(meta_sf_all)
   
   if (input$taxo_group != "All") {
   
   otutab_sf_tax <- otutab_sf %>% filter(Divisionlong %in% input$taxo_group)
   } else {
     otutab_sf_tax <- otutab_sf
   }
   
   
   
   otutab_sf_tax_num <- otutab_sf_tax %>% select(-Divisionlong, -newOTUid_wgen)
   #hmm <- names(otutab_sf_tax) #"Divisionlong is last column
   numotus <- length(which(rowSums(otutab_sf_tax_num)>0))
   rownames(otutab_sf_tax) <- otutab_sf_tax$newOTUid_wgen
   
   
   #adonis
   if (input$sizefract1 == "sf0.4.3") {
     admp1 <- adonis2(t(otutab_sf_tax_num[,c(1:8)]) ~ depthbin, data = meta_sf[c(1:8),])
     admp2 <- adonis2(t(otutab_sf_tax_num[,c(9:18)]) ~ depthbin, data = meta_sf[c(9:18),])
     admp3 <- adonis2(t(otutab_sf_tax_num[,c(19:28)]) ~ depthbin, data = meta_sf[c(19:28),])
     admp4 <- adonis2(t(otutab_sf_tax_num[,c(29:39)]) ~ depthbin, data = meta_sf[c(29:39),])
     admp5 <- adonis2(t(otutab_sf_tax_num[,c(40:44)]) ~ depthbin, data = meta_sf[c(40:44),])
     admp <- c(admp1$'Pr(>F)'[1], admp2$'Pr(>F)'[1], admp3$'Pr(>F)'[1], admp4$'Pr(>F)'[1], admp5$'Pr(>F)'[1])
   } else { if (input$sizefract1 == "sf3.180") {
     admp1 <- adonis2(t(otutab_sf_tax_num[,c(1:8)]) ~ depthbin, data = meta_sf[c(1:8),])
     admp2 <- adonis2(t(otutab_sf_tax_num[,c(9:18)]) ~ depthbin, data = meta_sf[c(9:18),])
     admp <- c(admp1$'Pr(>F)'[1], admp2$'Pr(>F)'[1])
   } else { if (input$taxo_group == "Stramenopiles_X" && input$sizefract1 == "sf3.10") {
     admp3 <- adonis2(t(otutab_sf_tax_num[,c(1:10)]) ~ depthbin, data = meta_sf[c(1:10),])
     admp4 <- adonis2(t(otutab_sf_tax_num[,c(11:21)]) ~ depthbin, data = meta_sf[c(11:21),])
     admp5 <- adonis2(t(otutab_sf_tax_num[,c(22:25)]) ~ depthbin, data = meta_sf[c(22:25),])
     admp <- c(admp3$'Pr(>F)'[1], admp4$'Pr(>F)'[1], admp5$'Pr(>F)'[1])
   } else {
     admp3 <- adonis2(t(otutab_sf_tax_num[,c(1:10)]) ~ depthbin, data = meta_sf[c(1:10),])
     admp4 <- adonis2(t(otutab_sf_tax_num[,c(11:21)]) ~ depthbin, data = meta_sf[c(11:21),])
     admp5 <- adonis2(t(otutab_sf_tax_num[,c(22:26)]) ~ depthbin, data = meta_sf[c(22:26),])
     admp <- c(admp3$'Pr(>F)'[1], admp4$'Pr(>F)'[1], admp5$'Pr(>F)'[1])
   }}}
   
   
   if (input$taxo_group == "Stramenopiles_X" && input$sizefract1 == "sf3.10") {
     nmds_tax <- mp_nmds(otutab_sf_tax[,-26] %>% select(-Divisionlong, -newOTUid_wgen))
     nmds_table <- data.frame(cbind(nmds_tax$points[,1], nmds_tax$points[,2], meta_sf[-26,]))
   } else {
   nmds_tax <- mp_nmds(otutab_sf_tax %>% select(-Divisionlong, -newOTUid_wgen)) #function mp_nmds transposes otu table
   nmds_table <- data.frame(cbind(nmds_tax$points[,1], nmds_tax$points[,2], meta_sf))
   }
   
   
   
   
   names(nmds_table)[c(1,2)] <- c("nmds_axis1", "nmds_axis2")
   
   #Hierarchical clustering of samples
   braydist <- vegdist(t(otutab_sf_tax_num), distance = "bray")
   otuclust <- hclust(braydist, "average")
  # plotclust <- plot(otuclust)
   
   
   if (input$sizefract1 == "sf0.4.3") {
     Month <-c("January", "March", "May", "August", "November")
   } else { if (input$sizefract1 == "sf3.180") {
   Month <- c("January", "March")
   } else {
   Month <- c("May", "August", "November")}}
   
   
   #mvaotutab <- mvabund(t)
   
   
   #mvabund
   if (input$manyglm == "no") {
     anovaja <- NULL
     anglmmp <- rep("NA",3)
     
     sigmp1 <- NULL
     sigmp2 <- NULL
     sigmp3 <- NULL
     sigmp4 <- NULL
     sigmp5 <- NULL
     #huh <- NULL
   } else {
   if (input$sizefract1 == "sf0.4.3" && length(input$taxo_group)==1) {
     glmmp1 <- manyglm(round(t(otutab_sf_tax[,c(1:8)])*40000,0) ~ meta_sf[c(1:8),]$depthbin, family = "negative.binomial")
     glmmp2 <- manyglm(round(t(otutab_sf_tax[,c(9:18)])*40000,0) ~ meta_sf[c(9:18),]$depthbin, family = "negative.binomial")
     glmmp3 <- manyglm(round(t(otutab_sf_tax[,c(19:28)])*40000,0) ~ meta_sf[c(19:28),]$depthbin, family = "negative.binomial")
     glmmp4 <- manyglm(round(t(otutab_sf_tax[,c(29:39)])*40000,0) ~ meta_sf[c(29:39),]$depthbin, family = "negative.binomial")
     glmmp5 <- manyglm(round(t(otutab_sf_tax[,c(40:44)])*40000,0) ~ meta_sf[c(40:44),]$depthbin, family = "negative.binomial")
     an.glm.mp1 <- anova.manyglm(glmmp1, p.uni = "adjusted")
     an.glm.mp2 <- anova.manyglm(glmmp2, p.uni = "adjusted")
     an.glm.mp3 <- anova.manyglm(glmmp3, p.uni = "adjusted")
     an.glm.mp4 <- anova.manyglm(glmmp4, p.uni = "adjusted")
     an.glm.mp5 <- anova.manyglm(glmmp5, p.uni = "adjusted")
     
     sigmp1 <- cbind.data.frame("OTUid" = attr(an.glm.mp1$uni.p, "dimnames")[[2]], "pval" = an.glm.mp1$uni.p[2,]) %>% filter(pval <= 0.050)
     sigmp2 <- cbind.data.frame("OTUid" = attr(an.glm.mp2$uni.p, "dimnames")[[2]], "pval" = an.glm.mp2$uni.p[2,]) %>% filter(pval <= 0.050)
     sigmp3 <- cbind.data.frame("OTUid" = attr(an.glm.mp3$uni.p, "dimnames")[[2]], "pval" = an.glm.mp3$uni.p[2,]) %>% filter(pval <= 0.050)
     sigmp4 <- cbind.data.frame("OTUid" = attr(an.glm.mp4$uni.p, "dimnames")[[2]], "pval" = an.glm.mp4$uni.p[2,]) %>% filter(pval <= 0.050)
     sigmp5 <- cbind.data.frame("OTUid" = attr(an.glm.mp5$uni.p, "dimnames")[[2]], "pval" = an.glm.mp5$uni.p[2,]) %>% filter(pval <= 0.050)
     
    # huh <- an.glm.mp3$uni.p
     
     anglmmp <- c(an.glm.mp1$table$'Pr(>Dev)'[2], an.glm.mp2$table$'Pr(>Dev)'[2], an.glm.mp3$table$'Pr(>Dev)'[2], 
                  an.glm.mp4$table$'Pr(>Dev)'[2], an.glm.mp5$table$'Pr(>Dev)'[2])
     anovaja <- "see above"
   } else {
     if (input$sizefract1 == "sf3.10" && length(input$taxo_group)==1) {
       glmmp3 <- manyglm(round(t(otutab_sf_tax[,c(1:10)])*40000,0) ~ meta_sf[c(1:10),]$depthbin, family = "negative.binomial")
       glmmp4 <- manyglm(round(t(otutab_sf_tax[,c(11:21)])*40000,0) ~ meta_sf[c(11:21),]$depthbin, family = "negative.binomial")
       glmmp5 <- manyglm(round(t(otutab_sf_tax[,c(22:26)])*40000,0) ~ meta_sf[c(22:26),]$depthbin, family = "negative.binomial")
  
       an.glm.mp3 <- anova.manyglm(glmmp3, p.uni = "adjusted")
       an.glm.mp4 <- anova.manyglm(glmmp4, p.uni = "adjusted")
       an.glm.mp5 <- anova.manyglm(glmmp5, p.uni = "adjusted")
       
       sigmp1 <- NULL
       sigmp2 <- NULL
       sigmp3 <- cbind.data.frame("OTUid" = attr(an.glm.mp3$uni.p, "dimnames")[[2]], "pval" = an.glm.mp3$uni.p[2,]) %>% filter(pval <= 0.050)
       sigmp4 <- cbind.data.frame("OTUid" = attr(an.glm.mp4$uni.p, "dimnames")[[2]], "pval" = an.glm.mp4$uni.p[2,]) %>% filter(pval <= 0.050)
       sigmp5 <- cbind.data.frame("OTUid" = attr(an.glm.mp5$uni.p, "dimnames")[[2]], "pval" = an.glm.mp5$uni.p[2,]) %>% filter(pval <= 0.050)
       
       # huh <- an.glm.mp3$uni.p
       
       anglmmp <- c(an.glm.mp3$table$'Pr(>Dev)'[2], 
                    an.glm.mp4$table$'Pr(>Dev)'[2], an.glm.mp5$table$'Pr(>Dev)'[2])
       anovaja <- "see above"
     } else {
       if (input$sizefract1 == "sf10.50" && length(input$taxo_group)==1) {
         glmmp3 <- manyglm(round(t(otutab_sf_tax[,c(1:10)])*40000,0) ~ meta_sf[c(1:10),]$depthbin, family = "negative.binomial")
         glmmp4 <- manyglm(round(t(otutab_sf_tax[,c(11:21)])*40000,0) ~ meta_sf[c(11:21),]$depthbin, family = "negative.binomial")
         glmmp5 <- manyglm(round(t(otutab_sf_tax[,c(22:26)])*40000,0) ~ meta_sf[c(22:26),]$depthbin, family = "negative.binomial")
         
         an.glm.mp3 <- anova.manyglm(glmmp3, p.uni = "adjusted")
         an.glm.mp4 <- anova.manyglm(glmmp4, p.uni = "adjusted")
         an.glm.mp5 <- anova.manyglm(glmmp5, p.uni = "adjusted")
         
         sigmp1 <- NULL
         sigmp2 <- NULL
         sigmp3 <- cbind.data.frame("OTUid" = attr(an.glm.mp3$uni.p, "dimnames")[[2]], "pval" = an.glm.mp3$uni.p[2,]) %>% filter(pval <= 0.050)
         sigmp4 <- cbind.data.frame("OTUid" = attr(an.glm.mp4$uni.p, "dimnames")[[2]], "pval" = an.glm.mp4$uni.p[2,]) %>% filter(pval <= 0.050)
         sigmp5 <- cbind.data.frame("OTUid" = attr(an.glm.mp5$uni.p, "dimnames")[[2]], "pval" = an.glm.mp5$uni.p[2,]) %>% filter(pval <= 0.050)
         
         # huh <- an.glm.mp3$uni.p
         
         anglmmp <- c(an.glm.mp3$table$'Pr(>Dev)'[2], 
                      an.glm.mp4$table$'Pr(>Dev)'[2], an.glm.mp5$table$'Pr(>Dev)'[2])
         anovaja <- "see above"
     
     } else {
       if (input$sizefract1 == "sf50.200" && length(input$taxo_group)==1) {
         glmmp3 <- manyglm(round(t(otutab_sf_tax[,c(1:10)])*10000,0) ~ meta_sf[c(1:10),]$depthbin, family = "negative.binomial")
         glmmp4 <- manyglm(round(t(otutab_sf_tax[,c(11:21)])*10000,0) ~ meta_sf[c(11:21),]$depthbin, family = "negative.binomial")
         glmmp5 <- manyglm(round(t(otutab_sf_tax[,c(22:26)])*10000,0) ~ meta_sf[c(22:26),]$depthbin, family = "negative.binomial")
         
         an.glm.mp3 <- anova.manyglm(glmmp3, p.uni = "adjusted")
         an.glm.mp4 <- anova.manyglm(glmmp4, p.uni = "adjusted")
         an.glm.mp5 <- anova.manyglm(glmmp5, p.uni = "adjusted")
         
         sigmp1 <- NULL
         sigmp2 <- NULL
         sigmp3 <- cbind.data.frame("OTUid" = attr(an.glm.mp3$uni.p, "dimnames")[[2]], "pval" = an.glm.mp3$uni.p[2,]) %>% filter(pval <= 0.050)
         sigmp4 <- cbind.data.frame("OTUid" = attr(an.glm.mp4$uni.p, "dimnames")[[2]], "pval" = an.glm.mp4$uni.p[2,]) %>% filter(pval <= 0.050)
         sigmp5 <- cbind.data.frame("OTUid" = attr(an.glm.mp5$uni.p, "dimnames")[[2]], "pval" = an.glm.mp5$uni.p[2,]) %>% filter(pval <= 0.050)
         
         # huh <- an.glm.mp3$uni.p
         
         anglmmp <- c(an.glm.mp3$table$'Pr(>Dev)'[2], 
                      an.glm.mp4$table$'Pr(>Dev)'[2], an.glm.mp5$table$'Pr(>Dev)'[2])
         anovaja <- "see above"
       } else {
         if (input$sizefract1 == "sf3.180" && length(input$taxo_group)==1) {
           glmmp1 <- manyglm(round(t(otutab_sf_tax[,c(1:8)])*88000,0) ~ meta_sf[c(1:8),]$depthbin, family = "negative.binomial")
           glmmp2 <- manyglm(round(t(otutab_sf_tax[,c(9:18)])*88000,0) ~ meta_sf[c(9:18),]$depthbin, family = "negative.binomial")
           an.glm.mp1 <- anova.manyglm(glmmp1, p.uni = "adjusted")
           an.glm.mp2 <- anova.manyglm(glmmp2, p.uni = "adjusted")
           
           sigmp1 <- cbind.data.frame("OTUid" = attr(an.glm.mp1$uni.p, "dimnames")[[2]], "pval" = an.glm.mp1$uni.p[2,]) %>% filter(pval <= 0.050)
           sigmp2 <- cbind.data.frame("OTUid" = attr(an.glm.mp2$uni.p, "dimnames")[[2]], "pval" = an.glm.mp2$uni.p[2,]) %>% filter(pval <= 0.050)
           sigmp3 <- NULL
           sigmp4 <- NULL
           sigmp5 <- NULL
           
           # huh <- an.glm.mp3$uni.p
           
           anglmmp <- c(an.glm.mp1$table$'Pr(>Dev)'[2], an.glm.mp2$table$'Pr(>Dev)'[2])
           anovaja <- "see above"
         } else {
     anovaja <- NULL
     anglmmp <- rep("NA",3)
     sigmp1 <- NULL
     sigmp2 <- NULL
     sigmp3 <- NULL
     sigmp4 <- NULL
     sigmp5 <- NULL
     #huh <- NULL
   }}}}}}
   
   
   difftab = as.data.frame(cbind(Month, admp, anglmmp))
   names(difftab) <- c("Month", "p-value adon", "p-value glm")

   ###Plotting NMDS
   plot_nmds <- nmdsplot(nmds_table)+
     xlim(min(nmds_table$nmds_axis1)-.1,max(nmds_table$nmds_axis1)+.1)+
     ylim(min(nmds_table$nmds_axis2)-.1,max(nmds_table$nmds_axis2)+.1)
     
     #coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
   nmdsplotly <- ggplotly(plot_nmds, tooltip = "text")
   
    
    list(nmdsplotly = nmdsplotly, stress = nmds_tax$stress, months = months, 
         difftab = difftab, numotus = numotus, anovaja = anovaja, sigmp1 = sigmp1, sigmp2 = sigmp2,
         sigmp3 = sigmp3, sigmp4 = sigmp4, sigmp5 = sigmp5, otutab = otutab_sf_tax,
         otuclust = otuclust)
  })
  
  #####Rankabundance
  
   rankabund <- reactive({
     
     if (input$taxo_group_ra != "All") {
       
       otutab_ra <- otutab_sfsep %>% filter(Divisionlong %in% input$taxo_group_ra)
     } else {
       otutab_ra <- otutab_sfsep
     }

  dim(otutab_ra)
  otutab_ra_num <- otutab_ra %>% select_if(is.numeric)
  otutab_ra_tax0 <- otutab_ra %>% select_if(negate(is.numeric))
  otutab_ra_tax <- apply(otutab_ra_tax0, 2, as.character)
  
  rownames(otutab_ra_num) <- otutab_ra$newOTUid_wgen 
  rownames(otutab_ra_tax) <- otutab_ra$newOTUid_wgen
  
  #Pick out OTUs > 5% of their taxo group
  #otutab_ra_5perc <- matrix(, nrow = dim(otutab_ra_num)[1], ncol = dim(otutab_ra_num)[2])
  #for (i in c(1:dim(otutab_ra_num)[1])) {
   # for (j in c(1:dim(otutab_ra_num)[2])) {
  #     ifelse(otutab_ra_num[i,j]/sum(otutab_ra_num[,j]) >= 0.05 & !is.na(otutab_ra_num[i,j]/sum(otutab_ra_num[,j])), otutab_ra_5perc[i,j] <- TRUE, otutab_ra_5perc[i,j] <- FALSE)
  #   }}
  # 
  # otunames_xperc <- list()
  # for (i in c(1:ncol(otutab_ra_num))) {
  #   otunames_xperc[[i]] <- rownames(otutab_ra_num)[otutab_ra_5perc[,i]]
  # }
  
 
  
  
  ###
  
  rankmat <- matrix(, nrow = dim(otutab_ra_num)[1], ncol = dim(otutab_ra_num)[2])
  for (i in c(1:dim(otutab_ra_num)[2])) {
    rankmat[,i] <- sort(otutab_ra_num[,i], decreasing = F)
  }
  
  
  
  #Transform into data frame
  rankmatdf <- as.data.frame(rankmat)
  names(rankmatdf) <- names(otutab_ra_num)
  
  rankotus <- matrix(, dim(otutab_ra_num)[1], ncol = dim(otutab_ra_num)[2])
  for (i in c(1:dim(otutab_ra_num)[2])) {
    rankotus[,i] <- rownames(otutab_ra_num[order(otutab_ra_num[,i]),, drop = FALSE])
  }
  
  notus <- 100 #number of OTUs to include in the "rank-abundance barplot"
  
  #Select only the 100 (notu) highest abundance values, and vectorize for ggplot. 
  #Remember that they are sorted in increasing order, so the highest values are at the bottom of the matrix
  rankvals100 <- c(rankmat[c(((dim(rankmat)[1]-notus)+1):dim(rankmat)[1]),])
  
  rankphylum <- matrix(, nrow = notus, ncol = dim(otutab_ra_num)[2])
  for (i in c(1:notus)) {
    for (j in c(1:dim(otutab_ra_num)[2])){
      rankphylum[i,j] <- as.character(otutab_ra_tax[rankotus[dim(otutab_ra)[1]-(notus-i),j],input$taxlevel_ra])
    }
  }
  
  #Vectorise this table for ggplot
  rankphyls100 <- c(rankphylum)
  
  #Group the "Divisionlongs" that do not have a representative among the 5 most abundant OTUs in any sample in "Other" category
  rankphyls100.2 <- character()
  for (i in 1:length(rankphyls100)){
    ifelse(rankphyls100[i] %in% unique(c(rankphylum[c(96:100),])), rankphyls100.2[i] <- rankphyls100[i], rankphyls100.2[i] <- "Other")
  }
  
  #Dummy matrix with sample names for ggplot
  samples <- matrix(, nrow = notus, ncol = dim(otutab_ra_num)[2])
  for (i in c(1:dim(otutab_ra_num)[2])) {
    samples[,i] <- rep(names(otutab_ra_num)[i],notus)
  }
  
  #Vectorise sample names for ggplot
  samplesvec <- c(samples) 
  
  rep100 <- function(x) {
    rep(x, 100)
  }
  
  rankotus_notus <- as.matrix(rankotus[c((dim(otutab_ra_num)[1]-99):dim(otutab_ra_num)[1]),]) #18841:18940
  rankotus_notus2 <- rankotus_notus
  #names(rankotus_notus2) <- names(OTUtab_merged_mean2)
  
  rankotus_notus_vec <- c(rankotus_notus)
  
  
  #Combine abundance values, taxonomic level and dummy variables into data frame
  dfra <- cbind.data.frame(rankvals100, rankphyls100.2, samplesvec, rankotus_notus_vec)
  
  
  ##MP1 and MP2
  dfra_mp12 <- dfra %>% filter(samplesvec %in% c(mp1names, mp2names))
  
  dfra_mp12$stdep <- rep(c(as.character(unlist(lapply(StationDepthsfsep_tb[c(1:nmp1)], rep100))), 
                           as.character(unlist(lapply(StationDepthsfsep_tb[c((nmp1+1):(nmp1+nmp2))], rep100)))), 2)
  
  dfra_mp12$sf <- c(rep("sf0.4.3", 100*(nmp1+nmp2)), rep("sf3.180", 100*(nmp1+nmp2)))
  dfra_mp12$month <- rep(c(rep("Jan", 100*nmp1), rep("Mar", 100*nmp2)), 2)
  
  seqra_12 <- c(rep(c(1:notus),(nmp1+nmp2)))
  
  dfra_mp12$seqra <- seqra_12
  
  dfra_mp12$stdep <- factor(dfra_mp12$stdep, ordered = T, levels =  unique(dfra_mp12$stdep))
  
  #MP3,4,5
  
  dfra_mp345 <- dfra %>% filter(samplesvec %in% c(mp3names, mp4names, mp5names))
  
  dfra_mp345$stdep <- rep(c(as.character(unlist(lapply(StationDepthsfsep_tb[c((nmp1+nmp2+1):(nmp1+nmp2+nmp3))], rep100))), 
                            as.character(unlist(lapply(StationDepthsfsep_tb[c((nmp1+nmp2+nmp3+1):(nmp1+nmp2+nmp3+nmp4))], rep100))), 
                            as.character(unlist(lapply(StationDepthsfsep_tb[c((nmp1+nmp2+nmp3+nmp4+1):(nmp1+nmp2+nmp3+nmp4+nmp5))], rep100)))), 
                          2)
  dfra_mp345$stdep <- factor(dfra_mp345$stdep, ordered = T, levels =  unique(dfra_mp345$stdep))
  
  dfra_mp345$sf <- c(rep("sf0.4.3", 100*(nmp3+nmp4+nmp5)), rep("sf3.10", 100*(nmp3+nmp4+nmp5)), 
                     rep("sf10.50", 100*(nmp3+nmp4+nmp5)), rep("sf50.200", 100*(nmp3+nmp4+nmp5)))
  
  dfra_mp345$sf <- factor(dfra_mp345$sf, ordered = T, levels = unique(dfra_mp345$sf))
  
  dfra_mp345$month <- factor(rep(c(rep("May", 100*nmp3), rep("Aug", 100*nmp4), rep("Nov", 100*nmp5)), 4), 
                             ordered = T, levels = c("May", "Aug", "Nov"))
  
  
  
  seqra_345 <- c(rep(c(1:notus),(nmp3+nmp4+nmp5)))
  
  dfra_mp345$seqra <- seqra_345
  
  ra12 <- ggplot(dfra_mp12, aes(x = reorder(stdep, desc(stdep)), y = rankvals100, fill =rankphyls100.2, group = seqra, text = sprintf("OTU: %s<br>Sample: %s<br>Percent: %s<br>Phylum: %s", rankotus_notus_vec, stdep, round(rankvals100*100,3), rankphyls100.2)))+
    geom_bar(stat = "identity", position = "stack")+
    ylim(0,1)+
    facet_grid(rows = vars(month), cols = vars(sf), scale = "free_y", space = "free_y")+
    theme(axis.title.y = element_blank())+
    coord_flip()
  
  
  ra12ly <- ggplotly(ra12, tooltip = "text")
  
  
  ra345 <- ggplot(dfra_mp345, aes(x = reorder(stdep, desc(stdep)), y = rankvals100, fill =rankphyls100.2, group = seqra, text = sprintf("OTU: %s<br>Sample: %s<br>Percent: %s<br>Phylum: %s", rankotus_notus_vec, stdep, round(rankvals100*100,3), rankphyls100.2)))+
    geom_bar(stat = "identity", position = "stack")+
    ylim(0,1)+
    facet_grid(rows = vars(month), cols = vars(sf), scale = "free_y", space = "free_y")+
    theme(axis.title.y = element_blank())+
    coord_flip()
  
  
  ra345ly <- ggplotly(ra345, tooltip = "text")
  
  list(ra12ly = ra12ly, ra345ly = ra345ly)
   })
  
   
###Betadiversity
   beta_div <- reactive({
     sf <- sfnames[[input$sizefract2]]
     otutab_sf <- otutab_sfsep[c(as.character(sf), "Divisionlong", "newOTUid_wgen")]
     sampnames_meta_sf <- meta_table[input$sizefract2]
     meta_sf_all <- meta_table[match(colnames(otutab_sf %>% select(-Divisionlong, -newOTUid_wgen)),as.character(sampnames_meta_sf[,1])),]
     
     meta_sf <- droplevels(meta_sf_all)
     
     if (input$taxo_group_bdiv != "All") {
       
       otutab_sf_tax <- otutab_sf %>% filter(Divisionlong %in% input$taxo_group_bdiv)
     } else {
       otutab_sf_tax <- otutab_sf
     }
     
     otutab_sf_tax_num <- otutab_sf_tax %>% select_if(is.numeric)
     
     
     beta2 <- beta.div(t(otutab_sf_tax_num), method = "hellinger", nperm = 999)
     summary(beta2)
     
     
     BDtotal <- beta2$beta
     
     beta2$LCBD
     beta2$p.LCBD
     
     p.adj.lcbd <- p.adjust(beta2$p.LCBD, "holm")
     
     betatab_sf_tax <- data.frame(lcbd = beta2$LCBD, p.lcbd = p.adj.lcbd, sample = names(beta2$LCBD))
     
     
     betatab_sf_tax2 <- full_join(betatab_sf_tax, depthbinsfsep_tb, by = "sample") %>%  na.omit()
     
     
     if (input$sizefract2 == "sf0.4.3") {
     stdep <- StationDepthsfsep_tb[1:44]
     cruise <- factor(cruisesfsep_tb$cruise, ordered = T)[1:44]
     } else {
       if (input$sizefract2 == "sf3.10") {
         stdep <- StationDepthsfsep_tb[45:70]
         cruise <- factor(cruisesfsep_tb$cruise, ordered = T)[45:70]
       } else {
         if (input$sizefract2 == "sf10.50") {
           stdep <- StationDepthsfsep_tb[71:96]
           cruise <- factor(cruisesfsep_tb$cruise, ordered = T)[71:96]
         } else {
           if (input$sizefract2 == "sf50.200") {
             stdep <- StationDepthsfsep_tb[97:122]
             cruise <- factor(cruisesfsep_tb$cruise, ordered = T)[97:122]
           } else {
               stdep <- StationDepthsfsep_tb[123:140]
               cruise <- factor(cruisesfsep_tb$cruise, ordered = T)[123:140]
           }
         }
       }
     }
     
     betatab_sf_tax2$cruise <- cruise
     betatab_sf_tax2$stdep <- stdep
     
     str(betatab_sf_tax2)
     #betahmm <- str(cruise)
     
     lcbdplot <- ggplot(betatab_sf_tax2, aes(x=cruise, y = lcbd, colour = depthbin, text = sprintf("LCBD: %s<br>Sample: %s<br>p.val: %s ", lcbd, stdep, p.lcbd)))+
       geom_point( size = bpps)
     
     lcbdplotly <- ggplotly(lcbdplot, tooltip = "text")
     
     list(lcbdplotly = lcbdplotly, BDtotal = BDtotal)
     
     
     
   })
   
   ###manyglm season
   ###Epipelagic
   manyglm_season <- reactive({
     sf <- sfnames[[input$sizefract3]]
     zone <- input$zone
     otutab_sf <- otutab_sfsep[c(as.character(sf), "Divisionlong", "newOTUid_wgen")]
     sampnames_meta_sf <- meta_table[input$sizefract3]
     meta_sf_all <- meta_table[match(colnames(otutab_sf %>% select(-Divisionlong, -newOTUid_wgen)),as.character(sampnames_meta_sf[,1])),]
     
     meta_sf <- droplevels(meta_sf_all)
     
     if (input$taxo_group4 != "All") {
       
       otutab_sf_tax <- otutab_sf %>% filter(Divisionlong %in% input$taxo_group4)
     } else {
       otutab_sf_tax <- otutab_sf
     }
    # str(otutab_sf_tax)
     
     limfun_mvabund <- function(x) {
       ifelse(x>=0.01,1,0)
     }
     
     otutab_sf_tax_num <- otutab_sf_tax %>% select_if(is.numeric)
     otutab_sf_tax_num2 <- sweep(otutab_sf_tax_num, 2, colSums(otutab_sf_tax_num), FUN = "/")
     
     otutab_sf_tax_num2_mat <- as.matrix(otutab_sf_tax_num2)
     otutab_sf_tax_num2_bin <- apply(otutab_sf_tax_num2_mat, 2, limfun_mvabund)
     rownames(otutab_sf_tax_num2_bin) <- otutab_sf_tax$newOTUid_wgen
     
     
     
     if (is.vector(otutab_sf_tax_num2_bin)) {
       otutab_sf_tax_num2_bin2 <- as.data.frame(as.list(otutab_sf_tax_num2_bin))
     } else {
       otutab_sf_tax_num2_bin2 <- as.data.frame(otutab_sf_tax_num2_bin)}
     #str(otutab_sf_tax_num2_bin2)
     
     otutab_sf_tax_num2_bin_sums <- otutab_sf_tax_num2_bin2 %>% select_if(is.numeric) %>% mutate(rowsumm = rowSums(.))
     otutab_sf_tax_num2_bin_sums$otuID <- rownames(otutab_sf_tax_num2_bin2)
     otutab_sf_tax_num2_bin_sums_yes <- filter(otutab_sf_tax_num2_bin_sums, rowsumm >0) %>% select(otuID)
     str(otutab_sf_tax_num2_bin_sums_yes)
     
     otutab_sf_tax_select <- otutab_sf_tax %>% filter(newOTUid_wgen %in% otutab_sf_tax_num2_bin_sums_yes$otuID)
     rownames(otutab_sf_tax_select) <- otutab_sf_tax_select$newOTUid_wgen
     
     
    
   sf_zone <- meta_sf %>% filter(depthbin %in% zone) %>% select(input$sizefract3)
   sf_zone[,1]
   
   otutab_sf_tax_select_zone <- otutab_sf_tax_select[as.character(sf_zone[,1])]
   dim(otutab_sf_tax_select_zone)
   meta_sf_zone <- meta_sf %>% filter(depthbin %in% input$zone)
   #str(otutab_sf_tax_select)
   
   if (input$manyglm2 == "no") {
     glm_tax_month_1_zone = NULL
     an_tax_month_1_zone = NULL
     tax_month_sf_sign_1_zone = NULL
   } else {
  
   glm_tax_month_1_zone <- manyglm(round(t(otutab_sf_tax_select_zone[,c(1:length(sf_zone[,1]))])*40000,0) ~ meta_sf_zone$Cruise, family = "negative.binomial")
   an_tax_month_1_zone <- anova.manyglm(glm_tax_month_1_zone, p.uni = "adjusted")
   tax_month_sf_sign_1_zone <- cbind.data.frame("OTUid" = attr(an_tax_month_1_zone$uni.p, "dimnames")[[2]], "pval" = an_tax_month_1_zone$uni.p[2,]) %>% filter(pval <= 0.050)
   }
   
   list(tax_month_sf_sign_1_zone = tax_month_sf_sign_1_zone)
   })

   ###Indicator values   
   indval <- reactive({
     #sfsamp <- sfnames[[input$sizefract_iva]]
     zone <- input$zone_iva
     tax <- input$taxo_group_iva
     
     sjekk1 <- "test"
     sjekk2 <- input$taxo_group_iva
    # sjekk2 <- get(input$zone_iva)
     
     
     otutab_tax <- otutab_sfsep %>% filter(Divisionlong %in% tax) %>% droplevels()
     # 
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
     iva.tax.zone.list <- list()
     for (i in 1:length(sfnames)) {
       iva.tax.zone.list[[i]] <- multipatt(otutab_tax_num_zone_list[[i]], zone_cruise_list[[i]]$Cruise,
                                           func = "r.g",
                                           max.order = 2,
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
     
     mulpat_forbp_list <- list()
     for (i in 1:length(sfnames)) {
       mulpat_forbp_list[[i]] <- mulpat_tab_list[[i]] %>% group_by(Family) %>% count(index)
     }
     
     season_mulpatplot_list <- list()
     for (i in 1:length(sfnames)) {
       season_mulpatplot_list[[i]] <- ggplotly(ggplot(mulpat_forbp_list[[i]], aes(x=index, y=n, fill = Family, text = sprintf("Taxon: %s<br>n OTUs: %s", Family, n)))+
                                                 geom_bar(stat = "identity", position = "stack")  )
     }
     

          list(iva.tax.zone.list = iva.tax.zone.list, indicplot_tax_list = indicplot_tax_list, season_mulpatplot_list = season_mulpatplot_list, sjekk1 = sjekk1, sjekk2 = sjekk2)
   })
   
   
  
  output$distPlot <- renderPlotly({
   sftest()$nmdsplotly
    })

  
  output$stress <- renderText({
    paste("Stress value: ", round(sftest()$stress,2))
  })
  
  output$plotclust <- renderPlot({
    plot(sftest()$otuclust)
  })
  
  


output$difftab <- renderTable({sftest()$difftab}, 
                           caption = "p-values from adonis test
                           of differential OTU composition between depths.",
                           caption.placement = getOption("xtable.caption.placement", "bottom"), 
                           caption.width = getOption("xtable.caption.width", NULL))
output$numotus <- renderText({ paste("Number of OTUs:", sftest()$numotus)})

output$anovaja <- renderText({sftest()$anovaja})

output$sigmp1 <- renderTable({sftest()$sigmp1})
output$sigmp2 <- renderTable({sftest()$sigmp2})
output$sigmp3 <- renderTable({sftest()$sigmp3})
output$sigmp4 <- renderTable({sftest()$sigmp4})
output$sigmp5 <- renderTable({sftest()$sigmp5})

output$tax_month_sf_sign_1_zone <- renderTable({manyglm_season()$tax_month_sf_sign_1_zone})

output$propPlot1 <- renderPlotly({
  prop_tax()$MP12plotly
})

output$propPlot2 <- renderPlotly({
  prop_tax()$MP345plotly
})

output$propPlotTot <- renderPlotly({
  prop_tax()$totalplotly
})



output$propPlot_pa1 <- renderPlotly({
  prop_tax()$MP12plotly_pa
})

output$propPlot_pa2 <- renderPlotly({
  prop_tax()$MP345plotly_pa
})

output$propPlotTot_pa <- renderPlotly({
  prop_tax()$totalplotly_pa
})



output$bp_tax <- renderPlot({
  prop_tax()$bp_tax
})

output$richsjekk <- renderTable({prop_tax()$richsjekk})

###Chloroplast

 output$propPlotchl <- renderPlotly({
   prop_tax_chl()$chloroplotly
 })

#output$event <- renderPrint({
#  d <- event_data("plotly_hover")
#  if (is.null(d)) "Hover on a point" else d
#})

#Rankabundance

output$raPlot1 <- renderPlotly({
  rankabund()$ra12ly
  
})

output$raPlot2 <- renderPlotly({
  rankabund()$ra345ly
  
})



#Alpha diversity
##Richness
output$richPlot1 <- renderPlotly({
  alpha_div()$rich12ly
  
})

output$richPlot2 <- renderPlotly({
  alpha_div()$rich345ly
  
})

##Evenness
output$evenPlot1 <- renderPlotly({
  alpha_div()$even12ly
  
})

output$evenPlot2 <- renderPlotly({
  alpha_div()$even345ly
  
})

#Betadiv
output$lcbdPlott <- renderPlotly({
  beta_div()$lcbdplotly
})
#output$sjekkmeltchl <- renderTable({prop_tax_chl()$sjekkmeltchl})

output$BDtotal <- renderText({
  paste("BDtotal:", round(beta_div()$BDtotal, 2))})
#Download table with p-values
output$downloadData_glmp3 <- downloadHandler(
  filename = function() {
    paste("pvalues_mp3", ".csv", sep = "")
  },
  content = function(file) {
    write.csv(sftest()$sigmp3, file, row.names = FALSE)
  }
)

output$downloadData_glmp4 <- downloadHandler(
  filename = function() {
    paste("pvalues_mp4", ".csv", sep = "")
  },
  content = function(file) {
    write.csv(sftest()$sigmp4, file, row.names = FALSE)
  }
)

output$downloadData_glmp5 <- downloadHandler(
  filename = function() {
    paste("pvalues_mp5", ".csv", sep = "")
  },
  content = function(file) {
    write.csv(sftest()$sigmp5, file, row.names = FALSE)
  }
)

output$downloadOTUtab <- downloadHandler(
  filename = function() {
    paste("otutab", ".csv", sep = "")
  },
  content = function(file) {
    write.csv(sftest()$otutab, file, row.names = FALSE)
  }
)

output$downloadnmdstab <- downloadHandler(
  filename = function() {
    paste("nmdstab", ".csv", sep = "")
  },
  content = function(file) {
    write.csv(sftest()$nmds_table, file, row.names = FALSE)
  }
)

output$iva <- renderPrint({summary(indval()$iva.tax.zone.list[[1]])})
output$sjekk1 <- renderText({indval()$sjekk1})
output$sjekk2 <- renderText({indval()$sjekk2})
output$iva0.4.3 <- renderPlotly({indval()$season_mulpatplot_list[[1]]})
output$iva3.180 <- renderPlotly({indval()$season_mulpatplot_list[[2]]})
output$iva3.10 <- renderPlotly({indval()$season_mulpatplot_list[[3]]})
output$iva10.50 <- renderPlotly({indval()$season_mulpatplot_list[[4]]})
output$iva50.200 <- renderPlotly({indval()$season_mulpatplot_list[[5]]})

}
# Run the application 
shinyApp(ui = ui, server = server)

