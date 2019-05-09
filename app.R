#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(here)
library(ggplot2)
library(tidyverse)
library(vegan)
library(cowplot)
library(ggpubr)

#otutab <- read.table("otutab.txt", header = T, sep = "\t")
#envtab <- read.table("envtab.txt", header = T, sep = "\t")

otutab_sfsep <- read.table(here("data", "OTUtab_nozooembr_minsize5_prop_wtax_wnewseqid_14032019_point.txt"), header = T, sep = "\t")
pico <- readLines(here("data","pico_descnames.txt"))
three10 <- readLines(here("data","three10_descnames.txt"))
three180 <- readLines(here("data","three180_descnames.txt"))
ten50 <- readLines(here("data","ten50_descnames.txt"))
fifty200 <- readLines(here("data","fifty200_descnames.txt"))
net_all <- readLines(here("data","net_all_descnames.txt"))
sfnames <- list("sf0.4.3" = pico, "sf3.10" = three10, "sf10.50" = ten50, "sf50.200" = fifty200, "sf3.180" = three180)

meta_table <- read.table(here("data", "MP_CB_profiles_20170802_wsampnames.txt"), header = T, sep = "\t", fill = NA)

source(here("src", "nmdsplot.R"))
source(here("src", "mp_nmds.R"))

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("NMDS of MicroPolar samples at the OTU level"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     sidebarPanel(
       radioButtons("sizefract", label = "Size fraction", 
                    choices = c("sf0.4.3", "sf3.180",  "sf3.10", "sf10.50", "sf50.200"), selected = NULL),
       checkboxGroupInput("taxo_group", label = "Division", 
                    choices = c("Ochrophyta_Bacillariophyta", "Dinoflagellata_Syndiniales", 
                                "Dinoflagellata_Dinophyceae", "Haptophyta", "Ciliophora", "Chlorophyta", "All"), 
                    selected = "All")
     ),
      
      # Show a plot of the generated distribution
     mainPanel(
       
       #tableOutput("view"),
       plotOutput(outputId = "distPlot"),
       textOutput(outputId = "stress"),
       #textOutput(outputId = "months"),
       tableOutput(outputId = "adon"),
       textOutput(outputId = "numotus")
     )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  sftest <- reactive({
    sf <- sfnames[[input$sizefract]]
    otutab_sf <- otutab_sfsep[c(as.character(sf), "Divisionlong")]
   sampnames_meta_sf <- meta_table[input$sizefract]
   meta_sf_all <- meta_table[match(colnames(otutab_sf[,-dim(otutab_sf)[2]]),as.character(sampnames_meta_sf[,1])),]
   
   meta_sf <- droplevels(meta_sf_all)
   
   if (input$taxo_group != "All") {
   
   otutab_sf_tax <- otutab_sf %>% filter(Divisionlong %in% input$taxo_group)
   } else {
     otutab_sf_tax <- otutab_sf
   }
   
   otutab_sf_tax_num <- otutab_sf_tax %>% select(-Divisionlong)
   #hmm <- names(otutab_sf_tax) #"Divisionlong is last column
   numotus <- length(which(rowSums(otutab_sf_tax_num)>0))
   
   if (input$sizefract == "sf0.4.3") {
     admp1 <- adonis2(t(otutab_sf_tax[,c(1:8)]) ~ depthbin, data = meta_sf[c(1:8),])
     admp2 <- adonis2(t(otutab_sf_tax[,c(9:18)]) ~ depthbin, data = meta_sf[c(9:18),])
     admp3 <- adonis2(t(otutab_sf_tax[,c(19:28)]) ~ depthbin, data = meta_sf[c(19:28),])
     admp4 <- adonis2(t(otutab_sf_tax[,c(29:39)]) ~ depthbin, data = meta_sf[c(29:39),])
     admp5 <- adonis2(t(otutab_sf_tax[,c(40:44)]) ~ depthbin, data = meta_sf[c(40:44),])
     admp <- c(admp1$'Pr(>F)'[1], admp2$'Pr(>F)'[1], admp3$'Pr(>F)'[1], admp4$'Pr(>F)'[1], admp5$'Pr(>F)'[1])
   } else { if (input$sizefract == "sf3.180") {
     admp1 <- adonis2(t(otutab_sf_tax[,c(1:8)]) ~ depthbin, data = meta_sf[c(1:8),])
     admp2 <- adonis2(t(otutab_sf_tax[,c(9:18)]) ~ depthbin, data = meta_sf[c(9:18),])
     admp <- c(admp1$'Pr(>F)'[1], admp2$'Pr(>F)'[1])
   } else {
     admp3 <- adonis2(t(otutab_sf_tax[,c(1:10)]) ~ depthbin, data = meta_sf[c(1:10),])
     admp4 <- adonis2(t(otutab_sf_tax[,c(11:21)]) ~ depthbin, data = meta_sf[c(11:21),])
     admp5 <- adonis2(t(otutab_sf_tax[,c(22:26)]) ~ depthbin, data = meta_sf[c(22:26),])
     admp <- c(admp3$'Pr(>F)'[1], admp4$'Pr(>F)'[1], admp5$'Pr(>F)'[1])
   }}
   
   
   
   nmds_tax <- mp_nmds(otutab_sf_tax[,-dim(otutab_sf_tax)[2]])
   
   nmds_table <- data.frame(cbind(nmds_tax$points[,1], nmds_tax$points[,2], meta_sf))
   
   
   names(nmds_table)[c(1,2)] <- c("nmds_axis1", "nmds_axis2")
   
   if (input$sizefract == "sf0.4.3") {
     Month <-c("January", "March", "May", "August", "November")
   } else { if (input$sizefract == "sf3.180") {
   Month <- c("January", "March")
   } else {
   Month <- c("May", "August", "November")}}
   
   
   adonistab = as.data.frame(cbind(Month, admp))
   names(adonistab) <- c("Month", "p-value")
   
   
   

   
    
    list(nmds_table = nmds_table, stress = nmds_tax$stress, months = months, 
         adonistab = adonistab, numotus = numotus)
  })
  
  
   
  
  
  
  
  
  
  #output$view <- renderTable({ sftest()})
  
  output$distPlot <- renderPlot({
   nmdsplot(sftest()$nmds_table)+
      xlim(min(sftest()$nmds_table$nmds_axis1)-.1,max(sftest()$nmds_table$nmds_axis1)+.1)+
      ylim(min(sftest()$nmds_table$nmds_axis2)-.1,max(sftest()$nmds_table$nmds_axis1)+.1)
    

      
    })
  output$stress <- renderText({
    paste("Stress value: ", round(sftest()$stress,2))
  })
  
  output$months <- renderText({sftest()$months})


output$adon <- renderTable({sftest()$adonistab}, 
                           caption = "p-values from adonis test",
                           caption.placement = getOption("xtable.caption.placement", "bottom"), 
                           caption.width = getOption("xtable.caption.width", NULL))
output$numotus <- renderText({ paste("Number of OTUs:", sftest()$numotus)})

}
# Run the application 
shinyApp(ui = ui, server = server)

