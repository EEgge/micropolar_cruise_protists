unassign <- read.table(text = gsub(";", "\t", readLines(here("data", "assigned_unclass.txt"))), header = F)

unassign <- read.table(text = gsub(";", "\t", readLines(here("data", "assigned_unclass.txt"))), header = F)
names(unassign) <- c("Accno", "Kingdom", "Supergroup", "Division","Class", "Order", "Family", "Genus", "Species")
genus_unassign <- count(unassign, Genus)
species_dist_unassign <- count(unassign, Species)
