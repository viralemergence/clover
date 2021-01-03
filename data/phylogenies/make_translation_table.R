# Matching CLOVER harmonized host names to mammal phylogenies

require(ape)

clover <- read.csv("./output/Clover_v1.0_NBCIreconciled_20201218.csv")

# subset to mammals
clover <- clover[clover$HostClass=="Mammalia",]


# Import the Upham et al. node-dated phylogeny for mammals
# Use a single tree for now
# tree_upham <- read.nexus("./vertlife_mammals_node_dated.nex")[[666]]
# write.nexus(tree_upham, file="./upham_tree_666.nex")
tree_upham <- read.nexus("./data/phylogenies/upham_tree_666.nex")

tree_upham$tip.label <- gsub("_"," ", tree_upham$tip.label)

clover$Host_Upham <- clover$Host

# % match
length(intersect(tree_upham$tip.label, clover$Host_Upham))/length(unique(clover$Host_Upham))# 95.8%
setdiff(clover$Host_Upham, tree_upham$tip.label)

# Collapsing subspecies (removes everything after second " ")
# Assumes all names in binomial are separated by " "
collapse.subsp <- function(host)
				{return(sub("(.*? .*?) .*", "\\1", host))}

clover$Host_Upham <- collapse.subsp(clover$Host_Upham)

length(intersect(tree_upham$tip.label, clover$Host_Upham))/length(unique(clover$Host_Upham))# 96.3%
setdiff(clover$Host_Upham, tree_upham$tip.label)

# Using Upham Master Taxonomy to translate NCBI names to IUCN / Upham
IUCN_2_NCBI <- read.csv("./data/phylogenies/Upham_S1_Data/Upham_IUCN_to_NCBI.csv")

IUCN_2_NCBI$NCBI_SciName <- gsub("_"," ", IUCN_2_NCBI$NCBI_SciName)
IUCN_2_NCBI$MasterTax_SciName <- gsub("_"," ", IUCN_2_NCBI$MasterTax_SciName)

lookup <- setNames(IUCN_2_NCBI$MasterTax_SciName, IUCN_2_NCBI$NCBI_SciName)
to_swap <- lookup[setdiff(clover$Host_Upham, tree_upham$tip.label)]
to_swap <- to_swap[!is.na(to_swap)]
clover$Host_Upham=plyr::revalue(clover$Host_Upham, to_swap)

length(intersect(tree_upham$tip.label, clover$Host_Upham))/length(unique(clover$Host_Upham))# 98.1%
setdiff(clover$Host_Upham, tree_upham$tip.label)

# Using ASM Mammal Diversity Database Taxonomy v1.3 (https://www.mammaldiversity.org/)
MDD <- read.csv("./data/phylogenies/MDD_v1.3_6513species.csv")

# First pass mapping from MSW3
MDD$MSW3_sciName <- gsub("_"," ", MDD$MSW3_sciName)
MDD$sciName <- gsub("_"," ", MDD$sciName)

lookup <- setNames(MDD$MSW3_sciName,MDD$sciName)
to_swap <- lookup[setdiff(clover$Host_Upham, tree_upham$tip.label)]
to_swap <- to_swap[!is.na(to_swap)]
clover$Host_Upham=plyr::revalue(clover$Host_Upham, to_swap)

length(intersect(tree_upham$tip.label, clover$Host_Upham))/length(unique(clover$Host_Upham))# 99.4%
setdiff(clover$Host_Upham, tree_upham$tip.label)

# six manual changes
unique(clover[!clover$Host_Upham %in% tree_upham$tip.label, c("Host","HostSynonyms")])
# looking at synonyms and manually searching notes in MDD and Upham MasterTaxonomy

clover$Host_Upham=plyr::revalue(clover$Host_Upham,
                           c("Bubalus carabanensis"="Bubalus arnee", ## B. carabensis is B. bubalis, which is the domestic form of B. arnee
                           	 "Dermanura cinerea"="Dermanura cinereus",
                             "Dermanura tolteca"="Dermanura toltecus",
                             "Lissonycteris angolensis"="Myonycteris angolensis",
                             "Macronycteris commersoni"="Hipposideros commersoni",
                             "Piliocolobus semlikiensis"="Procolobus badius"))# CLOVER synonym is Colobus badius semlikiensis, 
								# MDD reports some Piliocolobus were Colobus,
								# and Upham reports Procolobus badius as extinct MSW3 synonym of Piliocolobus badius

length(intersect(tree_upham$tip.label, clover$Host_Upham))/length(unique(clover$Host_Upham))# 100%

upham_translation_table <- data.frame(Host=clover$Host, Host_Upham=clover$Host_Upham)
write.csv(upham_translation_table,"./data/phylogenies/mammal_phylo_translations.csv")