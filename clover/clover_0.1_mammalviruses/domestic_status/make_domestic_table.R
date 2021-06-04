

# ====================== Table of domesticated host species ====================

# Source: EID2 and supplemented with additional species (for Gibb et al 2020, Nature)
# Combine into table
load("./source_data/domestic_species_EID2.R")
domestic = data.frame(Host = domestic, IsDomestic = 1, Source = "EID2")
domestic = rbind(domestic, data.frame(Host = c("canis familiaris", 
                                               "bos frontalis", 
                                               "bos grunniens", 
                                               "bos taurus indicus",
                                               "bos taurus primigenius",
                                               "bubalus bubalis",
                                               "bubalus carabanensis",
                                               "lama glama guanicoe",
                                               "vicugna pacos"),
                                      IsDomestic = 1, 
                                      Source = "Gibb"))
domestic$Host = Hmisc::capitalize(as.vector(domestic$Host))
write.csv(domestic, /HostLookup_Domestic.csv", row.names=FALSE)
