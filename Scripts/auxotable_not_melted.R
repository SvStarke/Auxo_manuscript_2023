##################    create a dataframe with genome infos #####################
#merged but not melted
#read metadata file with information about the genomes
Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")

#count = of rows
Auxotrophy$count <- rowSums(Auxotrophy == 0)

#merge files
Auxotrophy_12 <- merge(Auxotrophy, Metadata, by.x = "Genomes",
                       by.y = "HRGM name")
#change column name
names(Auxotrophy_12) [names(Auxotrophy_12) == "Completeness (%)"] <- "Completeness"