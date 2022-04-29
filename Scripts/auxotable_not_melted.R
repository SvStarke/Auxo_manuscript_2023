##################    create a dataframe with genome infos #####################
#merged but not melted
#read metadata file with information about the genomes
Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")

#count = of rows
Auxotrophy$count <- rowSums(Auxotrophy == 0)
o <- unique(Auxotrophy$count)
#merge files
#for analysis of completeness pathways
Auxo_info <- merge(Auxotrophy, Metadata, by.x = "Genomes",
                   by.y = "HRGM name")
#for any other analysis
Auxotrophy_12 <- merge(Auxotrophy, Metadata, by.x = "Genomes",
                       by.y = "HRGM name")
#change column name
names(Auxotrophy_12) [names(Auxotrophy_12) == "Completeness (%)"] <- "Completeness"

Auxotrophy_12$Status <- ifelse(Auxotrophy_12$count == 0, 1, 0)
Auxotrophy_12[, phylum := str_match(`GTDB Taxonomy`, "p__.*;c__")[,1]]
Auxotrophy_12[, phylum := gsub("p__|;c__","", phylum)]


Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Firmicutes_A"] <- "Firmicutes"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Firmicutes_B"] <- "Firmicutes"
Auxotrophy_12$phylum[Auxotrophy_12$phylum== "Firmicutes_G"] <- "Firmicutes"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Firmicutes_I"] <- "Firmicutes"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Bdellovibrionota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Campylobacterota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Cyanobacteria"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Desulfobacterota_A"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Elusimicrobiota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Eremiobacterota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Fibrobacterota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Halobacterota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Myxococcota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Patescibacteria"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Synergistota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Spirochaetota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Thermoplasmatota"] <- "Other"
Auxotrophy_12$phylum[Auxotrophy_12$phylum == "Verrucomicrobiota"] <- "Other"

Auxotrophy_13 <- Auxotrophy_12
Auxotrophy_13[, class := str_match(`GTDB Taxonomy`, "c__.*;o__")[,1]]
Auxotrophy_13[, class := gsub("c__|;o__","", class)]
Auxotrophy_13[, order := str_match(`GTDB Taxonomy`, "o__.*;f__")[,1]]
Auxotrophy_13[, order := gsub("o__|;f__","", order)]
Auxotrophy_13[, family := str_match(`GTDB Taxonomy`, "f__.*;g__")[,1]]
Auxotrophy_13[, family := gsub("f__|;g__","", family)]
Auxotrophy_13[, genus := str_match(`GTDB Taxonomy`, "g__.*;s__")[,1]]
Auxotrophy_13[, genus := gsub("g__|;s__","", genus)]
Auxotrophy_13[, species := str_match(`GTDB Taxonomy`, "s__")[,1]]
Auxotrophy_13[, species := gsub("s__|","", species)]





