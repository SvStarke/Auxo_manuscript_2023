metadata_REIMAGINE <- fread("/home/svenja/Downloads/SraRunTable_3.txt")
fwrite(metadata_REIMAGINE, file = "Metadata_Reimagine.csv")
describe(metadata_REIMAGINE_regions$BioSample)
describe(metadata_REIMAGINE_regions$Alias)
describe(metadata_REIMAGINE_regions$`Sample Name`)
View(metadata_REIMAGINE)         


metadata_REIMAGINE_regions <- fread("/home/svenja/Documents/SraRunTable_regions.csv")
describe(metadata_REIMAGINE_regions$Isolation_source)
View(metadata_REIMAGINE_regions)   
