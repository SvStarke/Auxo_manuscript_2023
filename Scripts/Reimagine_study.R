metadata_REIMAGINE <- fread("/home/svenja/Downloads/SraRunTable_3.txt")
fwrite(metadata_REIMAGINE, file = "Metadata_Reimagine.csv")
describe(metadata_REIMAGINE$BioSample)
describe(metadata_REIMAGINE$Alias)
describe(metadata_REIMAGINE$`Sample Name`)
View(metadata_REIMAGINE)         


metadata_REIMAGINE_regions <- fread("/home/svenja/Documents/SraRunTable_regions.csv")
describe(metadata_REIMAGINE_regions$Isolation_source)
View(metadata_REIMAGINE_regions)   
