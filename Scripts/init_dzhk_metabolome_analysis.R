#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
# Correlations dzhk metabolome and auxotrophy frequencies #
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
library(MicrobiomeGS2)
completeness_cuttoff <- 85
contamination_cutoff <- 2

# load HRGM model data
Metadata <- fread("data/REPR_Genomes_metadata.tsv")

# abundance data
dzhk_hrgm_abun <- t(read.table("data/mgx_abundances/dzhk_abun.csv"))

rel_models <- intersect(names(which(apply(dzhk_hrgm_abun,1, function(x) any(x > 0)))),
                        Metadata[`Completeness (%)`>= completeness_cuttoff & `Contamination (%)` <= contamination_cutoff, `HRGM name`])

dzhk_hrgm_abun <- dzhk_hrgm_abun[rel_models,]

# Predict auxotrophies
dzhk_models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/",
                                      IDs = rel_models)
dzhk_auxos <- predict_auxotrophies(dzhk_models, min.growth = 1e-12, min.growth.fraction = 1e-12)
dzhk_auxos <- lapply(dzhk_auxos, FUN = function(x) {
  tmp <- data.table(aa = names(x),
                    prototrophy = x)
})
dzhk_auxos <- rbindlist(dzhk_auxos, idcol = "model")
dzhk_auxos[is.na(prototrophy), prototrophy := 1]

# Calculate realtive abundancy table
dzhk_relabun <- data.table(as.table(dzhk_hrgm_abun))
setnames(dzhk_relabun, c("model","sample","prop"))
dzhk_relabun <- dzhk_relabun[model %in% rel_models]
dzhk_relabun[, prop := prop/sum(prop), by = sample]


# calculate auxotrophie abundancies
dzhk_auxofreq <- list()
for(spl_i in unique(dzhk_relabun$sample)) {
  tmp <- copy(dzhk_relabun[sample == spl_i])
  tmp <- merge(tmp, dzhk_auxos)
  
  dzhk_auxofreq[[spl_i]] <- tmp[, .(protofreq = sum(prop*prototrophy)), by = .(sample, aa)]
}
dzhk_auxofreq <- rbindlist(dzhk_auxofreq)

# add meta info (DZHK)
dzhk_spl_info <- merge(fread("data/meta/Metagenome_DZHK_NGS_EMGE_sampleID_delete_EMGE173.tsv",
                             header = FALSE, col.names = c("EMGE","sample")),
                       fread("data/meta/DZHK_finaler_export_v5_mod.tsv"),
                       by = "EMGE")

dzhk_auxofreq <- merge(dzhk_auxofreq, dzhk_spl_info[, .(EMGE, sample)])


#––––––––––––––––––––––#
# Load metabolome data #
#––––––––––––––––––––––#
source("Scripts/readBiocrates.R")
library(missForest)
library(tidyxl)
# serum metabolome data
BCdat <- readBiocrates(file = "data/dzhk_metabolome/DZHK_Serum_Ergbnisse- 2021-05-07_UKSH_Serum_Februar2021.xlsx",
                       sheet = "Data Export",
                       incl_val_status = c("Valid", "< LLOQ"))

# impute using random forests
dat_impute <- BCdat$dat[,1:630]
nspls <- nrow(dat_impute)

# Next, we use random forests to impute some missing values. But: Only missing
# values for metabolites, which have less than 20% missing values across the
# data set are imputed.
set.seed(24118)
dat_impute <- dat_impute[,which(apply(dat_impute, 2, function(x) sum(is.na(x))/nspls) < 0.2)]
dat_impute <- missForest(dat_impute, ntree = 100, maxiter = 10)
rm(.Random.seed, envir=globalenv())

dzhk_metabolome <- cbind(BCdat$sample_meas[,.(EMGE = `Sample Identification`)],
                         dat_impute$ximp)
setkey(dzhk_metabolome,"EMGE")

