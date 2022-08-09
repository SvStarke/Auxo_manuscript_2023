#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
# Correlations dzhk metabolome and auxotrophy frequencies #
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
library(MicrobiomeGS2)
completeness_cuttoff <- 85
contamination_cutoff <- 2

# get QC for MAGs
dzhk_MAGQC <- fread("/mnt/nuuk/2022/DZHK_MGX/atlas/completeness.tsv")
dzhk_rel_mags <- dzhk_MAGQC[Completeness >= completeness_cuttoff & contamination_cutoff <= contamination_cutoff, `Bin Id`]

# Predict auxotrophies
dzhk_models <- fetch_model_collection("/mnt/nuuk/2022/DZHK_MGX/models/",
                                      IDs = dzhk_rel_mags)
dzhk_auxos <- predict_auxotrophies(dzhk_models)
dzhk_auxos <- lapply(dzhk_auxos, FUN = function(x) {
  tmp <- data.table(aa = names(x),
                    prototrophy = x)
})
dzhk_auxos <- rbindlist(dzhk_auxos, idcol = "model")
dzhk_auxos[is.na(prototrophy), prototrophy := 1]

# load median coverage of MAGs
dzhk_mag_abun <- fread("/mnt/nuuk/2022/DZHK_MGX/atlas/median_coverage_genomes.tsv")
setnames(dzhk_mag_abun, "V1","sample")
dzhk_mag_abun <- melt(dzhk_mag_abun, id.vars = "sample", variable.name = "model", value.name = "cov")
dzhk_mag_abun <- dzhk_mag_abun[model %in% dzhk_rel_mags]
dzhk_mag_abun[, cov.norm := cov/sum(cov), by = "sample"]

# calculate auxotrophie abundancies
dzhk_auxofreq <- list()
for(spl_i in unique(dzhk_mag_abun$sample)) {
  tmp <- copy(dzhk_mag_abun[sample == spl_i])
  tmp <- merge(tmp, dzhk_auxos)
  
  dzhk_auxofreq[[spl_i]] <- tmp[, .(protofreq = sum(cov.norm*prototrophy)), by = .(sample, aa)]
}
dzhk_auxofreq <- rbindlist(dzhk_auxofreq)

# add meta info (DZHK)
dzhk_spl_info <- merge(fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/Metagenome_DZHK_NGS_EMGE_sampleID_delete_EMGE173.tsv",
                             header = FALSE, col.names = c("EMGE","sample")),
                       fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/DZHK_finaler_export_v5_mod.csv"),
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
