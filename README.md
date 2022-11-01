# Metabolic modeling of Auxotrophies in the human gut microbiota

## Project idea

We already know that humans are auxotrophic for some amino acids. We cannot synthesize these amino acids and cover our demand for these essential amino acids by diet. This incapability of organism to synthesize vital nutrients is called auxotrophy and results in a high dependency on the nutritional environment. However, the nuritional preferences of our human gut microbiome is still not known. This project is about finding out which kind of auxotrophies can be found in human gut microbiota and what is the influence of auxotrophic gut microbiota on the host health. Metagenomes from the HRGM cataologue are used for a general characterzation of auxotrophies in the human gut microbiome. The impact of auxotrophies on microbial ecology in the gut was analyzed with data from two cohorts based in Germany. Metabolic genome-scale models are reconstructed from genomic files with gapseq and the auxotrophies are predicted with metabolic modeling using flux-balance-analysis. The figure below supplies a visualization for the definition of auxotrophies and gives an overview about the workflow with the used datasets in the study.
![image alt ><](/output/plots/Definition_Auxo_readme.png)
img[alt$="><"] {
  display: block;
  max-width: 100%;
  height: auto;
  margin: auto;
  float: none!important;
}


## Instructions for running the scripts

The following README-file gives an overview about the order of running the scripts. It always starts first with the loading of the models. It is recommended to clear the environment before loading a new dataset of models. The instructions in this README includes the analysis for all datasets. The code of every analysis and resulting figures can be found in this repository. An overview about the running order is given by the flowchart below. Each colorful line represents the order of scripts for every figure. Detailed instructions are given under the flowchart. The scripts can be run in R (4.1.2). Further information about the R packages can be found in the first "load models" Scripts. 

![click here for an overview of the order for running the scripts](output/plots/Overview_Order_running_scripts.png)


### 1) Completeness of the amino acid biosynthesis pathways

##### Load models (Completeness >=85%, Contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_not_melted.R")
```

###### Analyze and then visualize the completeness of the amino acid biosynthesis pathways

```R
source("Scripts/Completeness_pathways.R")
```

### 2)Proportions of auxotrophic to prototrophic MAGs  per phylum

##### Load models (Completeness >=85%, Contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_not_melted.R")
```

###### Analyze and then visualize the proportions of auxotrophic to prototrophic MAGs per phylum

```R
source("Scripts/Auxo_Proto_data_per_phylum.R")
```

### 3) Distribution of the number of auxotrophies per phylum

##### Load models (completeness >= 85% and contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_not_melted.R")
```

###### Analyze and visualize the distribution of the number of auxotrophies in every phylum

```R
source("Scripts/number_auxo_per_phylum.R")
```

### 4) Scatterplot about the correlation of the completeness and the found number of auxotrophies in all genomes

##### Load all models

```R
source("Scripts/init_models.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_not_melted.R")
```

###### Visualize the correlation of Completeness of Genomes and Number of Auxotrophies

```R
source("Scripts/Scatterplot_Corr_NumbAuxos_Completeness.R")
```

### 5) Abundancies of amino acid auxotrophies

###### Load models (completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze amino acid auxotrophies

```R
source("Scripts/Abundancies.R")
```

### 6) Fermentation by-products(e.g.SCFA)

###### Load models (completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze the production of by products with statistical analysis

```R
source("Scripts/byproduct_production.R")
```

### 7) Effect of auxotrophic bacteria on the diversity of microbial communities

###### Load models from DZHK cohorte

```R
source("Scripts/DZHK_data_init.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Effect of auxotrophic bateria on the diversity 

```R
source("Scripts/diversity_Auxos.R")
```

### 8) Association of auxotrophies with Health markers in DZHK cohorte

###### Load models from DZHK cohorte

```R
source("Scripts/DZHK_data_init.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze associations with and health markers

```R
source("Scripts/DZHK_Healthmarkers.R")
```

### 9) Abundancies of auxotrophic bacteria in DZHK cohorte

###### Load models

```R
source("Scripts/DZHK_data_init.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Get abundancies 

```R
source("Scripts/Abundancies_Gut_DZHK.R")
```

### 10) Rasch Sampler

###### Load models() completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Rasch sampler

```R
source("Scripts/Rasch Sampler.R")
```

### 11) Occurence of Auxotrophies together

###### Load models() completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/Occurence_Auxos_together.R")### 
```


### 12) Correlation between  the number of auxotrophies and the diversity

###### Load models 

```R
source("Scripts/DZHK_data_init.R.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```
###### Load all scripts from flow number 8 

###### Correlation analysis

```R
source("numb_auxos_div.R")
```

### 13) Crossfeeding observations by calculating hamming distance

###### Load models 

```R
source("Scripts/DZHK_data_init.R.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```
###### Load all scripts from flow number 8 

###### Calculating hamming distance for correlation of diversity and hamming distance

```R
source("hamming_distance.R")
```

### 14) Metabolome levels and frequency of auxotrophic bacteria

###### Load models 

```R
source("Scripts/DZHK_data_init.R.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Association between the frequency of auxotrophic bacteria and serum levels

```R
source("dzhk_metabolome.R")
```

### 15) Influence of auxotrophic bacteria on the stability

###### Load models 

```R
source("Scripts/DZHK_data_init.R.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Auxotrophic bacteria and stability 

```R
source("Popgen_stability.R")
```

### 16) Intake of amino acids on the frequency of amino acid auxotrophic bacteria

###### Load models 

```R
source("Scripts/DZHK_data_init.R.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Auxotrophic bacteria and dietary intake of amino acids

```R
source("AA_intake_popgen_auxos.R")
```

### Further scripts that may be used in the future but are yet not displayed in the flowchart

### Completeness of bile acid metabolism pathways

###### Load models and filter them to get only the models with a completeness >=85% and a contamination <=2

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze bile acid metabolism pathways with statistical analysis

```R
source("Scripts/bile_acid_metabolism_trp_auxo.R")
```

### H2S Production by trp auxotrophic microbiota

###### Load models and filter them to get only the models with a completeness >=85% and a contamination <=2

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze H2S production with statistics

```R
source("Scripts/H2S_production_trp_auxo_statistics.R")
```

### Abundancies in the gut (FoCus cohorte)

###### Load models (completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze the abundance of auxotrophies in the gut by data from the FoCus cohorte

```R
source("Scripts/Abundancies_gut.R")
```

### Associations of auxotrophic bacteria and diseases (FoCus cohorte)

###### Load models (completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze associations with diseases or health factors(BMI, weight, age)

```R
source("Scripts/FoCus_new.R")
```

### Intake of amino acids and abundance of auxotrophic gut bacteria

###### Load models() completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Spearman correlation

```R
source("Scripts/Nutrition_Auxos_FoCus.R")
```

### Frequency of auxotrophic bacteria in chronic inflammatory diseases (TrypCID)

###### Load models (completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/auxotable_melted_merged.R")
```

###### Analyze association with diseases

```R
source("Scripts/TrypCID.R")
```

