# Metabolic modeling of Auxotrophies in the human gut microbiota

## Project idea

Auxotrophic bacteria are not able for the synthesis of essential nutrients. So they are highly dependent on their nutritional environment. This project is about finding out which kind of auxotrophies can be found in human gut microbiota and what is the influence of auxotrophic gut microbiota on the host. The analyzed genomes are taken from the HRGM dataset including new genomes from Korea, India, Japan. With a  modeling approach, the metabolic function of the auxotrophic microbiota is studied. The production and uptake of nutrients and the association with different kind of auxotrophies is observed.  The following readme file gives an overview about the order of running the scripts. It always starts first with loading the models, predicting the auxotrophies and adding the information about the genomes but depending on the visualization aspects the created dataframes differ. 

## Instructions for running the scripts

For getting the figures as a result, the scripts need to be ran in a specific order. An overview about the running order is given by the flowchart below. Each colorful line represents the order of scripts for every figure. The instructions in the flowchart are given in text below to C+P the code.

![click here for an overview of the order for running the scripts](output/plots/Overview_Order_running_scripts.png)



### 1) Proportions of auxotrophic to prototrophic MAGs  per phylum

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
source("auxotable_not_melted.R")
```

###### Analyze and then visualize the proportions of auxotrophic to prototrophic MAGs per phylum

```R
source("Auxo_Proto_data_per_phylum.R")
```

### 2) Distribution of the number of auxotrophies per phylum
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
source("auxotable_not_melted.R")
```

###### Analyze and visualize the distribution of the number of auxotrophies in every phylum

```R
source("number_auxo_per_phylum.R")
```

### 3) Scatterplot about the correlation of the completeness and the found number of auxotrophies in all genomes

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
source("auxotable_not_melted.R")
```

###### Visualize the correlation of Completeness of Genomes and Number of Auxotrophies

```R
source("Scatterplot_Corr_NumbAuxos_Completeness.R")
```

### 4) Abundancies of amino acid auxotrophies

###### Load models (completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("predict_auxos.R")
```

###### Add information about the genomes

```R
source("auxotable_melted_merged.R")
```

###### Analyze amino acid auxotrophies

```R
source("Abundancies.R")
```

### 5) Fermentation products(e.g.SCFA)

###### Load models (completeness >=85% and a contamination <=2)
```R
source("init_models_filtered.R")
```
###### Predict auxotrophies 
```R
source("predict_auxos.R")
```

###### Add information about the genomes
```R
source("auxotable_melted_merged.R")
```

###### Analyze the production of by products with statistical analysis
```R
source("byproduct_production.R")
```

### 6) Abundancies in the gut

###### Load models (completeness >=85% and a contamination <=2)
```R
source("init_models_filtered.R")
```
###### Predict auxotrophies 
```R
source("predict_auxos.R")
```

###### Add information about the genomes
```R
source("auxotable_melted_merged.R")
```

###### Analyze the abundance of auxotrophies in the gut by data from the FoCus cohorte
```R
source("Abundancies_gut.R")
```

### 7) Occurence of Auxotrophies together

###### Load models() completeness >=85% and a contamination <=2)
```R
source("init_models_filtered.R")
```

###### Predict auxotrophies 
```R
source("predict_auxos.R")
```

###### Add information about the genomes
```R
source("Occurence_Auxos_together.R")
```

### Further scripts that may be used in the future but are yet not displayed in the flowchart

### Completeness of bile acid metabolism pathways

###### Load models and filter them to get only the models with a completeness >=85% and a contamination <=2
```R
source("init_models_filtered.R")
```

###### Predict auxotrophies 
```R
source("predict_auxos.R")
```

###### Add information about the genomes
```R
source("auxotable_melted_merged.R")
```

###### Analyze bile acid metabolism pathways with statistical analysis
```R
source("bile_acid_metabolism_trp_auxo.R")
```

### H2S Production by trp auxotrophic microbiota

###### Load models and filter them to get only the models with a completeness >=85% and a contamination <=2
```R
source("init_models_filtered.R")
```

###### Predict auxotrophies 
```R
source("predict_auxos.R")
```

###### Add information about the genomes
```R
source("auxotable_melted_merged.R")
```

###### Analyze H2S production with statistics
```R
source("H2S_production_trp_auxo_statistics.R")
```



