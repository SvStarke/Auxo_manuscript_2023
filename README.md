# Amino acid auxotrophies in human gut bacteria are linked to higher microbiome diversity and long-term stability

## Project idea

We already know that humans are incapable to synthesize some amino acids. We cover our demand for these essential amino acids by diet. This incapability of organism to synthesize essential nutrients is called auxotrophy and results in a high dependency on the nutritional environment. However, little is still known  about the nuritional preferences of our human gut microbiome. This project is about the abundance of amino acid auxotrophies in the human gut microbiota and their influence on microbial ecology. 
Genomes from the HRGM catalogue (doi: 10.1186/s13073-021-00950-7) are used for general characterization of auxotrophies in the human gut microbiome. The impact of auxotrophies on microbial ecology in the gut was analyzed with data from three cohorts. Metabolic models are reconstructed from genomes with gapseq (doi: 10.1186/s13059-021-02295-1) and the auxotrophies are predicted with genome-scale metabolic modeling using flux-balance-analysis. The figure below supplies a visualization for the definition of auxotrophies (A) and gives an overview of the study workflow (B).

<p align="center" width = "100%" >
  <img width = "100%" src="/output/plots/Definition_Auxo_readme.png" />
</p>

Free available icons were taken from www.flaticon.com (creators: photo3idea_studio, Freepik, surang, Eucalyp, Voysla, juicy_fish, smashingstocks, SBTS2018).

## Instructions for running the scripts

The following README-file gives an overview about the running order of scripts. First, information about the metadata is loaded. It is recommended to clear the environment before loading a new dataset. The instructions in this README includes the analysis for all datasets. The code of every analysis and resulting figures can be found in this repository. An overview about the running order is given by the flowchart below. Each colorful line represents the order of scripts for every figure. Detailed instructions are given under the flowchart. The scripts can be run in R (4.1.2). Further information about the R packages can be found in the first "load models" Scripts. 

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

###### Analyze and visualize the completeness of the amino acid biosynthesis pathways

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

###### Analyze and  visualize the proportions of auxotrophic to prototrophic MAGs per phylum

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

### 5) Abundancies of amino acid auxotrophies (HRGM)

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

###### Analyze and visualize the abundance of amino acid auxotrophies

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

###### Analyze and visualize the production of by products with statistical analysis

```R
source("Scripts/byproduct_production.R")
```

### 7) Effect of auxotrophic bacteria on the diversity of microbial communities

###### Load models from discovery cohorte

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

### 8) Association of auxotrophies with Health markers in discovery cohorte

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

###### Analyze and visualize associations health markers

```R
source("Scripts/DZHK_Healthmarkers.R")
```

### 9) Frequencies of auxotrophic bacteria in discovery cohorte

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

###### Get frequencies and visualize  

```R
source("Scripts/Abundancies_Gut_DZHK.R")
```

### 10) Rasch Sampler

###### Load models(completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```
###### Analyze and visualize co-occurencies of auxotrophies

```R
source("Scripts/Occurence_Auxos_together.R") 
```

###### Rasch sampler analysis

```R
source("Scripts/Rasch_Sampler.R")
```

### 11) Occurence of Auxotrophies together

###### Load models(completeness >=85% and a contamination <=2)

```R
source("Scripts/init_models_filtered.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Analyze and visualize co-occurencies of auxotrophies

```R
source("Scripts/Occurence_Auxos_together.R") 
```


### 12) Correlation between  the number of auxotrophies and the diversity

###### Load models 

```R
source("Scripts/DZHK_data_init.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```
###### Load all scripts from flow number 7 

###### Correlation analysis

```R
source("numb_auxos_div.R")
```

### 13) Crossfeeding observations by calculating hamming distance

###### Load models 

```R
source("Scripts/DZHK_data_init.R")
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

###### Association between the frequency of auxotrophic bacteria and serum metabolite levels

```R
source("dzhk_metabolome.R")
```

### 15) Influence of auxotrophic bacteria on  stability

The workflow in the following script calculates the distribution of auxotrophies in the two longitudinal cohorts (Chen *et al.* 2021 (a.k.a. Lifelines), and *Troci et al.* 2022(a.k.a. Popgen)).

```R
source("stability.R")
```

### 16) Intake of amino acids on the frequency of amino acid auxotrophic bacteria

###### Load models 

```R
source("Scripts/Popgen_data_init.R")
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
### 17) Hamming distance and stability (Popgen Cohorte)

###### Load models 

```R
source("Scripts/Popgen_data_init.R")
```

###### Predict auxotrophies

```R
source("Scripts/predict_auxos.R")
```

###### Add information about the genomes

```R
source("Scripts/Popgen_stability.R")
```

###### Auxotrophic bacteria and dietary intake of amino acids

```R
source("Scripts/Distance_Stability.R")
```

