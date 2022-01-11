# Metabolic modeling of Auxotrophies in the human gut microbiota

## Project idea

Auxotrophic bacteria are not able for the synthesis of essential nutrients. So they are highly dependent on their nutritional environment. This project is about finding out which kind of auxotrophies can be found in human gut microbiota and what is the influence of auxotrophic gut microbiota on the host. The analyzed genomes are taken from the HRGM dataset including new genomes from Korea, India, Japan. With a  modeling approach, the metabolic function of the auxotrophic microbiota is studied. The production and uptake of nutrients and the association with different kind of auxotrophies is observed.  The following readme file gives an overview about the order of running the scripts. It always starts first with loading the models, predicting the auxotrophies and adding the information about the genomes but depending on the visualization aspects the created dataframes differ. 

## Instructions for running the scripts
![alt text] (https://cau-git.rz.uni-kiel.de/AEF/nutriinformatik/svenja/auxotrophies_hrgm/-/blob/75adb296b928357ead5676ee3a47ea4dce2d825a/output/plots/Overview_Order_running_scripts.png)

### For getting the scatterplot(Completeness, Number Auxos) as a result

###### Load models

init_models.R

###### Predict auxotrophies

predict_auxos.R

###### Add information about the genomes

auxotable_not_melted.R

###### Visualize the correlation of Completeness of Genomes and Number of Auxotrophies

Scatterplot_Corr_NumbAuxos_Completeness.R



### Abundancies of amino acid auxotrophies

###### Load models and filter them to get only the models with a completeness >=90% and a contamination <=2

init_models_filtered.R

###### Predict auxotrophies 

predict_auxos.R

###### Add information about the genomes

auxotable_melted_merged.R

###### Analyze amino acid auxotrophies

Abundancies_AA_auxotrophies_HRGM.R



### SCFA production by trp auxotrophic bacteria

###### Load models and filter them to get only the models with a completeness >=90% and a contamination <=2

init_models_filtered.R

###### Predict auxotrophies 

predict_auxos.R

###### Add information about the genomes

auxotable_melted_merged.R

###### Analyze SCFA production with statistical analysis

SCFA_production_trp_auxo_statistics.R



### Completeness of bile acid metabolism pathways

###### Load models and filter them to get only the models with a completeness >=90% and a contamination <=2

init_models_filtered.R

###### Predict auxotrophies 

predict_auxos.R

###### Add information about the genomes

auxotable_melted_merged.R

###### Analyze bile acid metabolism pathways with statistical analysis

bile_acid_metabolism_trp_auxo.R



### H2S Production by trp auxotrophic microbiota

###### Load models and filter them to get only the models with a completeness >=90% and a contamination <=2

init_models_filtered.R

###### Predict auxotrophies 

predict_auxos.R

###### Add information about the genomes

auxotable_melted_merged.R

###### Analyze bile acid metabolism pathways with statistical analysis

H2S_production_trp_auxo_statistics.R





