# Hypothalamic-DNA-Methylation-Aging-Neurodegeneration
This repository contains code from a study on age-related DNA methylation in the hypothalamus of male mice, comparing hippocampus and olfactory brain regions. Key findings include significant methylation differences, involvement of OXT and GnRH pathways, and potential interventions for neurodegenerative diseases like Alzheimer’s.

HC_YvsO_dmr.R : This R script utilizes the methylKit package to identify differentially methylated regions (DMRs) between young and middle-aged mouse hippocampus samples. The script works with methylation call files in text format. 

HCvsOB_12M_DMC.R: This script identifies differentially methylated cytosines (DMCs) between the hippocampus and olfactory bulb of 12-month-old mice based on methylKit package. The script uses SAM files as input to perform the analysis.

HCvsOB_2M_DMC.R : This script identifies differentially methylated cytosines (DMCs) between the hippocampus and olfactory bulb of 2-month aged mice using the SAM files as input based on methylKit package.

HT_YvsO_dmr.R: This R script uses the methylKit package to identify differentially methylated regions (DMRs) between young and middle-aged mouse hypothalamus samples. The script utilizes methylation call files in text format for the analysis.

HTvsHC_12M_DMC.R: This script identifies differentially methylated cytosines (DMCs) between the hypothalamus and hippocampus of 12-month-old mice. The analysis is performed using SAM files as input and methylKit package.

HTvsHC_2M_DMC.R: This R script uses the methylKit package to identify differentially methylated cytosines (DMCs) between the hypothalamus and hippocampus of 2-month-old mice. The analysis is performed using SAM files as input data.

HTvsOB_12M_DMC.R: This R script uses the methylKit package to identify differentially methylated cytosines (DMCs) between the hypothalamus and olfactory bulb of 12-month-old mice. The analysis is conducted using SAM files as input data.

HTvsOB_2M_DMC.R: This R script uses the methylKit package to identify differentially methylated cytosines (DMCs) between the hypothalamus and olfactory bulb of 2-month-old mice. The analysis is conducted using SAM files as input data.

OB_YvsO_dmr.R: This R script uses the methylKit package to identify differentially methylated regions (DMRs) between young and middle-aged mouse olfactory bulb samples. The analysis is conducted using methylation call files in text format.

bar&dot_plot.py: This Python script visualizes KEGG pathway data by generating bar and dot plots. It plots the number of genes involved in each KEGG pathway and the -log10(p-value) to illustrate the significance of each pathway.

lollipop.py: This Python script creates lollipop plots to visualize data. It typically represents a set of values for methylation difference, with each data point displayed as a lollipop with a line and a marker. 

manhattan.py: This Python script generates Manhattan plots to visualize differential methylation differences in brain regions. Specifically, it depicts DMCs with methylation differences greater than 25% between different tissue comparisons. The plot shows each differentially methylated cytosine (DMC) by chromosomal location on the x-axis and -log10(q-value) on the y-axis.

pathway_kegg_dotplot.py: This Python script generates dot plots to visualize KEGG pathway data. It displays pathways with the number of associated genes and the -log10(p-value) to highlight the significance of each pathway.

signaling_dotplot.py: This Python script creates a dot plot to visualize genes involved in OXT and GnRH signaling pathways. The plot displays each gene along with its corresponding methylation values. 

target_dotplot.py: This Python script creates a dot plot to visualize gene methylation data from the target sequening. It focuses on the genes involved in specific treatments and their corresponding methylation values as comapred to vehicle. 

venn.py: This Python script generates Venn diagrams to visualize the overlap and differences between multiple sets of data.

volcano.py: This Python script generates volcano plots to visualize differentially methylated cytosines (DMCs) between different tissue comparisons. The plot displays the methylation differences on the x-axis and -log10(q-values) on the y-axis. Points are colored based on significance thresholds to highlight DMCs with significant changes in methylation. 

