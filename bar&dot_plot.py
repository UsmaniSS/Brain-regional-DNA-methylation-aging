import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Data
data = {
    'Circadian rhythm': (4, 0.0001243),
    'Cholinergic synapse': (4, 0.0168),
    'Glutamatergic synapse': (4, 0.01729),
    'GABAergic synapse': (3, 0.04326),
    'Serotonergic synapse': (4, 0.02782),
    'Dopaminergic synapse': (4, 0.02988),
    'GnRH signaling pathway': (4, 0.007759),
    'Oxytocin signaling pathway': (8, 0.00005953),
    'Aldosterone synthesis and secretion': (6, 0.0002604),
    'Cortisol synthesis and secretion': (3, 0.02197),
    'Thyroid hormone synthesis': (3, 0.02544),
    'Parathyroid hormone synthesis, secretion and action': (5, 0.002388),
    'Aldosterone-regulated sodium reabsorption': (2, 0.04268),
    'Insulin secretion': (3, 0.03863),
    'Nicotinate and nicotinamide metabolism': (2, 0.03867),
    'Choline metabolism in cancer': (5, 0.001698),
    'Proteoglycans in cancer': (5, 0.03187),
    'Calcium signaling pathway': (5, 0.02442),
    'Gastric acid secretion': (3, 0.02635),
    'Salivary secretion': (3, 0.03017),
    'Hypertrophic cardiomyopathy (HCM)': (5, 0.0009044),
    'Dilated cardiomyopathy (DCM)': (5, 0.00111),
    'Arrhythmogenic right ventricular cardiomyopathy (ARVC)': (3, 0.02455),
    'Vascular smooth muscle contraction': (4, 0.0335),
    'Adrenergic signaling in cardiomyocytes': (4, 0.0398),
    'Amphetamine addiction': (3, 0.02115),
    'Morphine addiction': (3, 0.04568),
    'Fc gamma R-mediated phagocytosis': (4, 0.006895),
    'Inflammatory mediator regulation of TRP channels': (4, 0.02459),
    'Tight junction': (6, 0.00336),
    'Adherens junction': (4, 0.003521),
    'MAPK signaling pathway': (8, 0.004153),
    'Gap junction': (4, 0.006622),
    'Phospholipase D signaling pathway': (5, 0.009624),
    'Rap1 signaling pathway': (6, 0.009793),
    'Glioma': (3, 0.02728),
    'Retrograde endocannabinoid signaling': (4, 0.04148)
    }

# Extract data for plotting
genes = list(data.keys())
numbers = [d[0] for d in data.values()]
p_values = [d[1] for d in data.values()]
neg_log_p_values = [-np.log10(p) for p in p_values]

# Create figure and axes with adjusted layout
fig, ax = plt.subplots(figsize=(12, 6))

# Set minimum gap between y-axis and the first bar
#ax.set_xlim(left=-0.5)

# Plotting the bar plot for number of occurrences
ax.bar(genes, numbers, width = 0.4, color='blue', alpha=0.5, label='Number of genes')

# Create a secondary y-axis
ax2 = ax.twinx()

# Plotting the scatter plot for negative logarithm on the secondary y-axis
ax2.scatter(genes, neg_log_p_values, color='red', marker='o', label='-log(p-value)')

# Adding labels and title
ax.set_ylabel('Number of genes')
ax2.set_ylabel('-log(p-value)')
ax.set_title('Gene Occurrences and Negative Log')

# Remove x-axis labels
ax.set_xticklabels([])  ### it can be done here, but I used ppt for better size and fit on the fifure panel. ####

# Show legend
ax.legend(loc='upper left')
ax2.legend(loc='upper right')

# Adjust layout to make room for x-axis labels
plt.tight_layout()

# Display the plot
plt.show()    
