import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

os.chdir("path_to_directory_containing_files")
df = pd.read_excel("file_name.xlsx", sheet_name='Sheet3')
df

# Set the style of seaborn
sns.set(style="white")

# Create a scatter plot
plt.figure(figsize=(2.5, 15))
scatter = sns.scatterplot(x='Tissue_1', 
                          y='Gene', 
                          hue='Meth.type', 
                          size='Meth.diff', 
                          sizes=(200, 500), 
                          data=df_2,
                          palette={'Hypermethylated': 'green', 'Hypomethylated': 'red'}
                         )

# Add a legend outside the plot
scatter.legend(title='', bbox_to_anchor=(1.2, 1), loc='upper left')
# Add vertical lines
plt.axvline(x=1.5, color='gray', linestyle='--')  # Vertical line at x=1.5
plt.axvline(x=2.5, color='gray', linestyle='--')  # Vertical line at x=2.5

# Add horizontal lines
for gene in df_2['Gene'].unique():
    plt.axhline(y=gene, color='gray', linestyle=':', linewidth=0.5)  # Horizontal line at each gene
    
# Adjust the x-axis limits for better visibility with a bit of space on the left side
#plt.xlim(-0.5,len(df_1['Tissue_1'].unique()) + 0.5)
plt.xlim(0.5, 3.5)
#plt.ylim(-1, len(df['Gene'].unique()) + 0.2)
#Add a little space on the top and bottom
#plt.ylim(df_2['Gene'].min() - 0.5, df_2['Gene'].max() + 0.5)

# Remove plot labels
plt.xlabel('')
plt.ylabel('')

# Show the plot
plt.show()
