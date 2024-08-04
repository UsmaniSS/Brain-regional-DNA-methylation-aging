import os 
import pandas as pd
import numpy as np
os.chdir("path_to_directory_contaning_files")
df_ht = pd.read_table("file_name.txt", sep='\t')
df_ht_sig = df_ht.loc[df_ht['P-value']<= 0.05]
df_ht_sig['Number of genes'] = df_ht_sig['Genes'].str.split(';').map(len)
df_ht_sig['-log10_pvalue'] = -(np.log10(df_ht_sig['P-value']))
import plotly.express as px
fig = px.scatter(df_ht_sig, 
                 x="-log10_pvalue", 
                 y="Term", 
                 color="Number of genes", 
                 size='Combined Score')

fig.update_layout(height=400, width=600)  # Set the height and width of the plot

# Increase font size of axis labels and tick marks
fig.update_xaxes(title_font=dict(size=12), tickfont=dict(size=8))  # Set x-axis label and tick font size
fig.update_yaxes(title_font=dict(size=12), tickfont=dict(size=8))  # Set y-axis label and tick font size

# Update grid options
#fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')  # Show and customize x-axis grid
#fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')  # Show and customize y-axis grid

fig.show()

#### for panel B#### 
df_hip = pd.read_table("file_name.txt", sep='\t')
df_hip_sig = df_hip.loc[df_hip['P-value']<= 0.05]
df_hip_sig['Number of genes'] = df_hip_sig['Genes'].str.split(';').map(len)
df_hip_sig['-log10_pvalue'] = -(np.log10(df_hip_sig['P-value']))
df_hip_sig
fig = px.scatter(df_hip_sig, 
                 x="-log10_pvalue",
                 y="Term", 
                 color="Number of genes", 
                 size='Combined Score')
fig.update_layout(height=400, width=600)  # Set the height and width of the plot

# Increase font size of axis labels and tick marks
fig.update_xaxes(title_font=dict(size=12), tickfont=dict(size=8))  # Set x-axis label and tick font size
fig.update_yaxes(title_font=dict(size=12), tickfont=dict(size=6))  # Set y-axis label and tick font size

# Update grid options
#fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')  # Show and customize x-axis grid
#fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')  # Show and customize y-axis grid

fig.show()

#### for panel C ####
df_ob= pd.read_table("Kfile_name.txt", sep='\t')
df_ob_sig = df_ob.loc[df_ob['P-value']<= 0.05]
df_ob_sig['Number of genes'] = df_ob_sig['Genes'].str.split(';').map(len)
df_ob_sig['-log10_pvalue'] = -(np.log10(df_ob_sig['P-value']))
df_ob_sig
fig = px.scatter(df_ob_sig, 
                 x="-log10_pvalue", 
                 y="Term", 
                 color="Number of genes", 
                 size='Combined Score')
fig.update_layout(height=400, 
                  width=650)  # Set the height and width of the plot

# Increase font size of axis labels and tick marks
fig.update_xaxes(title_font=dict(size=12), tickfont=dict(size=8))  # Set x-axis label and tick font size
fig.update_yaxes(title_font=dict(size=12), tickfont=dict(size=8))  # Set y-axis label and tick font size

# Update grid options
#fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')  # Show and customize x-axis grid
#fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')  # Show and customize y-axis grid

fig.show()



