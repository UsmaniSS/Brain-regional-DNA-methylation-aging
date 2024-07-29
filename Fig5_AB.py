import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

df1 = pd.read_excel("Fig_5_oxt_gnrh_data.xlsx", sheet_name='HT_Oxt_gene_promo')
df1
#Prepare data for plotting
my_range = range(1, len(df1.index) + 1)

# Create a new figure
plt.figure(figsize=(3.5, 10))

# Determine colors for points
colors = ['red' if val < 0 else 'green' for val in df1['meth.diff']]

# Plot the horizontal lines for lollipops
plt.hlines(y=my_range, xmin=0, xmax=df1['meth.diff'], color=colors)

# Plot the points
plt.scatter(df1['meth.diff'], my_range, c=colors)

# Add a vertical line at x = 0 for better visualization
plt.axvline(x=0, color='grey', linestyle='--')

# Set the font properties for the tick labels
font_properties = FontProperties(family='Arial', size=8)

# Add titles and axis names
plt.yticks(my_range, df1['Pos'], fontproperties=font_properties)
plt.xticks(range(-30, 21, 10), fontproperties=font_properties)
#plt.title("A vertical lollipop plot", loc='left')
#plt.xlabel('Value of the variable', fontproperties=font_properties)
#plt.ylabel(fontproperties=font_properties)

# Remove gaps on the y-axis by adjusting the limits
plt.ylim(0.5, len(df1.index) + 0.5)

# Show the plot
plt.show()

##### for panel B ####

df = pd.read_excel("Fig_5_oxt_gnrh_data.xlsx", sheet_name='HT_Gnrh_gene_promo')

#Prepare data for plotting
my_range = range(1, len(df.index) + 1)

# Create a new figure
plt.figure(figsize=(3.5, 10))

# Determine colors for points
colors = ['red' if val < 0 else 'green' for val in df['meth.diff']]

# Plot the horizontal lines for lollipops
plt.hlines(y=my_range, xmin=0, xmax=df['meth.diff'], color=colors)

# Plot the points
plt.scatter(df['meth.diff'], my_range, c=colors)

# Add a vertical line at x = 0 for better visualization
plt.axvline(x=0, color='grey', linestyle='--')

# Set the font properties for the tick labels
font_properties = FontProperties(family='Arial', size=8)

# Add titles and axis names
#plt.yticks(my_range, df['Pos'], fontproperties=font_properties)
plt.yticks([i for i in my_range if i % 2 != 0], [df['Pos'][i-1] for i in my_range if i % 2 != 0], fontproperties=font_properties)
plt.xticks(range(-20, 11, 10), fontproperties=font_properties)
#plt.xticks(fontproperties=font_properties)
#plt.title("A vertical lollipop plot", loc='left')
#plt.xlabel('Value of the variable', fontproperties=font_properties)
#plt.ylabel(fontproperties=font_properties)

# Remove gaps on the y-axis by adjusting the limits
plt.ylim(0.5, len(df.index) + 0.5)

# Show the plot
plt.show()

