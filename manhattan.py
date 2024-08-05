import os
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

os.chdir("path_to_directory_containing_files")
df_ht_hip_y_C = pd.read_csv("file_name")

# Define a custom sorting key function
def custom_sort_key(chr_name):
    if chr_name == 'chrX':
        return float('inf') - 1  # Make chrX come second last
    elif chr_name == 'chrY':
        return float('inf')  # Make chrY come last
    return int(chr_name[3:])

# Sort the dataframe using the custom key
df_ht_hip_y_C= df_ht_hip_y_C.sort_values(by='chr', key=lambda col: col.map(custom_sort_key))

# Reset the index for a clean display
df_ht_hip_y_C.reset_index(drop=True, inplace=True)

# Display the sorted dataframe
df_ht_hip_y_C

# Assuming your dataframe is named df_ht_hip_y_C
df_ht_hip_y_C_filtered_hyper = df_ht_hip_y_C[(df_ht_hip_y_C['qvalue'] < 0.05) & (df_ht_hip_y_C['meth.diff'] > 25)]
# Add a new column -log10qvalue
df_ht_hip_y_C_filtered_hyper['-log10qvalue'] = -np.log10(df_ht_hip_y_C_filtered_hyper['qvalue'])

df_ht_hip_y_C_filtered_hyper['ind'] = range(len(df_ht_hip_y_C_filtered_hyper))
df_ht_hip_y_C_filtered_hyper_grouped = df_ht_hip_y_C_filtered_hyper.groupby(('chr'))

# Create the figure and axes
fig = plt.figure(figsize=(5.3, 2.3))
ax = fig.add_subplot(111)

# Plot the data
colors = ['darkblue', 'orange']
x_labels = []
x_labels_pos = []

for num, (name, group) in enumerate(df_ht_hip_y_C_filtered_hyper_grouped):
    group.plot(kind='scatter', x='ind', y='-log10qvalue', color=colors[num % len(colors)], ax=ax, s=10, alpha=1)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

# Set the x-ticks and labels
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df_ht_hip_y_C_filtered_hyper)])
ax.set_xlabel('Chromosome')

# Set the y-axis ticks and range
y_ticks = range(0, 35, 10)  # Define the y-axis tick positions
ax.set_yticks(y_ticks)
ax.set_ylim([0, 35])  # Set the y-axis range

# Remove the frame by hiding the spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Optionally, adjust the appearance of the x and y axes
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

# Show the plot
plt.show()

### lower panel in figures, for hypomethylation ####
# Assuming your dataframe is named df_ht_hip_y_C
df_ht_hip_y_C_filtered_hypo = df_ht_hip_y_C[(df_ht_hip_y_C['qvalue'] < 0.05) & (df_ht_hip_y_C['meth.diff'] < -25)]
# Add a new column +log10qvalue
df_ht_hip_y_C_filtered_hypo['+log10qvalue'] = np.log10(df_ht_hip_y_C_filtered_hypo['qvalue'])

df_ht_hip_y_C_filtered_hypo['ind'] = range(len(df_ht_hip_y_C_filtered_hypo))
df_ht_hip_y_C_filtered_hypo_grouped = df_ht_hip_y_C_filtered_hypo.groupby(('chr'))
# Create the figure and axes
fig = plt.figure(figsize=(5.3, 2.3))
ax = fig.add_subplot(111)

# Plot the data
colors = ['darkblue', 'orange']
x_labels = []
x_labels_pos = []

for num, (name, group) in enumerate(df_ht_hip_y_C_filtered_hypo_grouped):
    group.plot(kind='scatter', x='ind', y='+log10qvalue', color=colors[num % len(colors)], ax=ax, s=10, alpha=1)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

# Set the x-ticks and labels
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df_ht_hip_y_C_filtered_hypo)])
ax.set_xlabel('Chromosome')

# Set the y-axis ticks and range
y_ticks = range(0, -35, -10)  # Define the y-axis tick positions
ax.set_yticks(y_ticks)
ax.set_ylim([-35, 0])  # Set the y-axis range

# Remove the frame by hiding the spines
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)

# Optionally, adjust the appearance of the x and y axes
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('top')

# Show the plot
plt.show()
