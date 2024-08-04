import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
os.chdir("path_to_directory_containing_files")
df_adcy7_reg2 = pd.read_excel("adcy7_plot.xlsx", sheet_name='adcy7_region2_dotplot')

# Create a new column based on the values of 'existing_column'
df_adcy7_reg2['type'] = df_adcy7_reg2['meth'].apply(lambda x: 'hyper' if x > 0 else 'hypo')

# Function to assign significance based on conditions
def assign_significance(row):
    if 'Oxt+GnRH' in row['treatment']:
        return 'sig'
    elif row['qvalue'] == 0 or row['qvalue'] > 0.05:
        return 'non-sig'
    else:
        return 'sig'

# Applying the function to create the 'significance' column
df_adcy7_reg2['significance'] = df_adcy7_reg2.apply(assign_significance, axis=1)

# Convert a column to text (string) type
df_adcy7_reg2['start'] = df_adcy7_reg2['start'].astype(str)

# Define a function to populate the new column based on 'type' and 'significance'
def new_column_value(row):
    if row['significance'] == 'sig':
        return row['type']
    else:
        return 'no_sig_change'

# Apply the function to create the new column
df_adcy7_reg2['color_type'] = df_adcy7_reg2.apply(new_column_value, axis=1)


# Define a function to populate the new column based on 'meth.abs' and 'significance'
def new_column_value(row):
    if row['significance'] == 'sig':
        return row['meth.abs']
    else:
        return 0.1

# Apply the function to create the new column
df_adcy7_reg2['value_type'] = df_adcy7_reg2.apply(new_column_value, axis=1)
df_adcy7_reg2

sns.set(style="whitegrid")

# Create a scatter plot
plt.figure(figsize=(6, 1.8))
scatter = sns.scatterplot(x='start', y='treatment', hue='color_type', size='value_type', sizes=(20,500),
                          style ='significance', data=df_adcy7_reg2,
                          palette={'hyper': 'green', 'hypo': 'red', 'no_sig_change': 'orange'}, legend=False)

# Add a legend outside the plot
scatter.legend(title='', bbox_to_anchor=(1.03, 1), loc='upper left')

# Adjust the x-axis limits for better visibility
plt.ylim(-0.5, len(df_adcy7_reg2['treatment'].unique()) - 0.5)
#plt.xlim(-1, len(df_adcy5['Site'].unique()) + 0.2)

# Rotate x-axis labels and set font size
plt.xticks(rotation=45, fontsize=12)

# Set font size for y-axis tick labels
plt.yticks(fontsize=12)

# Remove plot labels
plt.xlabel('')
plt.ylabel('')

# Remove axes
#plt.axis('off')

# Show the plot
plt.show()

### similarly, other panels in the Fig 9 was plotted ####
