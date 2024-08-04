import os
import pandas as pd
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
# depict venn diagram
venn3(subsets=(757, 1125, 180, 1708, 354, 423, 173), 
      set_labels=('HT', 'HC', 'OB'), 
      set_colors=("orange", "blue", "red"), alpha=0.3)
  
# outline of circle line style and width
venn3_circles(subsets=(757, 1125, 180, 1708, 354, 423, 173),
              linestyle="dashed", linewidth=1)
  
# title of the venn diagram
#plt.title("Venn Diagram in geeks for geeks")
plt.show()
