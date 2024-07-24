import os
import pandas as pd
import bioinfokit
from bioinfokit import analys, visuz

### code for plotting Fig1A left panel i.e. volcano plot of DMCs, HT vs HC, 2 Months ########

os.chdir("path_to_working_directory")
df_ht_hip_y_CG = pd.read_csv("path_to_DMCG_directory/myDiff.csv")
df_ht_hip_y_CHG = pd.read_csv("path_to_DMCHG_directory/myDiff.csv")
df_ht_hip_y_CHH = pd.read_csv("path_to_DMCHH_directroy/myDiff.csv")

frames = [df_ht_hip_y_CG, df_ht_hip_y_CHG, df_ht_hip_y_CHH]
df_ht_hip_y_C = pd.concat(frames)
df_ht_hip_y_C.shape

visuz.GeneExpression.volcano(df=df_ht_hip_y_C,
                             lfc='meth.diff',
                             lfc_thr=(25, 25),
                             pv='qvalue',
                             pv_thr=(0.05, 0.05), 
                             dotsize=10,
                             sign_line=True,
                             axtickfontsize=12,
                             axlabelfontsize=15,
                             axxlabel= 'Methylation difference (%)',
                             axylabel='-log10(qvalue)',
                             figname='HT_vs_Hip_Y_C'
                            )

### code for plotting Fig1A right panel i.e. volcano plot of DMCs, HT vs HC, 12 Months ########
### change the path to directory, having DMCs results file between HT vs HC, 12 Months #####
os.chdir("path_to_working_directory")

df_ht_hip_o_CG = pd.read_csv("path_to_DMCG_directory/myDiff.csv")
df_ht_hip_o_CHG = pd.read_csv("path_to_DMCHG_directory/myDiff.csv")
df_ht_hip_o_CHH = pd.read_csv("path_to_DMCHH_directroy/myDiff.csv")

frames = [df_ht_hip_o_CG, df_ht_hip_o_CHG, df_ht_hip_o_CHH]
df_ht_hip_o_allC = pd.concat(frames)
df_ht_hip_o_allC.shape

visuz.GeneExpression.volcano(df=df_ht_hip_o_allC,
                             lfc='meth.diff',
                             lfc_thr=(25, 25), 
                             pv='qvalue', 
                             pv_thr=(0.05, 0.05), 
                             dotsize=10, 
                             sign_line=True, 
                             axtickfontsize=12, 
                             axlabelfontsize=15,
                             axxlabel= 'Methylation difference (%)', 
                             axylabel='-log10(qvalue)',
                             figname='HT_vs_Hip_o_allC'
                            )

### in the same way, plotted rest of the volcano plots, reading the correponsing DMCs, and using above libraries. ####
