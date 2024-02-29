#%%
import seaborn as sns
from matplotlib.pyplot import gcf
import matplotlib.pyplot as plt
import pandas as pd
# fig, axs = plt.subplots()

networks = sns.load_dataset("brain_networks", index_col=0, header=[0, 1, 2])

# Label 1
network_labels = networks.columns.get_level_values("network")
network_pal = sns.cubehelix_palette(network_labels.unique().size, light=.9, dark=.1, reverse=True, start=1, rot=-2)
network_lut = dict(zip(map(str, network_labels.unique()), network_pal))
network_colors = pd.Series(network_labels, index=networks.columns).map(network_lut)

# Label 2
node_labels = networks.columns.get_level_values("node")
node_pal = sns.cubehelix_palette(node_labels.unique().size)
node_lut = dict(zip(map(str, node_labels.unique()), node_pal))
node_colors = pd.Series(node_labels, index=networks.columns).map(node_lut)

# Label 3
lab3_labels = networks.columns.get_level_values("node")
lab3_pal = sns.color_palette("hls", lab3_labels.unique().size)
lab3_lut = dict(zip(map(str, lab3_labels.unique()), lab3_pal))
lab3_colors = pd.Series(lab3_labels, index=networks.columns, name='lab3').map(lab3_lut)

# Label 4
lab4_labels = networks.columns.get_level_values("node")
lab4_pal = sns.color_palette("husl", lab4_labels.unique().size)
lab4_lut = dict(zip(map(str, lab4_labels.unique()), lab4_pal))
lab4_colors = pd.Series(lab4_labels, index=networks.columns, name='lab4').map(lab4_lut)

network_node_colors = pd.DataFrame(network_colors).join(pd.DataFrame(node_colors)).join(pd.DataFrame(lab3_colors)).join(pd.DataFrame(lab4_colors))

g = sns.clustermap(networks.corr(),
    row_cluster=False, col_cluster=False,
    row_colors = network_node_colors,
    col_colors = network_node_colors,
    linewidths=0,
    xticklabels=False, yticklabels=False,
    center=0, cmap="vlag")


# add legends
for label in network_labels.unique():
    g.ax_col_dendrogram.bar(0, 0, color=network_lut[label], label=label, linewidth=0);
l1 = g.ax_col_dendrogram.legend(title='Network', loc="center", ncol=5, bbox_to_anchor=(0.35, 0.89), bbox_transform=gcf().transFigure)

for label in node_labels.unique():
    g.ax_row_dendrogram.bar(0, 0, color=node_lut[label], label=label, linewidth=0);
l2 = g.ax_row_dendrogram.legend(title='Node', loc="center", ncol=2, bbox_to_anchor=(0.66, 0.89), bbox_transform=gcf().transFigure)

# create a list for the bar plot patches
xx = []
for label in lab3_labels.unique():
    x = g.ax_row_dendrogram.bar(0, 0, color=lab3_lut[label], label=label, linewidth=0)
    xx.append(x)

# add the legend
legend3 = plt.legend(xx, lab3_labels.unique(), loc="center", title='lab3', bbox_to_anchor=(.78, 0.89), bbox_transform=gcf().transFigure)


# create a list for the bar plot patches
yy = []
for label in lab4_labels.unique():
    y = g.ax_row_dendrogram.bar(0, 0, color=lab4_lut[label], label=label, linewidth=0)
    yy.append(y)

    
# add the second legend
legend4 = plt.legend(yy, lab4_labels.unique(), loc="center", title='lab4', ncol=2, bbox_to_anchor=(.9, 0.89), bbox_transform=gcf().transFigure)
plt.gca().add_artist(legend3)
#%%