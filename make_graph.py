import pandas as pd
import seaborn as sns

# constants
dataset = "performance.csv"
output = "output.png"

# import dataset
data = pd.read_csv(dataset)

# define labels
hue_order = ["bit_vector", "bit_set"]
x_order = ["emp::vector<T>", "T"]
graph_titles = ['CountOnes_Mixed', 'NOT', 'AND', 'OR']
legend_labels = ["BitVector", "BitSet"]

# define pretty colors
colors = ['#4a73ab', '#8a9b26']
bg_color = '#ede4d5'

# set graph style
sns.set_style("dark")
sns.set(rc={"axes.facecolor" : bg_color, "figure.facecolor" : bg_color,
            "savefig.facecolor": "#fefefe", "grid.color": "#fefefe"},
        font="Arial")

# generate palette
palette = {k : v for k, v in zip(hue_order, colors)}

# plot the thingy
plot = sns.FacetGrid(data, col="operation", col_wrap=2)
plot = plot.map(sns.barplot, "container", "time", "treatment", 
                order=x_order, hue_order=hue_order, log=True, dodge=True,
                palette=palette, edgecolor=bg_color)
plot.add_legend()

# fix graphs labels
graphs = plot.axes.flatten()

# fix y-axis labels
graphs[0].set_ylabel("time (ns)")
graphs[2].set_ylabel("time (ns)")

# remove x-axis labels
graphs[2].set_xlabel('')
graphs[3].set_xlabel('')

# fix legend labels
for t, l in zip(plot._legend.texts, legend_labels): t.set_text(l)

# set graph titles
for i in range(len(graphs)):
        graphs[i].set_title(graph_titles[i])

# save to disk
plot.savefig(output)