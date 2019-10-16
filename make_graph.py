import matplotlib.pyplot as py
import numpy as np
import pandas as pd
import seaborn as sns


# constants
dataset = "performance2.csv"
output = "output2.png"

# import dataset
data = pd.read_csv(dataset)

# set seaborn style
sns.set()

# plot the thingy
plot = sns.boxplot(data=data, x="operation", y="time", hue="treatment")
plot.set(ylim=(0, 400))

# save to disk
fig = plot.get_figure()
fig.savefig(output)