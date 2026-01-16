import matplotlib.pyplot as plt

# font sizes (may move to top of the whatever code)
titleSize = 32 
labelSize = 25
legendSize = 18
axisSize = 15
colors6Same = ["#F7BD2B", "#C78F01", "#5A440B", "#42DBDB", "#019E9E", "#002929"]
colors2 = ["#FFC20A", "#0C7BDC"] # davidmathlogic.com/colorblind/
colors6 = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

fig, ax = plt.subplots(figsize=(20, 8))

# titles / axis labels
ax.set_title(f'Graph Title',fontsize=titleSize)
ax.set_xlabel(f'X Axis Label', fontsize=labelSize)
ax.set_ylabel(f'Y Axis Label', fontsize=labelSize)

# plot data
xAxisList = [1,2,3]
yAxisList = [2,4,6]

# creates a line. https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html 
ax.plot(xAxisList, yAxisList, label=f"Data Label", color=colors6[0], linestyle="-", marker=".")
# places dots. https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
ax.scatter(xAxisList, yAxisList, label=f"Data Label", color=colors6[1], marker=".") 


# set x-ticks
xTicks = range(0,len(xAxisList)+1,1)
ax.set_xticks(ticks=xTicks)

# set y-ticks
yTicks = range(0,max(yAxisList)+1,1)
ax.set_yticks(ticks=yTicks)
ax.tick_params(direction='out', length=6, width=2, colors='black', grid_color='black', grid_alpha=1, labelsize=axisSize)

# x / y limits
ax.set_ylim(bottom=0, top=7)
ax.set_xlim(left=0, right=4)


# legend
ax.legend(loc='upper right', shadow=True, ncol=4, fontsize=legendSize)

# fun annotations (may want to make a seperate function for special annotations)
ax.axvspan(1, 2, color='#7CAED1', alpha=0.3) # highlights a section
ax.annotate('Annotation',xy=(1,1),horizontalalignment='left', verticalalignment='top',fontsize=axisSize) # puts text in graph

plt.style.context("seaborn-talk")
plt.savefig(f'file_name.png', format='png', bbox_inches='tight')
