import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import importlib
import matplotlib.patches as patches


COLORL_grey = "#F0F0F0"
COLORL_g = "#B2D3C5"
COLORL_b = "#B2ECFF"
COLORL_r = "#F8BBC4"
COLORL_y = "#FFFCB2"
COLOR_grey = "#F0F0F0"
COLOR_g = "#5EA285"
COLOR_b = "#3D9BCF"#"#13C7FF"
COLOR_r = "#CE545D"#E8223E"
COLOR_y = "#FFFCB2"
dir = "/mnt/c/Users/myrth/Desktop/coding/raptor/build/bin/evaluation/"

# Adjust line and tick mark sizes
mpl.rcParams['lines.linewidth'] = 0.8
mpl.rcParams['lines.markersize'] = 3.0
myarrow = patches.ArrowStyle("Fancy",  head_length=0.05,head_width=0.05)


# Figure and axes
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

# First plot: Insertions
ax1.set_ylabel('Time (s)', color='black')
ax1.tick_params(axis='y', colors='black')

ax1_twin = ax1.twinx()
ax1_twin.set_ylabel('Memory (MB)', color='black')
ax1_twin.tick_params(axis='y', colors='black')

# Second plot: Queries
ax2.set_ylabel('Time (s)', color='black')
ax2.tick_params(axis='y', colors='black')
ax2_twin = ax2.twinx()
ax2_twin.set_ylabel('Memory (MB)', color='black')
ax2_twin.tick_params(axis='y', colors='black')

# Third plot: Index size
ax3_twin = ax3.twinx()
ax3_twin.set_ylabel('Memory (MB)')
ax3_twin.set_xlabel('UB insertions')

# plotting data
for method, colour in {"naive":'black', "find_ibf_idx_traverse_by_similarity": COLOR_b, "find_ibf_idx_ibf_size": COLOR_g, "find_ibf_idx_traverse_by_fpr": COLOR_r}.items() :
    mod = importlib.import_module(  method )

    x = np.arange(len(mod.time_query))
    ax1.plot(x[:-1] + 0.5, mod.time_insertion, color=colour, marker='o', linestyle='-', label='Time; ' + method) # Shifted x-values for the second and third plots
    #ax1_twin.plot(x[:-1] + 0.5,  [t/1000000 for t in mod.memory_insertion], color='None', markeredgecolor=colour, marker='o', linestyle='--', label='Memory ' + method)
    ax2.plot(x, mod.time_query, color=colour, marker='o', linestyle='-', label='Time')
    ax2_twin.plot(x, [t/1000000 for t in mod.memory_query], color='None', markeredgecolor=colour, marker='o', linestyle='--', label='Memory')
    ax3_twin.plot(x, [t/1000000 for t in mod.size_index], color='None', markeredgecolor=colour, marker='o', linestyle='--', label='Size')


    # Annotate a point with partial rebuild
    if method == "find_ibf_idx_traverse_by_fpr":
        first_label = False
        for i in range(len(mod.rebuild)):
            if mod.rebuild[i] == 1:
                if first_label == False:
                    ax1.annotate('partial rebuild', xy=(i + 0.5, mod.time_insertion[i]), xytext=(0, 20), textcoords='offset points', ha='center',
                                 arrowprops=dict(arrowstyle='->'))
                    first_label = True
                else:
                    ax1.annotate('', xy=(i + 0.5, mod.time_insertion[i]), xytext=(0, 20), textcoords='offset points', ha='center',
                                 arrowprops=dict(arrowstyle='->'))
    # ax3_twin.annotate('partial rebuild', xy=(label_rebuilds[0] + 0.5, size_index[label_rebuilds[0]]), xytext=(0, 20), textcoords='offset points', ha='center',
    #                   arrowprops=dict(arrowstyle='-')#, connectionstyle='arc3,rad=0.5')
    #                   )

# Set subplot and y-labels.
#ax1.annotate(' ', xy=(-5, 20000), xytext=(-5,20000),arrowprops=dict(arrowstyle='-', ls="--"), va='center', rotation=90,  annotation_clip=False)

# Set subplot and y-labels.
ax1.yaxis.set_label_coords(-0.15, 0.5)
ax1.set_ylabel('Insertions \n\n Time (s)', va='bottom', ha='center') #, ha='right')
ax2.yaxis.set_label_coords(-0.15, 0.5)
ax2.set_ylabel('Queries \n\n Time (s)', va='bottom', ha='center')
ax3.yaxis.set_label_coords(-0.15, 0.5)
ax3.set_ylabel('Index Size \n\n ', va='bottom', ha='center')

# change x axis 
ax3.set_xticklabels([])
ax1.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
ax2.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
ax3.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
#ax1.set_ylim(0, max(time_insertion)*1.5)
ax1.set_xlim(0, len(mod.time_insertion))
ax2.set_ylim(0)

ax3.set_yticklabels([])

# legend
ax1.legend(edgecolor="white", bbox_to_anchor=(0.2, 1), loc='upper left',)

# Adjust layout
plt.tight_layout()

# Display the plots
plt.show("results/mixed_bins.pdf")
plt.show()


