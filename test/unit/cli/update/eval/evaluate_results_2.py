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
COLOR_grey = "#bfbfbf"
COLOR_g = "#5EA285"
COLOR_b = "#3D9BCF"#"#13C7FF"
COLOR_r = "#CE545D"#E8223E"
COLOR_y = "#FFFCB2"
dir = "/mnt/c/Users/myrth/Desktop/coding/raptor/build/bin/evaluation/"

# Adjust line and tick mark sizes
mpl.rcParams['lines.linewidth'] = 0.8
mpl.rcParams['lines.markersize'] = 3.0
myarrow = patches.ArrowStyle("Fancy",  head_length=0.05,head_width=0.05)

def ub_insertions():
    # Figure and axes
    fig, (  (ax3_twin, ax3), (ax1_twin, ax1 ), (ax2_twin, ax2)) = plt.subplots(3, 2, figsize=(10, 10), sharex=True)

    # First plot: Insertions
    ax1.set_ylabel('Time (s)', color='black')
    ax1.tick_params(axis='y', colors='black')
    ax1_right_axis = ax1.twinx()
    ax1_twin.set_ylabel('Memory (MB)', color='black')
    ax1_twin.tick_params(axis='y', colors='black')
    ax1_right_axis.set_ylabel('Cumulative time (min)')
    ax1.set_zorder(2)
    ax1.patch.set_visible(False)
    ax1_right_axis.set_zorder(1)

    # Second plot: Queries
    ax2.set_ylabel('Average time (s)', color='black')
    ax2.tick_params(axis='y', colors='black')
    ax2_twin.set_ylabel('Memory (MB)', color='black')
    ax2_twin.tick_params(axis='y', colors='black')

    # Third plot: Index size
    ax3_twin.set_ylabel('Memory (MB)')

    # plotting data

    for method, colour in {"mantis_results":"black", "naive":COLOR_grey, "find_ibf_idx_traverse_by_similarity": COLOR_b, "find_ibf_idx_ibf_size": COLOR_g, "find_ibf_idx_traverse_by_fpr": COLOR_r}.items() :
        mod = importlib.import_module(  method )
        if method == "naive": method = "Naive"
        if method == "find_ibf_idx_traverse_by_similarity": method = "Similarity"
        if method == "find_ibf_idx_ibf_size": method = "IBF size"
        if method == "find_ibf_idx_traverse_by_fpr": method = "FPR"
        if method == "mantis_results":
            mod.time_query=mod.time_query + [mod.time_query[-1]]
            mod.memory_query=mod.memory_query + [mod.memory_query[-1]]
            mod.size_index=np.array(mod.size_index + [0])*1000
            method = "Mantis"
        x = np.arange(len(mod.time_query))
        # Compute the cumulative sum of memory data
        cumulative_memory = np.cumsum(mod.time_insertion)/60
        ax1_right_axis.fill_between(x[:-1] + 0.5, 0, cumulative_memory, color=colour, alpha=0.2, zorder=0)


        ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0]) + 0.5,
                 [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 0], zorder=10, color=colour, marker='o',   linestyle='', label=method, alpha=1) # Shifted x-values for the second and third plots
        if method == "FPR":
            label = "Full rebuild"
        else:
            label = ""
        ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1]) + 0.5,
                 [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 1], zorder=10,   color=colour, marker='x', markersize=6, linestyle='', label=label, alpha=1) # Shifted x-values for the second and third plots
        ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2]) + 0.5,
                 [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 2], zorder=10,   color=colour, marker='+', markersize=6, linestyle='', alpha=1) # Shifted x-values for the second and third plots

        ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0]) + 0.5,
                 [mod.memory_insertion[i]/1000000  for i in x[:-1] if mod.rebuild[i] == 0], zorder=10, color=colour, marker='o',   linestyle='',  alpha=1) # Shifted x-values for the second and third plots
        ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1]) + 0.5,
                 [mod.memory_insertion[i]/1000000  for i in x[:-1] if mod.rebuild[i] == 1], zorder=10,   color=colour, marker='x', markersize=6, linestyle='', alpha=1) # Shifted x-values for the second and third plots
        ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2]) + 0.5,
                 [mod.memory_insertion[i]/1000000  for i in x[:-1] if mod.rebuild[i] == 2], zorder=10,   color=colour, marker='+', markersize=6, linestyle='', alpha=1) # Shifted x-values for the second and third plots


        #ax1_twin.plot(x[:-1] + 0.5,  [t/1000000 for t in mod.memory_insertion],  color=colour, marker='o', linestyle=''  )
        ax2.plot(x, mod.time_query, color=colour, marker='o', linestyle='')
        ax2_twin.plot(x, [t/1000000 for t in mod.memory_query], color=colour, marker='o', linestyle= '' )
        print(len(x), len(mod.size_index))
        ax3_twin.plot(x[:len(mod.size_index)], [t/1000000 for t in mod.size_index], color=colour, marker='o', linestyle='' )



        # Annotate a point with partial rebuild
        # if method == "find_ibf_idx_traverse_by_fpr":
        #     first_label = False
        #     for i in range(len(mod.rebuild)):
        #         if mod.rebuild[i] == 1:
        #             if first_label == False:
        #                 ax1.annotate('full rebuild', xy=(i + 0.5, mod.time_insertion[i]), xytext=(0, 20), textcoords='offset points', ha='center',
        #                              arrowprops=dict(arrowstyle='->')) #Partial rebuilds (another shape); partial rebuild, triangle or empty star, full rebuild full star.
        #                 first_label = True
        #             else:
        #                 ax1.annotate('', xy=(i + 0.5, mod.time_insertion[i]), xytext=(0, 20), textcoords='offset points', ha='center',
        #                              arrowprops=dict(arrowstyle='->'))
        #         elif mod.rebuild[i] == 2:
        #             print("x")
        #             ax1.annotate('', xy=(i + 0.5, mod.time_insertion[i]), xytext=(0, 50), textcoords='offset points', ha='center',
        #                          arrowprops=dict(arrowstyle='->'))
        # ax3_twin.annotate('partial rebuild', xy=(label_rebuilds[0] + 0.5, size_index[label_rebuilds[0]]), xytext=(0, 20), textcoords='offset points', ha='center',
        #                   arrowprops=dict(arrowstyle='-')#, connectionstyle='arc3,rad=0.5')
        #                   )

    # Set subplot and y-labels.
    #ax1.annotate(' ', xy=(-5, 20000), xytext=(-5,20000),arrowprops=dict(arrowstyle='-', ls="--"), va='center', rotation=90,  annotation_clip=False)

    # Set subplot and y-labels.
    ax1_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax1_twin.set_ylabel('Insertions \n\n Memory (MB)', va='bottom', ha='center') #, ha='right')
    ax2_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax2_twin.set_ylabel('Queries \n\n Memory (MB)', va='bottom', ha='center')
    ax3_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax3_twin.set_ylabel('Index Size \n\n Memory (MB)', va='bottom', ha='center')

    # change x axis
    ax2_twin.set_xlabel('Number of inserted samples')
    ticks = ax3_twin.get_xticks()
    tick_labels = ax3_twin.get_xticklabels()
    ax2.set_xlabel('Number of inserted samples')
    ax2.set_xticks(ticks, tick_labels)
    ax1.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
    ax2.tick_params(axis='x', bottom=True)  # Hide x-axis tick labels
    ax1_twin.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
    ax3_twin.tick_params(axis='x', bottom=False)

    # Set limits
    ax1.set_xlim(0, len(mod.time_query))#, len(mod.time_insertion))
    ax2.set_ylim(0)
    ax2.set_ylim(0)
    ax1.set_ylim(0)
    ax3_twin.set_ylim(0)
    ax1_twin.set_xlim(0, len(mod.time_query))#, len(mod.time_insertion))
    ax2_twin.set_xlim(0, len(mod.time_query))
    ax3_twin.set_xlim(0, len(mod.time_query))
    ax1_twin.set_ylim(0)
    ax2_twin.set_ylim(0)
    ax1_right_axis.set_ylim(0)


    # delete one of the subplots.
    fig.delaxes(ax3)

    # legend
    #ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
    fig.legend(loc='upper center', bbox_to_anchor=(0.70, 0.97), ncol=2, frameon=False)
    #fig.legend.get_frame().set_facecolor('none')


    #fig.text(0.75, 0, 'Time', ha='center')
    fig.text(0.28, 0.985, 'Memory', ha='center')
    ax1.set_title('Time', fontsize =10) # fontsize =3

    # Adjust layout
    plt.tight_layout()

    # Display the plots
    plt.savefig("results/mixed_bins.pdf")
    plt.show()

#################################
# SEQUENCE INSERTION

def seq_insertions():
    # Figure and axes
    fig, (  (ax3_twin, ax3), (ax1_twin, ax1 ), (ax2_twin, ax2)) = plt.subplots(3, 2, figsize=(10, 10), sharex=True)

    # First plot: Insertions
    ax1.set_ylabel('Time (s)', color='black')
    ax1.tick_params(axis='y', colors='black')
    ax1_right_axis = ax1.twinx()
    ax1_twin.set_ylabel('Memory (MB)', color='black')
    ax1_twin.tick_params(axis='y', colors='black')
    ax1_right_axis.set_ylabel('Cumulative time (min)')
    ax1.set_zorder(2)
    ax1.patch.set_visible(False)
    ax1_right_axis.set_zorder(1)

    # Second plot: Queries
    ax2.set_ylabel('Average time (s)', color='black')
    ax2.tick_params(axis='y', colors='black')
    ax2_twin.set_ylabel('Memory (MB)', color='black')
    ax2_twin.tick_params(axis='y', colors='black')

    # Third plot: Index size
    ax3_twin.set_ylabel('Memory (MB)')

    # plotting data
    for method, colour in {"results_insertsequences_with_FPRbuffer": COLOR_b, }.items() :
        mod = importlib.import_module(  "results_insertsequences_with_FPRbuffer" )
        x = np.arange(len(mod.time_query))

        index_full_rebuild = 0
        ax1_twin.axhline(mod.memory_insertion[index_full_rebuild]/1000000 , color=COLOR_grey, label="Naive", lw=3)
        ax1.axhline(mod.time_insertion[index_full_rebuild], color=COLOR_grey, lw=3)
        cumulative_memory_del = np.cumsum([mod.time_insertion[index_full_rebuild] for i in x[:-1]])/60
        ax1_right_axis.fill_between(x[:-1] + 0.5, 0, cumulative_memory_del, color=COLOR_grey, alpha=0.2, zorder=0)

    # Compute the cumulative sum of memory data
        cumulative_memory = np.cumsum(mod.time_insertion)/60
        ax1_right_axis.fill_between(x[:-1] + 0.5, 0, cumulative_memory, color=colour, alpha=0.2, zorder=0)

        ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0]) + 0.5,
                 [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 0], zorder=10, color=colour, marker='o',   linestyle='', label="Dynamic HIBF", alpha=1) # Shifted x-values for the second and third plots
        ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1]) + 0.5,
                 [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 1], zorder=10,   color=colour, marker='x', markersize=6, linestyle='', label="Full rebuild", alpha=1) # Shifted x-values for the second and third plots
        ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2]) + 0.5,
                 [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 2], zorder=10,   color=colour, marker='+', markersize=6, linestyle='',  label="Partial rebuild", alpha=1) # Shifted x-values for the second and third plots
        ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0]) + 0.5,
                 [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 0], zorder=10, color=colour, marker='o',   linestyle='',  alpha=1) # Shifted x-values for the second and third plots
        ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1]) + 0.5,
                 [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 1], zorder=10,   color=colour, marker='x', markersize=6, linestyle='',  alpha=1) # Shifted x-values for the second and third plots
        ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2]) + 0.5,
                 [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 2], zorder=10,   color=colour, marker='+', markersize=6, linestyle='',   alpha=1) # Shifted x-values for the second and third plots

        ax2.plot(x, mod.time_query, color=colour, marker='o', linestyle='')
        ax2_twin.plot(x, [t/1000000 for t in mod.memory_query], color=colour, marker='o', linestyle= '' )
        ax3_twin.plot(x, [t/1000000 for t in mod.size_index], color=colour, marker='o', linestyle='' )


    # Set subplot and y-labels.
    ax1_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax1_twin.set_ylabel('Insertions \n\n Memory (MB)', va='bottom', ha='center') #, ha='right')
    ax2_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax2_twin.set_ylabel('Queries \n\n Memory (MB)', va='bottom', ha='center')
    ax3_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax3_twin.set_ylabel('Index Size \n\n Memory (MB)', va='bottom', ha='center')

    # change x axis
    ax2_twin.set_xlabel('Batches of inserted sequences') # TODO perhaps compare with the time of 1 rebuild of the complete index, which I can easily get from the data itself.
    ticks = ax3_twin.get_xticks()
    tick_labels = ax3_twin.get_xticklabels()
    ax2.set_xlabel('Batches of inserted sequences')
    ax2.set_xticks(ticks, tick_labels)
    ax1.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
    ax2.tick_params(axis='x', bottom=True)  # Hide x-axis tick labels
    ax1_twin.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
    ax3_twin.tick_params(axis='x', bottom=False)

    # Set limits
    ax1.set_xlim(0, len(mod.time_query))#, len(mod.time_insertion))
    ax2.set_ylim(0)
    ax2.set_ylim(0)
    ax1.set_ylim(0)
    ax1_twin.set_xlim(0, len(mod.time_query))#, len(mod.time_insertion))
    ax2_twin.set_xlim(0, len(mod.time_query))
    ax3_twin.set_xlim(0, len(mod.time_query))
    ax1_twin.set_ylim(0)
    ax2_twin.set_ylim(0)
    ax1_right_axis.set_ylim(0)

    # delete one of the subplots.
    fig.delaxes(ax3)

    # legend
    fig.legend(loc='upper center', bbox_to_anchor=(0.70, 0.97), ncol=2, frameon=False)

    fig.text(0.28, 0.985, 'Memory', ha='center')
    ax1.set_title('Time', fontsize =10) # fontsize =3

    # Adjust layout
    plt.tight_layout()

    # Display the plots
    plt.savefig("results/insertions.pdf")
    plt.show()

#################################
# DELETIONS
def deletions():
    del_or_insert =[]
    number_of_operations = 0
    while number_of_operations < 30:
        #if user_bin_filenames_counter + number_of_operations < len(user_bin_filenames):
        for _ in range(number_of_operations):
            del_or_insert.append(1)
        for _ in range(number_of_operations):
            del_or_insert.append(0)
        number_of_operations += 1


    # Figure and axes
    fig, (  (ax3_twin, ax3), (ax1_twin, ax1 ), (ax2_twin, ax2)) = plt.subplots(3, 2, figsize=(10, 10), sharex=True)

    # First plot: Insertions
    ax1.set_ylabel('Time (s)', color='black')
    ax1.tick_params(axis='y', colors='black')
    ax1_right_axis = ax1.twinx()
    ax1_twin.set_ylabel('Memory (MB)', color='black')
    ax1_twin.tick_params(axis='y', colors='black')
    ax1_right_axis.set_ylabel('Cumulative time (min)')
    ax1.set_zorder(2)
    ax1.patch.set_visible(False)
    ax1_right_axis.set_zorder(1)

    # Second plot: Queries
    ax2.set_ylabel('Average time (s)', color='black')
    ax2.tick_params(axis='y', colors='black')
    ax2_twin.set_ylabel('Memory (MB)', color='black')
    ax2_twin.tick_params(axis='y', colors='black')

    # Third plot: Index size
    ax3_twin.set_ylabel('Memory (MB)')

    # plotting data
    mod = importlib.import_module(  "deletions" )
    index_full_rebuild = -2
    ax1_twin.axhline(mod.memory_insertion[index_full_rebuild]/1000000 , color=COLOR_grey, label="Naive", lw=3)
    ax1.axhline(mod.time_insertion[index_full_rebuild], color=COLOR_grey, lw=3)

    mod = importlib.import_module(  "deletions" )
    x = np.arange(len(mod.time_query))

    del_or_insert = del_or_insert[:len(x)-1]

    colormap = np.array([COLOR_g, COLOR_r])
    # Compute the cumulative sum of memory data
    cumulative_memory_insert = np.cumsum([mod.time_insertion[i]  if del_or_insert[i] ==1 else 0 for i in x[:-1] ])/60
    ax1_right_axis.fill_between(x[:-1] + 0.5, 0, cumulative_memory_insert, color=colormap[1], alpha=0.2, zorder=0)
    cumulative_memory_del = np.cumsum([mod.time_insertion[i] if del_or_insert[i] ==0 else 0 for i in x[:-1]])/60
    ax1_right_axis.fill_between(x[:-1] + 0.5, 0, cumulative_memory_del, color=colormap[0], alpha=0.2, zorder=0)
    cumulative_memory_del = np.cumsum([mod.time_insertion[index_full_rebuild] if del_or_insert[i] ==0 else 0 for i in x[:-1]])/60
    ax1_right_axis.fill_between(x[:-1] + 0.5, 0, cumulative_memory_del, color=COLOR_grey, alpha=0.2, zorder=0)
    cumulative_memory_del = np.cumsum([mod.time_insertion[index_full_rebuild] for i in x[:-1]])/60
    ax1_right_axis.fill_between(x[:-1] + 0.5, 0, cumulative_memory_del, color=COLOR_grey, alpha=0.2, zorder=0)

    ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==0]) + 0.5,
             [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==0], zorder=10, c=colormap[0], marker='o',   linestyle='', label="Insertions", alpha=1) # Shifted x-values for the second and third plots
    ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==0]) + 0.5,
             [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==0], zorder=10,   color=colormap[0], marker='x', markersize=6, linestyle='', label="Full rebuild", alpha=1) # Shifted x-values for the second and third plots
    ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==0]) + 0.5,
             [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==0], zorder=10,   color=colormap[0], marker='+', markersize=6.5, linestyle='', label="Partial rebuild",  alpha=1) # Shifted x-values for the second and third plots

    ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==1]) + 0.5,
             [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==1], zorder=10, c=colormap[1], marker='o',   linestyle='', label="Deletions", alpha=1) # Shifted x-values for the second and third plots
    ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==1]) + 0.5,
             [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==1], zorder=10,   color=colormap[1], marker='x', markersize=6, linestyle='',alpha=1) # Shifted x-values for the second and third plots
    ax1.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==1]) + 0.5,
             [mod.time_insertion[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==1], zorder=10,   color=colormap[1], marker='+', markersize=6.5, linestyle='',  alpha=1) # Shifted x-values for the second and third plots

    ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==0]) + 0.5,
             [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==0], zorder=10, c=colormap[0], marker='o',   linestyle='',  alpha=1) # Shifted x-values for the second and third plots
    ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==0]) + 0.5,
             [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==0], zorder=10,   color=colormap[0], marker='x', markersize=6, linestyle='',  alpha=1) # Shifted x-values for the second and third plots
    ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==0]) + 0.5,
             [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==0], zorder=10,   color=colormap[0], marker='+', markersize=6.5, linestyle='',  alpha=1) # Shifted x-values for the second and third plots
    ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==1]) + 0.5,
             [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==1], zorder=10, c=colormap[1], marker='o',   linestyle='', alpha=1) # Shifted x-values for the second and third plots
    ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==1]) + 0.5,
             [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==1], zorder=10,   color=colormap[1], marker='x', markersize=6, linestyle='',alpha=1) # Shifted x-values for the second and third plots
    ax1_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==1]) + 0.5,
             [mod.memory_insertion[i]/1000000 for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==1], zorder=10,   color=colormap[1], marker='+', markersize=6.5, linestyle='',  alpha=1) # Shifted x-values for the second and third plots

    #ax1_twin.plot(x[:-1] + 0.5,  [t/1000000 for t in mod.memory_insertion][:-1],  color=COLOR_b, marker='o', linestyle=''  )
    ax2.plot(x, mod.time_query, color=COLOR_b, marker='o', linestyle='', label="Queries & index size")
    ax2_twin.plot(x, [t/1000000 for t in mod.memory_query], color=COLOR_b, marker='o', linestyle= '' )
    #ax3_twin.plot(x, [t/1000000 for t in mod.size_index], color=COLOR_b, marker='o', linestyle='' )

    ax3_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==0]) + 0.5,
                  [mod.size_index[i+1]/1000000 for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==0], zorder=10, c=colormap[0], marker='o',   linestyle='',  alpha=1) # Shifted x-values for the second and third plots
    ax3_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==0]) + 0.5,
                  [mod.size_index[i+1]/1000000 for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==0], zorder=10,   color=colormap[0], marker='x', markersize=6, linestyle='',  alpha=1) # Shifted x-values for the second and third plots
    ax3_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==0]) + 0.5,
                  [mod.size_index[i+1]/1000000 for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==0], zorder=10,   color=colormap[0], marker='+', markersize=6.5, linestyle='',  alpha=1) # Shifted x-values for the second and third plots
    ax3_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==1]) + 0.5,
                  [mod.size_index[i+1]/1000000 for i in x[:-1] if mod.rebuild[i] == 0 and del_or_insert[i]==1], zorder=10, c=colormap[1], marker='o',   linestyle='', alpha=1) # Shifted x-values for the second and third plots
    ax3_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==1]) + 0.5,
                  [mod.size_index[i+1]/1000000 for i in x[:-1] if mod.rebuild[i] == 1 and del_or_insert[i]==1], zorder=10,   color=colormap[1], marker='x', markersize=6, linestyle='',alpha=1) # Shifted x-values for the second and third plots
    ax3_twin.plot(np.array([x[i] for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==1]) + 0.5,
                  [mod.size_index[i+1]/1000000 for i in x[:-1] if mod.rebuild[i] == 2 and del_or_insert[i]==1], zorder=10,   color=colormap[1], marker='+', markersize=6.5, linestyle='',  alpha=1) # Shifted x-values for the second and third plots


    # Set subplot and y-labels.
    ax1_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax1_twin.set_ylabel('Insertions \n\n Memory (MB)', va='bottom', ha='center') #, ha='right')
    ax2_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax2_twin.set_ylabel('Queries \n\n Memory (MB)', va='bottom', ha='center')
    ax3_twin.yaxis.set_label_coords(-0.15, 0.5)
    ax3_twin.set_ylabel('Index Size \n\n Memory (MB)', va='bottom', ha='center')

    # change x axis
    ax2_twin.set_xlabel('Number of inserted/deleted samples')
    ticks = ax3_twin.get_xticks()
    tick_labels = ax3_twin.get_xticklabels()
    ax2.set_xlabel('Batches of inserted/deleted samples')
    ax2.set_xticks(ticks, tick_labels)
    ax1.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
    ax2.tick_params(axis='x', bottom=True)  # Hide x-axis tick labels
    ax1_twin.tick_params(axis='x', bottom=False)  # Hide x-axis tick labels
    ax3_twin.tick_params(axis='x', bottom=False)

    # Set limits
    ax1.set_xlim(0, len(mod.time_query))#, len(mod.time_insertion))
    ax2.set_ylim(0)
    ax2.set_ylim(0)
    ax1.set_ylim(0)
    ax1_twin.set_xlim(0, len(mod.time_query))#, len(mod.time_insertion))
    ax2_twin.set_xlim(0, len(mod.time_query))
    ax3_twin.set_xlim(0, len(mod.time_query))
    ax1_twin.set_ylim(0)
    ax2_twin.set_ylim(0)
    ax1_right_axis.set_ylim(0)

    # delete one of the subplots.
    fig.delaxes(ax3)

    # legend
    fig.legend(loc='upper center', bbox_to_anchor=(0.70, 0.97), ncol=2, frameon=False)

    fig.text(0.28, 0.985, 'Memory', ha='center')
    ax1.set_title('Time', fontsize =10) # fontsize =3

    # Adjust layout
    plt.tight_layout()

    # Display the plots
    plt.savefig("results/insertions.pdf")
    plt.show()

#################################
# BENCHMARK AGAINST MANTIS
# for method, colour in {"mantis":COLOR_grey, "find_ibf_idx_traverse_by_fpr": COLOR_r}.items() :
#     mod = importlib.import_module(  method )

ub_insertions()