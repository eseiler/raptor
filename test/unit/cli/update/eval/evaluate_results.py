import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import timeit
import math
import random
# from exhaustive import *
# from simulate_data import *
import time
import importlib

COLOR_exhaustive = "grey"
COLOR_h1 = "xkcd:azure"
COLOR_h2 = "teal"
COLOR_h3 = "green"

"""
EVALUATE C++ HEURISTICS
"""

def evaluate_cpp_probability_distribution(dir):
    i = importlib.import_module(dir + "_distribution")
    i.results_exhaustive
    bins = np.histogram(np.hstack(i.results_exhaustive), bins=20)[1]

    fig, axs = plt.subplots(3, 1, sharex=True)#, subplot_kw=dict(projection="polar"))


    ax_i=0
    for result, name, colour in (
            (i.results_h1, "heuristic 1, $k$=" + str(i.k)+ ", fraction found: " +
                           str(round(len(i.results_h1)/len(i.results_exhaustive)*100))+"%", COLOR_h1 ),
            (i.results_h2, "heuristic 2, $b$=" + str(i.b)+ ", fraction found: " +
                           str(round(len(i.results_h2)/len(i.results_exhaustive)*100))+"%", COLOR_h2),
            (i.results_h3, "heuristic 3, $k$=" + str(i.k) + ", fraction found: " +
                           str(round(len(i.results_h3)/len(i.results_exhaustive)*100))+"%" , COLOR_h3)
    ):
        if ax_i == 0: label = "exact, |{Strings}|=" + str(len(i.results_exhaustive))
        else: label = None
        axs[ax_i].hist(i.results_exhaustive, bins=bins, color="lightgrey", label=label)
        axs[ax_i].hist(result, label=name, bins=bins, alpha=1, color=colour)
        axs[ax_i].set_xlim(xmin=bins[0], xmax=max(i.results_exhaustive))
        axs[ax_i].set_ylabel("|{Strings}|")
        axs[ax_i].legend(edgecolor="white", bbox_to_anchor=(0.2, 1), loc='upper left',)
        ax_i+=1
    plt.xlabel("Log$_2$ probability")


    # plt.annotate('Most probable string', xy =(list(result_exhaustive.values())[0], 1),
    #             xytext =(list(result_exhaustive.values())[0], 1 + 1),
    #             arrowprops = dict(facecolor ='black', width=0.1),)
    # if reasonable, start at x=0
    #plt.xlim(xmax=0)
    #plt.axvline(x=max(i.result_exhaustive), color="k", linewidth=0.5)
    #plt.axvline(x=bins[0], color="k", linewidth=0.5)
    plt.savefig(dir.replace(".", "\\") + "_dist.pdf")
    fig.suptitle("Profile: " + i.file_name +
                 "\n$n$: " + str(i.n) +
                 ", $T$: " + str(i.T),
                 fontdict = {'fontsize': 11,'horizontalalignment': 'left'})
    plt.savefig(dir.replace(".", "\\") + "_dist.png")
    plt.show()

# N vs. TOP FRACTION VS SPEEDUP
def evaluate_cpp_vs_n(dir):
    i=importlib.import_module(dir + "_vsn")
    fig, ax = plt.subplots()
    ax2 = plt.twinx()

    for y1, y2, label, colour, annotations in (
            (i.found_fractions_h1, i.speedups_h1, "heuristic 1", COLOR_h1 ,["$k$="+str(round(k)) for k in i.ks]),
            (i.found_fractions_h2, i.speedups_h2, "heuristic 2, b=" + str(round(i.bs[0],0)), COLOR_h2, ["" for b in i.bs]),
            (i.found_fractions_h3, i.speedups_h3, "heuristic 3", COLOR_h3, ["$k$="+str(round(k)) for k in i.ks])):
        y1 = [item * 100 for item in y1]
        ax.plot(i.ns, y1, label=label, marker='.', color=colour)
        ax2.plot(i.ns, y2, label=label, marker='.', color=colour, ls="--")
        for label, n, y in zip(annotations, i.ns, y2):
            ax2.annotate(label, xy=(n, y), xytext=(n, y ))
    ax.annotate("$T$", xy=(0, 0), xytext=(10, 40), xycoords='figure pixels')
    ax.annotate("$t_{exact}$ (ms)", xy=(-1, 0), xytext=(10, 25), xycoords='figure pixels')
    ax.annotate("$|Strings|$", xy=(-1, 0), xytext=(10, 10), xycoords='figure pixels')
    #normalized_thresholds = np.array([math.log2(label) for label in i.thresholds])*np.array(i.ns)
    for label, n, y in zip(i.thresholds, i.ns, y1):
        ax.annotate(str(round(label)), xy=(n, 0), xytext=(n, -14))
    for label, n, y in zip(i.times_exhaustive, i.ns, y1):
        ax.annotate(str(round(label*100,4)), xy=(n, 0), xytext=(n , -20))
    for label, n, y in zip(i.sizes_exhaustive, i.ns, y1):
        ax.annotate( str(round((label))), xy=(n, 0), xytext=(n, -26))

    ax.set_ylim(ymin=0, ymax=100)
    ax.set_xlim(xmin=4, xmax=21)
    ax2.set_ylim(ymin=0)
    ax2.axhline(y=1, color="black", linewidth = 0.5 )
    ax.set_xlabel("$n$")
    ax2.set_ylabel("Speedup")
    ax.set_ylabel("Top fraction (%)")
    ax.legend(edgecolor="white")
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    plt.savefig(dir.replace(".", "\\") + "_vsn.pdf")
    plt.title("Profile: " + i.file_name,
              fontdict = {'fontsize': 11,'horizontalalignment': 'left'})
    plt.show()

# TOP FRACTION VS SPEEDUP
def evaluate_cpp_top_vs_speed(dir):
    i=importlib.import_module(dir + "_topvsspeed")
    for x, y, annotations, label, color in (
            (i.found_fractions_h1, i.speedups_h1, ["$k$="+str(k) for k in i.ks], "heuristic 1", COLOR_h1),
            (i.found_fractions_h2, i.speedups_h2, ["$b$="+str(b) for b in i.bs], "heuristic 2",COLOR_h2),
            (i.found_fractions_h3, i.speedups_h3, ["$k$="+str(k) for k in i.ks], "heuristic 3",COLOR_h3)
    ):
        x = [item*100 for item in x]
        plt.plot(x, y, label=label, marker='.', color=color)  # , color=color)
        for label, x_val, y_val in zip(annotations, x, y):
            plt.annotate(label, xy=(x_val, y_val), xytext=(x_val, y_val ))
    plt.ylim(ymin=0) #plt.ylim(ymin=0)
    plt.xlim(xmin=0, xmax=100)
    plt.axhline(y=1, color="black", linewidth = 0.75 )
    plt.annotate("$t_{exact}$=" + str(round(i.time_exhaustive*100, 4)) + " ms", xy=(1,1))
    plt.ylabel("Speedup ")
    plt.xlabel("Top strings found (%)")
    plt.legend(edgecolor="white")#,loc=4)
    plt.savefig(dir.replace(".", "\\") + "_topvsspeed.pdf")
    plt.title("Profile: " + i.file_name +
              "\n$n$: " + str(i.n) +
              ", $T$: " + str(i.T) +
              "\nTotal number of found strings: " + str(i.size_exhaustive),
              fontdict = {'fontsize': 11,'horizontalalignment': 'left'})
    plt.savefig(dir.replace(".", "\\") + "_topvsspeed_optimization.png")
    plt.show()

def evaluate_cpp_top_vs_speed_objv(dir):
    i=importlib.import_module(dir + "_topvsspeed")
    for x, y, annotations, label, color in (
            (i.found_fractions_h1, i.speedups_h1, ["$k$="+str(k) for k in i.ks], "heuristic 1", COLOR_h1),
            (i.found_fractions_h2, i.speedups_h2, ["$b$="+str(b) for b in i.bs], "heuristic 2",COLOR_h2),
            #(i.found_fractions_h3, i.speedups_h3, ["$k$="+str(k) for k in i.ks], "heuristic 3",COLOR_h3)
    ):
        x = [item*100 for item in x]
        y = -((i.time_exhaustive/np.array(y))-i.time_exhaustive)/i.time_exhaustive
        plt.plot(x, y, label=label, marker='.', color=color)  # , color=color)
        for label, x_val, y_val in zip(annotations, x, y):
            plt.annotate(label, xy=(x_val, y_val), xytext=(x_val, y_val ))
            # , xytext=(2, 2), arrowprops=dict(facecolor='black', shrink=0.05))
    plt.plot([50,100], [1,0], color='black', ls='--', linewidth=0.75, label="Objective function")
    plt.ylim(ymin=-0.10, ymax=1) #plt.ylim(ymin=0)
    plt.xlim(xmin=0, xmax=100)
    plt.axhline(y=0, color="black", linewidth = 0.75 )
    plt.annotate("$t_{exact}$=" + str(round(i.time_exhaustive*100, 4)) + " ms", xy=(0,0.005))
    plt.ylabel("-$\Delta t$/$t_{exact}$")
    plt.xlabel("Top strings found (%)")
    #plt.plot(x_exh, y_exh, label=label_exh, marker='.', color=color_exh)
    plt.legend(edgecolor="white")#,loc=4)
    plt.savefig(dir.replace(".", "\\") + "_topvsspeed_objv.pdf")
    plt.title("Profile: " + i.file_name +
              "\n$n$: " + str(i.n) +
              ", $T$: " + str(i.T) +
              "\nTotal number of found strings: " + str(i.size_exhaustive),
              fontdict = {'fontsize': 11,'horizontalalignment': 'left'})
    plt.savefig(dir.replace(".", "\\") + "_topvsspeed_optimization.png")
    plt.show()

#----------------------------------
folder = "data.DNA_profiles.results."
file_names = ("MA0007", "MA0528")#, "MA0659")
# folder = "data.RNA_profiles.results."
# file_names = ("RF0000", "RF0001", "RF0052")
# folder = "data.protein_profiles.results."
# file_names = ("PS5080", "PS5101", "PS5180")

for file_name in file_names:
    evaluate_cpp_vs_n(folder + file_name)
    # evaluate_cpp_probability_distribution(folder + file_name)
    # evaluate_cpp_top_vs_speed(folder + file_name)
    # evaluate_cpp_top_vs_speed_objv(folder + file_name)