import re
from collections import defaultdict
import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
import pickle




# Set Arial globally in Matplotlib configuration
plt.rcParams['font.family'] = 'DejaVu Sans'

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = 9
plt.rcParams['svg.fonttype'] = 'none' 
plt.rcParams['lines.linewidth'] = 0.5


color_dict = {'GEARS': '#FFCC99',
              'Geneformer': '#99CC99',
              'scGPT': '#FF9999',
              'scFoundation': '#D1B26F',
              'UCE': '#99FFCC',
              'scLong': '#CC99FF'}


with open('pkls/norman_mse_pod_pearson.pkl', 'rb') as file:
    value_dict = pickle.load(file)

z = np.random.normal()

# Plotting     
models = ['GEARS', 'Geneformer', 'scGPT', 'scFoundation',  'UCE','scLong']

fig, axs = plt.subplots(1,2)
metric_to_plot = list(value_dict.keys())
indices = np.arange(4) * 1.3
fig.set_size_inches(7.2, 1.8)
scenario_names = ['Seen 0/2', 'Seen 1/2', 'Seen 2/2', 'Seen 0/1']
for i, metric in enumerate(metric_to_plot):
    ax = axs[i]
    width = 1.2 / (len(models) + 2)
    for model_idx, model in enumerate(models):
        values = np.array(value_dict[metric][model])[[0, 2, 4, 6]]
        errors = np.array(value_dict[metric][model])[[1, 3, 5, 7]]
        samples = np.array(value_dict[metric][model])[8:].reshape((4, 5))
        label = model
        # draw bars
        bar_positions = indices + model_idx * width
        ax.bar(bar_positions, values, yerr=errors, label=label, width=width, color=color_dict[model])
        for k, pos in enumerate(bar_positions):
            y = samples[k]
            x = pos + np.linspace(-0.25, 0.25, len(y)) * width
            ax.scatter(x, y, color='black', s=1, zorder=3,marker='.')
    if i == 0:
        ax.set_ylim(0.1,0.32)
    if i == 1:
        ax.set_ylim(0.5,0.92)

    #ax.set_xlabel('Scenario', fontsize = 11.5)
    ax.set_ylabel(f'{metric}')
    if i == 0:
        ax.set_yticks(ticks = np.arange(2, 7) * 0.05)
        ax.legend(bbox_to_anchor=(1, 1.1), loc='center', ncol=6)
        ax.get_legend().draw_frame(False)
    if i == 1:
        ax.set_yticks(ticks = np.arange(5, 10) * 0.1)
    ax.set_xticks(ticks=indices + width * (len(models) - 1) /2)
    ax.set_xticklabels(scenario_names)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

#fig.subplots_adjust(hspace=0, wspace=0.5)
fig.subplots_adjust(hspace=0, wspace=0.3, bottom=0.15)

#fig.tight_layout()

fig.savefig(f'figs/norman_split_1_mse_pearson.svg', format='svg')
fig = None
plt.close()
