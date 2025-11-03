import re
from collections import defaultdict
import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt




# Set Arial globally in Matplotlib configuration
plt.rcParams['font.family'] = 'DejaVu Sans'

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = 5
plt.rcParams['svg.fonttype'] = 'none' 
plt.rcParams['lines.linewidth'] = 0.5


color_dict = {'w/o LEG': '#99CC99',
              'Random  LEG': '#FF9999',
              'scLong': '#CC99FF'}
import pickle
with open('../pkls/abl_low_dict.pkl', 'rb') as file:
    abl_low_dict = pickle.load(file)

value_dict = abl_low_dict['value_dict']
data1 = abl_low_dict['data1']
data2 = abl_low_dict['data2']
error1 = abl_low_dict['error1']
error2 = abl_low_dict['error2']
data3 = abl_low_dict['data3']


# Plotting     
models = ['w/o LEG', 'Random  LEG', 'scLong']

colors = [color_dict[model] for model in models]


fig = plt.figure()  # 设置整体图形大小
fig.set_size_inches(7.5, 1.8)
gs = fig.add_gridspec(1, 21)
axs = [fig.add_subplot(gs[0, 0:5]),
       fig.add_subplot(gs[0, 6:11]),
       fig.add_subplot(gs[0, 12:14]),
       fig.add_subplot(gs[0, 15:17]),
       fig.add_subplot(gs[0, 19:21])]

metric_to_plot = list(value_dict.keys())
indices = np.arange(4)
fig.set_size_inches(8.5, 2.2)
scenario_names = ['Seen 0/2', 'Seen 1/2', 'Seen 2/2', 'Seen 0/1']
for j in range(5):
    ax = axs[j]
    width = 1 / (len(models)+1)
    if j in {0,1}:
        i = j
        metric = metric_to_plot[i]
        for model_idx, model in enumerate(models):
            values = np.array(value_dict[metric][model])[[0, 2, 4, 6]]
            errors = np.array(value_dict[metric][model])[[1, 3, 5, 7]]
            label = model
            ax.bar(indices + model_idx * width, values, yerr=errors, label=label, width=width, color= color_dict[model]) #, color=colors[model_idx % len(colors)])
            samples = np.array(value_dict[metric][model])[8:].reshape((4, 5))
            for k, pos in enumerate(indices + model_idx * width):
                y = samples[k]
                x = pos + np.linspace(-0.25, 0.25, len(y)) * width
                ax.scatter(x, y, color='black', s=1, zorder=3,marker='.')
        if i == 0:
            ax.set_ylim(0.1,0.26)
            ax.set_title('Predicting transcriptional outcomes of genetic perturbation', fontsize = 7)
            ax.set_yticks(ticks = np.arange(2, 6) * 0.05)
            ax.legend(bbox_to_anchor=(1, -0.2), loc='center', ncol=6, fontsize=7)
            ax.get_legend().draw_frame(False)
        if i == 1:
            ax.set_ylim(0.5,0.95)
            ax.set_yticks(ticks = np.arange(5, 10) * 0.1)

        ax.set_ylabel(f'{metric}', fontsize = 7)
        ax.set_xticks(ticks=indices + width * (len(models) - 1) /2)
        ax.set_xticklabels(scenario_names, fontsize = 5)

    if j in {2,3,4}:
        i = j - 2
        indices = np.arange(len(models))
        values = [data2, data1, data3][i]
        errors = [error2, error1, error1 * 0][i]
        for model_idx, model in enumerate(models):
            if i < 2:
                print(i)
                ax.bar([width * model_idx], values[model_idx], yerr=errors[model_idx], width=width, color = colors[model_idx]) #, color=colors[model_idx % len(colors)])
                samples = [values[3 + model_idx * 5: 3 + (model_idx+1) * 5]]
                for k, pos in enumerate([width * model_idx]):
                    y = samples[k]
                    if len(y) < 5:
                        raise ValueError('len(samples) < 5')
                    x = pos + np.linspace(-0.25, 0.25, len(y)) * width
                    ax.scatter(x, y, color='black', s=1, zorder=3,marker='.')
            if i == 2:
                print(i)
                ax.bar([width * model_idx], values[model_idx], width=width, color = colors[model_idx])
        if i == 0:
            ax.set_ylim((1, 1.38))
            ax.set_yticks(ticks = [1.0, 1.1, 1.2, 1.3])
            ax.set_title('GRN inference', fontsize = 7)
        elif i == 1:
            ax.set_ylim((0.9, 1.28))
            ax.set_yticks(ticks = [0.9, 1.0, 1.1, 1.2])
        elif i == 2:
            ax.set_ylim((0.85, 0.98))
            ax.set_yticks(ticks = [0.85,0.9, 0.95])
            ax.set_title('Batch integration', fontsize = 7)
        ax.set_ylabel(['AUPR', 'EPR', 'Batch ASW'][i], fontsize = 7)
        ax.set_xticks(ticks= [width])
        ax.set_xticklabels([], fontsize = 7)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

fig.subplots_adjust(hspace=0, bottom=0.25)


fig.savefig(f'figs/all_abl_low.svg', format='svg')
fig = None
plt.close()
