import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
import matplotlib.gridspec as gridspec
import pickle


plt.rcParams['font.family'] = 'DejaVu Sans'

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = 7
plt.rcParams['svg.fonttype'] = 'none' 
plt.rcParams['lines.linewidth'] = 0.5

color_dict = {'GENIE3': '#FF99CC',
                'DeepSEM': '#FFCC99',
              'Geneformer': '#99CC99',
              'scGPT': '#FF9999',
              'scFoundation': '#D1B26F',
              'UCE': '#99FFCC',
              'scLong': '#CC99FF'}



with open('../pkls/grn_inference_dict.pkl', 'rb') as file:
    grn_inference_dict = pickle.load(file)
data1 = np.array(grn_inference_dict['data1'])
data2 = np.array(grn_inference_dict['data2'])
error1 = np.array(grn_inference_dict['error1'])
error2 = np.array(grn_inference_dict['error2'])

models = ['GENIE3', 'DeepSEM', 'Geneformer', 'scGPT', 'scFoundation', 'UCE', 'scLong']
colors = [color_dict[model] for model in models]
 
width = 1 / (len(models) + 1)

fig, axes = plt.subplots(1,2)
fig.set_size_inches(7.2, 2.4)

for i in range(2):
    ax = axes[i]
    indices = np.arange(len(models)) * 2.0
    values = [data2, data1][i]
    errors = [error2, error1][i]
    for model_idx, model in enumerate(models):
        ax.barh([-(indices * width)[model_idx]] , width=values[model_idx], xerr=errors[model_idx], height=width * 1.2, color = colors[model_idx]) #, color=colors[model_idx % len(colors)])
        bar_positions = [-(indices * width)[model_idx]]
        samples = [values[7 + model_idx * 5 : 7 + (model_idx+1) * 5],]
        for k, pos in enumerate(bar_positions):
            y = samples[k]
            if len(y) < 5:
                raise ValueError('len(samples) < 5')
            x = pos + np.linspace(-0.25, 0.25, len(y)) * width* 1.2
            ax.scatter(y, x, color='black', s=1, zorder=3,marker='.')
    if i == 0:
        ax.set_xlim((1, 1.4))
        ax.set_xticks(ticks = [1.0, 1.1, 1.2, 1.3, 1.4])
    elif i == 1:
        ax.set_xlim((1.0, 1.25))
        ax.set_xticks(ticks = [1.0, 1.05, 1.10, 1.15, 1.20, 1.25])
    ax.set_xlabel(['AUPR', 'EPR'][i], fontsize = 7)
    #plt.title(f'Histogram of {column} by Model')
    ax.set_yticks(ticks=-indices * width)
    ax.set_yticklabels(models, fontsize=7)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

fig.subplots_adjust(bottom=0.2)

fig.savefig(f'figs/grn_inference_barh.svg',format="svg")
fig = None
plt.close()