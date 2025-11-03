import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
import pickle

plt.rcParams['font.family'] = 'DejaVu Sans'

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = 9
plt.rcParams['svg.fonttype'] = 'none' 
plt.rcParams['lines.linewidth'] = 0.5


title = "Gene Profile Prediction"
dataset = "DeepCE"


color_dict = {'DeepCE': '#FFCC99',
              'Geneformer': '#99CC99',
              'scGPT': '#FF9999',
              'scFoundation': '#D1B26F',
              'UCE': '#99FFCC',
              'scLong': '#CC99FF'}

with open('../pkls/gene_profile_dict.pkl', 'rb') as file:
    gene_profile_dict = pickle.load(file)

data1 = np.array(gene_profile_dict['data1'])
samples1 = np.array(gene_profile_dict['data1'])
std_err1 = np.array(gene_profile_dict['std_err1'])

data2 = np.array(gene_profile_dict['data2'])
std_err2 = np.array(gene_profile_dict['std_err2'])

data3 = np.array(gene_profile_dict['data3'])
std_err3 = np.array(gene_profile_dict['std_err3'])

data4 = np.array(gene_profile_dict['data4'])
std_err4 = np.array(gene_profile_dict['std_err4'])

data5 = np.array(gene_profile_dict['data5'])
std_err5 = np.array(gene_profile_dict['std_err5'])

y1 = np.mean(data1, axis=1)
y2 = np.mean(data2, axis=1)
y3 = np.mean(data3, axis=1)
y4 = np.mean(data4, axis=1)
y5 = np.mean(data5, axis=1)

#ys = [y1, y2, y3, y4, y5]
ys = [y5, y1, y2, y3, y4]

print(y5.shape)



#std_errs = [std_err1, std_err2, std_err3, std_err4, std_err5]
std_errs = [std_err5, std_err1, std_err2, std_err3, std_err4]

models = ['DeepCE', 'Geneformer', 'scGPT', 'scFoundation', 'UCE', 'scLong']


#######################
#, 'Spearman', 'Pearson'


fig = plt.figure(figsize=(7.2, 1.8))  # 设置整体图形大小
gs = fig.add_gridspec(1, 8)
axs = [fig.add_subplot(gs[0, 0:2]),
       fig.add_subplot(gs[0, 2:5]),
       fig.add_subplot(gs[0, 5:8])]

scenario_names_sets = [['Root mean squared error'], ["Spearman", "Pearson"], ['Pos-P@100', 'Neg-P@100'], ]
idxes = [[0,1],[1, 3], [3, 5]]

for ax_i, ax in enumerate(axs):
    width = 1 / (len(models) + 1)
    if ax_i != 0:
        indices = np.arange(2)
    else:
        indices = np.arange(1)

    start = idxes[ax_i][0]
    end = idxes[ax_i][1]

    for model_idx, model in enumerate(models):
        values = np.array([y[model_idx] for y in ys[start: end]])
        errors = np.array([std_err[model_idx] for std_err in std_errs[start: end]])
        samples = np.array([y[6 + model_idx * 5 : 6 + (model_idx+1) * 5] for y in ys[start: end]])
        label = model

        bar_positions = indices + model_idx * width
        ax.bar(bar_positions, values, yerr=errors, label=label, width=width, color = color_dict[model]) #, color=colors[model_idx % len(colors)])
        
        for i, pos in enumerate(bar_positions):
            y = samples[i]
            if len(y) < 5:
                raise ValueError('len(samples) < 5')
            x = pos + np.linspace(-0.25, 0.25, len(y)) * width
            ax.scatter(x, y, color='black', s=1, zorder=3,marker='.')

    scenario_names = scenario_names_sets[ax_i]
    if ax_i == 1:
        ax.set_ylim((0.4, 0.48))
        ax.set_yticks(ticks = np.arange(20, 25)*0.02)
    elif ax_i == 2:
        ax.set_ylim((0.2, 0.32))
        ax.set_yticks(ticks = np.arange(10, 17)*0.02)
    else:
        ax.set_ylim((1.7, 1.78))
        ax.set_yticks(ticks = np.arange(85, 90)*0.02)

    ax.set_xticks(ticks=indices + width * (len(models) - 1) /2)
    ax.set_xticklabels(scenario_names, fontsize=9)

    #ax.legend(loc = 'upper left', fontsize=9,frameon=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    #fig.subplots_adjust(hspace=0, wspace=0.1)

#ax.set_title('Gene profile prediction', fontsize = 18)
fig.tight_layout()
fig.savefig(f'figs/gene_profile_all.svg', format='svg')
fig = None
plt.close()