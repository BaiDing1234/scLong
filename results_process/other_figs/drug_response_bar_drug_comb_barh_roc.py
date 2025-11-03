import matplotlib.pyplot as plt
import numpy as np
import re
import pickle


plt.rcParams['font.family'] = 'DejaVu Sans'

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = 6
plt.rcParams['svg.fonttype'] = 'none' 
plt.rcParams['lines.linewidth'] = 0.5

with open('../pkls/drug_and_drug_comb_dict.pkl', 'rb') as file:
    drug_and_drug_comb_dict = pickle.load(file)

curve_dir = '/home/ding.bai/ding-scfmv1-downstream/update/results_process/other_figs/drug_comb_roc'


fig = plt.figure()  # 设置整体图形大小
gs = fig.add_gridspec(1, 15)
axes = [fig.add_subplot(gs[0, 0:4]),
       fig.add_subplot(gs[0, 5:9]),
       fig.add_subplot(gs[0, 10:15])]

#fig, axes = plt.subplots(1,3)

fig.set_size_inches(7.2, 2.6)


for ax_i, ax in enumerate(axes):
    if ax_i == 0:

        
        data1 = np.array(drug_and_drug_comb_dict['drug']['data1'])
        errors = np.array(drug_and_drug_comb_dict['drug']['errors'])

        color_dict = {'DeepCDR': '#FFCC99',
              'Geneformer': '#99CC99',
              'scGPT': '#FF9999',
              'scFoundation': '#D1B26F',
              'UCE': '#99FFCC',
              'scLong': '#CC99FF'}


        models = ['DeepCDR', 'Geneformer','scGPT', 'scFoundation',  'UCE', 'scLong']

        colors = [color_dict[model] for model in models]
        indices = np.arange(len(models)) * 2

        width = 1 / (len(models) + 1)
        values = data1
        for model_idx, model in enumerate(models):
            #ax.bar([(indices * width)[model_idx]], width=values[model_idx], xerr=errors[model_idx], width=width, color = colors[model_idx]) #, color=colors[model_idx % len(colors)])
            ax.barh([-(indices * width)[model_idx]] , width=values[model_idx], xerr=errors[model_idx], height=width * 1.2 , color = colors[model_idx])

            bar_positions = [-(indices * width)[model_idx]]
            samples = [values[6 + model_idx * 5 : 6 + (model_idx+1) * 5],]
            for i, pos in enumerate(bar_positions):
                y = samples[i]
                if len(y) < 5:
                    raise ValueError('len(samples) < 5')
                x = pos + np.linspace(-0.25, 0.25, len(y)) * width* 1.2
                ax.scatter(y, x, color='black', s=1, zorder=3,marker='.')
        ax.set_xlim((0.82, 0.89))
        ax.set_xlabel(f'Pearson correlation', fontsize = 6)
        ax.set_xticks(ticks=0.02 * np.arange(41, 45))
        ax.set_yticks(ticks=-indices * width)
        ax.set_yticklabels(models, fontsize=6)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)


        ax.set_title('Single drug response prediction', fontsize = 9)
        
    elif ax_i == 1:
                
        data1 = np.array(drug_and_drug_comb_dict['drug_comb']['data1'])
        errors = np.array(drug_and_drug_comb_dict['drug_comb']['errors'])

        color_dict = {'DeepDDS': '#FFCC99',
              'Geneformer': '#99CC99',
              'scGPT': '#FF9999',
              'scFoundation': '#D1B26F',
              'UCE': '#99FFCC',
              'scLong': '#CC99FF'}


        models = ['DeepDDS', 'Geneformer','scGPT', 'scFoundation',  'UCE', 'scLong']
        colors = [color_dict[model] for model in models]
        
        indices = np.arange(len(models)) * 2

        width = 1 / (len(models) + 1)
        values = data1
        for model_idx, model in enumerate(models):
            #ax.bar([(indices * width)[model_idx]], values[model_idx], yerr = errors[model_idx], width=width, color = colors[model_idx]) 
            ax.barh([-(indices * width)[model_idx]] , width=values[model_idx], xerr=errors[model_idx], height=width * 1.2 , color = colors[model_idx])
            
            bar_positions = [-(indices * width)[model_idx]]
            samples = [values[6 + model_idx * 5 : 6 + (model_idx+1) * 5],]
            for i, pos in enumerate(bar_positions):
                y = samples[i]
                if len(y) < 5:
                    raise ValueError('len(samples) < 5')
                x = pos + np.linspace(-0.25, 0.25, len(y)) * width* 1.2
                ax.scatter(y, x, color='black', s=1, zorder=3,marker='.')
        
        ax.set_xlim((0.56, 0.66))
        ax.set_xticks(ticks = 0.02 * np.arange(28, 34))
        ax.set_xlabel(f'AUROC', fontsize = 6)
        ax.set_yticks(ticks=-indices * width)
        ax.set_yticklabels(models, fontsize = 6)

        
        ax.set_title('Drug combination response prediction', fontsize = 9)

        #ax.legend(loc = 'upper right', fontsize=10,frameon=False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)


    else:

        color_dict = {'DeepDDS': '#FFCC99',
              'Geneformer': '#99CC99',
              'scGPT': '#FF9999',
              'scFoundation': '#D1B26F',
              'UCE': '#99FFCC',
              'scLong': '#CC99FF'}

        arrays_dict = {}
        
        with open(f'{curve_dir}/fpr_tpr_deepdds_mean.pkl', 'rb') as f: 
            arrays_dict['DeepDDS'] = pickle.load(f)
        auroc = np.trapz(arrays_dict['DeepDDS'][1], arrays_dict['DeepDDS'][0])
        print(f"DeepDDS: AUROC:", auroc)

        with open(f'{curve_dir}/fpr_tpr_geneformer_mean.pkl', 'rb') as f: 
            arrays_dict['Geneformer'] = pickle.load(f)
        auroc = np.trapz(arrays_dict['Geneformer'][1], arrays_dict['Geneformer'][0])
        print(f"Geneformer: AUROC:", auroc)

        with open(f'{curve_dir}/fpr_tpr_scgpt_mean.pkl', 'rb') as f: 
            arrays_dict['scGPT'] = pickle.load(f)
        auroc = np.trapz(arrays_dict['scGPT'][1], arrays_dict['scGPT'][0])
        print(f"scGPT: AUROC:", auroc)

        with open(f'{curve_dir}/fpr_tpr_scfoundation_mean.pkl', 'rb') as f:
            arrays_dict['scFoundation'] = pickle.load(f)
        auroc = np.trapz(arrays_dict['scFoundation'][1], arrays_dict['scFoundation'][0])
        print(f'scFoundation: {auroc}')

        with open(f'{curve_dir}/fpr_tpr_uce_mean.pkl', 'rb') as f: 
            arrays_dict['UCE'] = pickle.load(f)
        auroc = np.trapz(arrays_dict['UCE'][1], arrays_dict['UCE'][0])
        print(f"UCE: AUROC:", auroc)

        with open(f'{curve_dir}/fpr_tpr_ours_mean.pkl', 'rb') as f: 
            arrays_dict['scLong'] = pickle.load(f)
        auroc = np.trapz(arrays_dict['scLong'][1], arrays_dict['scLong'][0])
        print(f"scLong: AUROC:", auroc)
        
        ax.plot((0,1), (0,1), 'r--')  # 45度红色虚线

        models = ['DeepDDS', 'Geneformer','scGPT', 'scFoundation',  'UCE', 'scLong']

        # 遍历字典中的每个 method 和对应的 array
        for method in models:
        #for method, data in arrays_dict.items():
            data =  arrays_dict[method]
            fpr = data[0]  # 第一行为 FPR
            tpr = data[1]  # 第二行为 TPR

            # 使用NumPy的trapz函数计算曲线下面积，即AUROC
            auroc = np.trapz(tpr, fpr)

            print(f"{method}: AUROC:", auroc)

            # 绘制 ROC 曲线
            ax.plot(fpr, tpr, label=method, color = color_dict[method])

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # 设置 x 和 y 轴标签
        ax.set_xlabel('False positive rate', fontsize = 6)
        ax.set_ylabel('True positive rate', fontsize = 6)

        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))

        ax.set_xticks(np.arange(0, 6)*0.2)
        ax.set_yticks(np.arange(0, 6)*0.2)

        # 设置标题

        # 添加图例
        ax.legend(loc = 'lower right', fontsize=6,frameon=False)
        ax.set_aspect('equal', adjustable='box')

fig.subplots_adjust(hspace=0, wspace=0.2, bottom=0.2)
#fig.tight_layout()
fig.savefig(f'figs/drug_and_drug_comb_res_all_barh.svg', format='svg')
fig = None
plt.close()