import random
import torch.nn.functional as F
import torch.nn as nn
from models.gat_r import GATNet
from models.gat_gcn_test import GAT_GCN
from models.gcn import GCNNet
from models.ginconv import GINConvNet
import math
from utils_test import *
from sklearn.metrics import roc_curve, confusion_matrix
from sklearn.metrics import cohen_kappa_score, accuracy_score, roc_auc_score, precision_score, recall_score, balanced_accuracy_score
from sklearn import metrics
import pandas as pd



# training function at each epoch
def train(model, device, drug1_loader_train, drug2_loader_train, optimizer, epoch):
    print('Training on {} samples...'.format(len(drug1_loader_train.dataset)))
    model.train()
    # train_loader = np.array(train_loader)
    for batch_idx, data in enumerate(zip(drug1_loader_train, drug2_loader_train)):
        data1 = data[0]
        data2 = data[1]
        data1 = data1.to(device)
        data2 = data2.to(device)
        y = data[0].y.view(-1, 1).float().to(device)
        y = y.squeeze(1)
        output = model(data1, data2)
        loss = loss_fn(output, y.unsqueeze(1))
        # print('loss', loss)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        if batch_idx % LOG_INTERVAL == 0:
            print('Train epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(epoch,
                                                                           batch_idx * len(data1.x),
                                                                           len(drug1_loader_train.dataset),
                                                                           100. * batch_idx / len(drug1_loader_train),
                                                                           loss.item()))


def predicting(model, device, drug1_loader_test, drug2_loader_test):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    total_prelabels = torch.Tensor()
    print('Make prediction for {} samples...'.format(len(drug1_loader_test.dataset)))
    with torch.no_grad():
        for data in zip(drug1_loader_test, drug2_loader_test):
            data1 = data[0]
            data2 = data[1]
            data1 = data1.to(device)
            data2 = data2.to(device)
            output = model(data1, data2)
            ys = output.to('cpu').data.numpy()
            # predicted_labels = list(map(lambda x: np.argmax(x), ys))
            total_preds = torch.cat((total_preds, torch.Tensor(ys)), 0)
            # total_prelabels = torch.cat((total_prelabels, torch.Tensor(predicted_labels)), 0)
            total_labels = torch.cat((total_labels, data1.y.view(-1, 1).cpu()), 0)
    return total_preds.numpy().flatten(), total_labels.numpy().flatten()


def shuffle_dataset(dataset, seed):
    np.random.seed(seed)
    np.random.shuffle(dataset)
    return dataset


def split_dataset(dataset, ratio):
    n = int(len(dataset) * ratio)
    dataset_1, dataset_2 = dataset[:n], dataset[n:]
    return dataset_1, dataset_2

def compute_mae_mse_rmse(target,prediction):
    error = []
    for i in range(len(target)):
        error.append(target[i] - prediction[i])
    squaredError = []
    absError = []
    for val in error:
        squaredError.append(val * val)  # target-prediction之差平方
        absError.append(abs(val))  # 误差绝对值
    mae=sum(absError)/len(absError)  # 平均绝对误差MAE
    mse=sum(squaredError)/len(squaredError)  # 均方误差MSE
    RMSE=sum(absError)/len(absError)
    return mae,mse,RMSE

def compute_rsquared(X, Y):
    xBar = np.mean(X)
    yBar = np.mean(Y)
    SSR = 0
    varX = 0
    varY = 0
    for i in range(0, len(X)):
        diffXXBar = X[i] - xBar
        diffYYBar = Y[i] - yBar
        SSR += (diffXXBar * diffYYBar)
        varX += diffXXBar ** 2
        varY += diffYYBar ** 2

    SST = math.sqrt(varX * varY)
    r2=round((SSR / SST) ** 2,3)
    return r2

# modeling = [GINConvNet, GATNet, GAT_GCN, GCNNet][int(sys.argv[2])]
# datasets = [['davis', 'kiba'][int(sys.argv[1])]]
# model_st = modeling.__name__
modeling = GATNet

TRAIN_BATCH_SIZE = 128
TEST_BATCH_SIZE = 128
LR = 0.0001
LOG_INTERVAL = 20
NUM_EPOCHS = 2000

print('Learning rate: ', LR)
print('Epochs: ', NUM_EPOCHS)
datafile = 'score_0_10'

# CPU or GPU

if torch.cuda.is_available():
    device = torch.device('cuda')
    print('The code uses GPU...')
else:
    device = torch.device('cpu')
    print('The code uses CPU!!!')

drug1_data = TestbedDataset(root='data', dataset=datafile + '_drug1')
drug2_data = TestbedDataset(root='data', dataset=datafile + '_drug2')

lenth = len(drug1_data)
pot = int(lenth/5)
print('lenth', lenth)
print('pot', pot)


random_num = random.sample(range(0, lenth), lenth)
for i in range(5):
    test_num = random_num[pot*i:pot*(i+1)]
    train_num = random_num[:pot*i] + random_num[pot*(i+1):]

    drug1_data_train = drug1_data[train_num]
    drug1_data_test = drug1_data[test_num]
    drug1_loader_train = DataLoader(drug1_data_train, batch_size=TRAIN_BATCH_SIZE, shuffle=None)
    drug1_loader_test = DataLoader(drug1_data_test, batch_size=TRAIN_BATCH_SIZE, shuffle=None)


    drug2_data_test = drug2_data[test_num]
    drug2_data_train = drug2_data[train_num]
    drug2_loader_train = DataLoader(drug2_data_train, batch_size=TRAIN_BATCH_SIZE, shuffle=None)
    drug2_loader_test = DataLoader(drug2_data_test, batch_size=TRAIN_BATCH_SIZE, shuffle=None)

    model = modeling().to(device)
    loss_fn = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)


    # model_file_name = 'data/result/GATNet(DrugA_DrugB)' + str(i) + '--model_' + datafile +  '.model'
    result_file_name = 'data/result/GATNet(DrugA_DrugB)_R' + str(i) + '--result_' + datafile + '.csv'
    file_AUCs = 'data/result/GATNet(DrugA_DrugB)_R' + str(i) + '--AUCs--' + datafile + '.txt'
    AUCs = ('Epoch\tepoch\tmse\tmae\trmse\tr2\tP')
    with open(file_AUCs, 'w') as f:
        f.write(AUCs + '\n')

    best_auc = 10000000000000
    for epoch in range(NUM_EPOCHS):
        train(model, device, drug1_loader_train, drug2_loader_train, optimizer, epoch + 1)
        S, T = predicting(model, device, drug1_loader_test, drug2_loader_test)
        # T is correct score
        # S is predict score
        # Y is predict label
        print(S)
        print(T)

        # compute preformence

        mae, mse, rmse = compute_mae_mse_rmse(S, T)
        
        # save data

        if best_auc > mse:
            r2 = compute_rsquared(S, T)
            best_auc = mse
            AUCs = [epoch, mse, mae, rmse, r2, math.sqrt(r2)]
            save_AUCs(AUCs, file_AUCs)
            # torch.save(model.state_dict(), model_file_name)
            independent_num = []
            # independent_num.append(test_num)
            independent_num.append(T)
            # independent_num.append(Y)
            independent_num.append(S)
            txtDF = pd.DataFrame(data=independent_num)
            txtDF.to_csv(result_file_name, index=False, header=False)

        print('best_mse', best_auc)
