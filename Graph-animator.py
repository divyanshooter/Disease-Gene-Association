import matplotlib.pyplot as plt
import numpy as np
indices = ['Jaccard','Simpson','Geometric','Cosine']
def plot_bar_x(k_data,auc_data,x):
    plt.bar(k_data,auc_data)
    plt.xlabel('k-data', fontsize=15)
    plt.ylabel('Accuracy of data', fontsize=15)
    plt.title(indices[x]+' Index')
    plt.show()
        
k_data = [5,6,7,8,9,10,11,12,13,14,15,16]
auc_data_jaccard = [42,53.3,69,77,83,85.5,87.21,89.32,88.5,86,83.9,81]
auc_data_simpson = [40.3,52.5,67.6,75.8,82.9,83.3,84.53,85.23,86.51,84.23,85.6,79.6]
auc_data_geometric = [41.6,51.65,66.517,73.694,81.456,81.369,82.542,85.653,88.253,83.369,84.695,80.453]
auc_data_cosine = [40.978,51.675,63.785,72.368,79.145,81.598,81.458,82.345,87.154,82.694,84.942,80.265]
plot_bar_x(k_data,auc_data_jaccard,0)
plot_bar_x(k_data,auc_data_simpson,1)
plot_bar_x(k_data,auc_data_geometric,2)
plot_bar_x(k_data,auc_data_cosine,3)

index = np.arange(4)
bar_width = 0.2
plt.bar(k_data,auc_data_jaccard,bar_width,color='b',label='Jaccard')
for i in range(12):
    k_data[i]+=bar_width
plt.bar(k_data,auc_data_simpson,bar_width,color='g',label='Simpson')
for i in range(12):
    k_data[i]+=bar_width
plt.bar(k_data,auc_data_geometric,bar_width,color='r',label='Geometric')
for i in range(12):
    k_data[i]+=bar_width
plt.bar(k_data,auc_data_cosine,bar_width,color=(0.2, 0.4, 0.6, 0.6),label='Cosine')
plt.xlabel('k-data', fontsize=15)
plt.ylabel('Accuracy of data', fontsize=15)
plt.title('Indices')
plt.show()

import matplotlib.pyplot as plt
import numpy as np

false_positive_values = np.array([0.2,0.35,0.5,0.65,0.8])
true_positive_values = np.array([0.3,0.59,0.69,0.8,0.75])
plt.plot(false_positive_values,true_positive_values)
plt.ylabel('True Positive Values', fontsize=15)
plt.xlabel('False Positive Values', fontsize=15)
plt.title('ROC')
plt.show()

