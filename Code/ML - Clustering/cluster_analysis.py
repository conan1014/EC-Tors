from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools
import pickle
import pandas as pd
from matplotlib.pyplot import cm





n_clusters = 12
filename_imp = 'clusters_all_large.pk'.format(i)
with open(filename_imp, 'rb') as fi:
    df = pickle.load(fi)
  
cols = ['a', 'b', 'No_1', 'No_3', 'No_9', 'Cluster']  

colors1 = iter(cm.rainbow(np.linspace(0, 1, n_clusters)))

for i in range(n_clusters):   
    c = next(colors1)
    cluster = df[['a', 'b', 'No_1', 'No_3', 'No_9', 'Cluster']].loc[df['Cluster']==i,:]
    plt.scatter(cluster['No_1'], cluster['No_3'], c=c, marker='.')

    
plt.show()
    

#group = 7

for ind in df.index:
    if df['No_3'][ind]>0.3 and df['No_1'][ind]>0.45:
            df['Cluster'][ind]=0
    elif df['Cluster'][ind]==7:
        if df['No_1'][ind]<0.1:
            df['Cluster'][ind]=2
        else:
            df['Cluster'][ind]=5


colors = ['y', 'b', 'r', 'c', 'm', 'g', 'k', 'k', 'k', 'k', 'k']
c=0


for i in range(n_clusters):   
    cluster = df[['a', 'b', 'No_1', 'No_3', 'No_9', 'Cluster']].loc[df['Cluster']==i,:]
    plt.scatter(cluster['No_1'], cluster['No_3'], c=colors[c], marker='.')
    if not cluster.empty:
        c+=1
    print(cluster[0:10])
        
plt.show()



means = []
lens = []
c=0
for i in range(n_clusters):
    cluster = df[['a', 'b', 'No_1', 'No_3', 'No_9', 'Cluster']].loc[df['Cluster']==i,:]
    means.append([len(cluster), cluster['No_1'].mean(),cluster['No_3'].mean(),cluster['No_9'].mean()])
    plt.plot(cluster['a'], cluster['b'],  '.', c=colors[c])
    if not cluster.empty:
        c+=1
    
for x in means:
    print(x)

    

'''
#METHODS TO MERGE/SPLIT RAW CLUSTER FILE

n_clusters = 10
filename_imp = 'clusters_all.pk'.format(i)
with open(filename_imp, 'rb') as fi:
    df = pickle.load(fi)
  
cols = ['a', 'b', 'No_1', 'No_3', 'No_9', 'Cluster']  

colors1 = iter(cm.rainbow(np.linspace(0, 1, n_clusters)))

for i in range(n_clusters):   
    c = next(colors1)
    cluster = df[['a', 'b', 'No_1', 'No_3', 'No_9', 'Cluster']].loc[df['Cluster']==i,:]
    plt.scatter(cluster['No_1'], cluster['No_3'], c=c, marker='.')

    
plt.show()
    

for ind in df.index:
    if df['No_3'][ind]>0.3 and df['No_1'][ind]>0.45:
            df['Cluster'][ind]=0
    elif df['Cluster'][ind]==7:
        if df['No_1'][ind]<0.35:
            df['Cluster'][ind]=9
'''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    