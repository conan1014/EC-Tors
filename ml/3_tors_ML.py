from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools
import pickle
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from matplotlib.pyplot import cm


P, n = 5000, 50
filename_imp = 'Distrib_data_P={0}_n={1}.pk'.format(P, n)

with open(filename_imp, 'rb') as fi:
    arr = pickle.load(fi)
    
    


plt.scatter(arr[0:,2], arr[:,3], marker='.')
plt.title("Distribution scatter graph for P = {0}, n ={1}".format(P, n))
plt.xlabel('#E[3] = 1 proportion')
plt.ylabel('#E[3] = 3 proportion')
plt.show()

df = pd.DataFrame()

df['a'] = arr[0:,0]
df['b'] = arr[0:,1]
df['No_1'] = arr [0:,2]
df['No_3'] = arr [0:,3]
df['No_9'] = arr [0:,4]



df.info()

cont_features = ['No_1', 'No_3', 'No_9']

df_filt = df[cont_features]

mms = MinMaxScaler().fit(df_filt)
X = mms.transform(df_filt)

Sum_of_squared_distances = []
K = range(1,15)
for k in K:
    km = KMeans(n_clusters=k, init="random", random_state=99)
    km = km.fit(X)
    Sum_of_squared_distances.append(km.inertia_)
    
plt.plot(K, Sum_of_squared_distances, 'bx-')
plt.xlabel('k')
plt.ylabel('Sum of Squared Distances')
plt.title('Elbow Method For Optimal k')
plt.show()

n_clusters =8
km = KMeans(n_clusters=n_clusters)
km = km.fit(X)

df['Cluster'] = km.fit_predict(df_filt)

df.head()
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
color = iter(cm.rainbow(np.linspace(0, 1, n_clusters)))


filename_exp = 'clusters_all_large.pk'
with open(filename_exp, 'wb') as fi:
    # dump your data into the file
    pickle.dump(df, fi)

for i in range(n_clusters):
    c = next(color)
    cluster = df[['a', 'b', 'No_1', 'No_3', 'No_9', 'Cluster']].loc[df['Cluster']==i,:]
    plt.scatter(cluster['No_1'], cluster['No_3'], c=c, marker='.')
    print(cluster.head())
    print()

    '''
    filename_exp = 'cluster_{0}.pk'.format(i)
    with open(filename_exp, 'wb') as fi:
        # dump your data into the file
        pickle.dump(cluster, fi)

    '''




























