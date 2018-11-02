import os
import csv
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA 
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.decomposition import TruncatedSVD


# input_files = [files for files in os.listdir('.')]
# print(files)
'''with open('pca_a.txt', 'r') as in_file:
	stripped = (line.strip() for line in in_file)
	lines = (line.split(",") for line in stripped if line)
	with open('pca_a.csv', 'w') as out_file:
		writer = csv.writer(out_file)
		# writer.writerow(('a1', 'a2','a3','a4','a5'))
		writer.writerows(lines)
        '''

pca_a = pd.read_csv("pca_a.csv",sep=',',header = None)
x = pca_a.iloc[:, 0 : -1 ].values

#pca = PCA(n_components = 2)
#x_trans = pca.fit_transform(x)
#d1 = {'pca1' : x_trans[:,0], 'pca2' : x_trans[:,1], 'Disease': pca_a.iloc[:,-1].values}
#df1 = pd.DataFrame(data = d1)
#ax = plt.axes()
#sns.scatterplot(x = 'pca1', y = 'pca2', hue = "Disease", data= df1,ax=ax)
#ax.set_title("Plot of Pca ")

x -= x.mean(axis = 0)
cov_mat = np.cov(x,rowvar=False)
eig_val_sc, eig_vec_sc = np.linalg.eig(cov_mat)
idx = eig_val_sc.argsort()
idx = idx[:: -1]
eig_val_sc = eig_val_sc[idx]
eig_vec_sc = eig_vec_sc[:,idx]
y = np.dot(x,eig_vec_sc)
eig_vec_sc = (eig_vec_sc.transpose())*-1
var = np.dot(x, eig_vec_sc[:2, :].transpose())
print(var)
d_pca = {'pca1' : var[:,0], 'pca2' : var[:,1], 'Disease': pca_a.iloc[:,-1].values}
df_pca = pd.DataFrame(data = d_pca)
ax = plt.axes()
sns.scatterplot(x = 'pca1', y = 'pca2', hue = "Disease", data= df_pca,ax=ax)
ax.set_title("Plot of PCA_C ")
plt.show()

