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

### for svd
svd = TruncatedSVD(n_components=2, n_iter=7, random_state=42)
x_svd = svd.fit_transform(x)
print(x_svd.shape)
d_svd = {'dim1' : x_svd[:,0], 'dim2' : x_svd[:,1], 'Disease': pca_a.iloc[:,-1].values}
df_svd = pd.DataFrame(data=d_svd)
ax = plt.axes()
sns.scatterplot(x = 'dim1', y = 'dim2', hue = "Disease", data= df_svd,ax=ax)