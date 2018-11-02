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

pca_a = pd.read_csv("pca_b.csv",sep=',',header = None)
x = pca_a.iloc[:, 0 : -1 ].values
tsne = TSNE(n_components=2,n_iter=300)
x_tsne = tsne.fit_transform(x)
d_tsne = {'dim1' : x_tsne[:,0], 'dim2' : x_tsne[:,1], 'Disease': pca_a.iloc[:,-1].values}
df_tsne = pd.DataFrame(data=d_tsne)
ax = plt.axes()
sns.scatterplot(x = 'dim1', y = 'dim2', hue = "Disease", data= df_tsne,ax=ax)
ax.set_title("Plot of tsne pca_a")

