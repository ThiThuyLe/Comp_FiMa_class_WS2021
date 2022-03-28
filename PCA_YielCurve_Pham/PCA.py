# ------------------------------------------------------------------------------
# Published in:   Comp_FiMa_class_WS2021
# ------------------------------------------------------------------------------
# Project title:  Principal Component Analysis for Yield Curve Modelling
# ------------------------------------------------------------------------------
# Description:    Applied PCA for daily US Treasure Yields 2018-2022
# ------------------------------------------------------------------------------
# Output:         List pictures, tables or other outputs of the script.
# ------------------------------------------------------------------------------
# Author :        Thi Ngan Ha Pham, 566022
# ------------------------------------------------------------------------------
# Dataset:        List datafiles necessary to reproduce results (Don't forget to
#                 upload them to the Project folder on GitHub)
# ------------------------------------------------------------------------------

"""These code based on the source code: https://colab.research.google.com/drive/1cT2-QftYsK3haS5kn5l6VylhD51qb3Wk#scrollTo=qw8eaBPC8Hzx """

from datetime import date, datetime, timedelta
import pandas as pd
import numpy as np
import scipy
import pandas_datareader.data as pdr
from pandas_datareader.fred import FredReader
import matplotlib.pyplot as plt
import seaborn as sns

# download daily treasury yield from Fred
# pull data from fred through the function pandas_datareader.data, we could receive data from fred website https://fred.stlouisfed.org
# each DSG3MO, DSG60,... means Market Yield on U.S. Treasury Securities at 3-Month, 6-Month,... Constant Maturity (DGS1MO), (DSG3Mo)...
# https://fred.stlouisfed.org/searchresults/?st=%20Market%20Yield%20on%20U.S.%20Treasury%20Securities%20at this website show the data

codes = ['DGS3MO', 'DGS6MO', 'DGS1', 'DGS2', 'DGS3', 'DGS5', 'DGS7', 'DGS10', 'DGS20']
start_date = datetime(2018, 3, 20)
end_date = datetime(2022, 3, 20)
df = pd.DataFrame()

# A Pandas DataFrame is a 2 dimensional data structure, like a 2 dimensional array, or a table with rows and columns.
# we create a dataframe, which combines of these 'codes' datasets to apply pca technique for this work by using pandas.Dataframe()

for code in codes:
    reader = FredReader(code, start_date, end_date)
    df0 = reader.read()
    df = df.merge(df0, how='outer', left_index=True, right_index=True, sort=False) # Merge data into DataFrame
    reader.close()
df.dropna(axis=0, how='any', inplace=True) # drop rows which contain missing value (NA)

# tail function helps show 5 last rows, default =5
print(df.tail())
print(df.shape)

# view the yield curve
plt.figure(figsize=(5, 5))
plt.plot(df)
plt.show()

""" Get mean and standard deviation """
means = {}
stds = {}
for code in codes:
    print("##### Calculate for " + code)
    current_column = df[code]
    means[code] = np.mean(current_column.values)
    print(np.round(means[code], 3))
    stds[code] = np.std(current_column.values)
    print(np.round(stds[code], 3))

"""Result mean and std respectively """
##### Calculate for DGS3MO: 1.05  1.0
##### Calculate for DGS6MO: 1.106  1.021
##### Calculate for DGS1: 1.156  1.035
##### Calculate for DGS2: 1.242  1.021
##### Calculate for DGS3: 1.32  0.989
##### Calculate for DGS5: 1.479  0.909
##### Calculate for DGS7: 1.658  0.839
##### Calculate for DGS10: 1.799  0.783
##### Calculate for DGS20: 2.153  0.643

""""barplot daily mean and standard deviation"""
# means:
# means = {'DGS3MO': 123, 'DGS6MO': 234, ...}
# means.keys() = ['DGS3MO', 'DGS6MO', 'DGS1', 'DGS2', 'DGS3', 'DGS5', 'DGS7', 'DGS10', 'DGS20']
# means.values() = [123, 234, ..]

x_pos = np.arange(len(means.keys()))
plt.bar(x_pos, means.values(), color=(0.5, 0.1, 0.5, 0.6))
plt.xticks(x_pos, means.keys())
plt.title('Daily Mean of US Treasure Yield 2018-2022')
plt.show()

"""standard deviation"""
y_pos = np.arange(len(stds.keys()))
plt.bar(y_pos, stds.values(), color=(0.5, 0.1, 0.5, 0.6))
plt.xticks(y_pos, stds.keys())
plt.title('Daily Volatilities of US Treasure Yield 2018-2022')
plt.show()

# Correlation between features

df_diff = df.diff()
df.diff()
df_diff.dropna(inplace=True)

print(df.shape)
print(df_diff.shape)

sns.heatmap(df_diff.corr(), cmap="YlGnBu", annot=True)

# displaying heatmap
plt.show()

# correlation among tenors
#sns.pairplot(df.)
#plt.show()

""""pca """

from sklearn.decomposition import PCA

pca_X = PCA().fit(df_diff)
pca_Y = PCA().fit(df)

print(pca_X.explained_variance_)        # eigenvalues
print(pca_X.explained_variance_ratio_)     # normalized eigenvalues (sum to 1)
print(np.cumsum(pca_X.explained_variance_ratio_))

""""Result"""
#Proportion of Variance:
#[0.79578655 0.10619035 0.05175344 0.01530251 0.01076603 0.00925707 0.00539401 0.0033451  0.00220494]
#Cumulative Proportion
# [0.79578655 0.9019769  0.95373034 0.96903285 0.97979888 0.98905595  0.99444996 0.99779506 1.        ]

plt.plot(pca_X.explained_variance_ratio_.cumsum())
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance')
plt.show()

# plot component
df_pca_level = pca_Y.transform(df)            # T or PCs
df_pca_level = pd.DataFrame(df_pca_level, columns=[f'PCA_{x+1}' for x in range(df_pca_level.shape[1])])  # np.array to dataframe
df_pca_level.index = df.index
plt.figure(figsize=(15, 8))
plt.plot(df_pca_level['PCA_1'], label='first component')
plt.plot(df_pca_level['PCA_2'], label='second component')
plt.plot(df_pca_level['PCA_3'], label='third component')
plt.legend()
plt.show()





