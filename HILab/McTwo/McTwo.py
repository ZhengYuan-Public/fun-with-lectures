#Import Libraries
import os
import sys
import platform
import numpy as np
import pandas as pd
import scipy as sp
import minepy
import sklearn
from minepy import MINE
print('-'*100)
print('Software Information:')
print('Python Version =',platform.python_version())  
print('Numpy Version =',np.__version__)
print('Pandas Version =',pd.__version__)
print('Minepy Version =',minepy.__version__)
print('Scipy Version =',sp.__version__)
print('Scikit-learn Version = ',sklearn.__version__)

'''----McOne----'''
print('-'*100)
print('McOne Started')
def fLoadDataMatrix(FilePath):
    with open(FilePath) as file_object:
        contents = pd.read_csv(FilePath,sep='\t')
    return (contents)
def minfo(matrix,name):
    '''Display Matrix Information, it receives two parameters(matrix,matrix_name)'''
    print(f'\nMatrix ---*{name}*--- Information')
    try:
        print("Type:",type(matrix))
    except Exception:
        pass
    try:
        print("Data type:",matrix.dtype)
    except Exception:
        pass
    try:
        print("Size:",matrix.size)
    except Exception:
        pass
    try:
        print("Shape:",matrix.shape)
    except Exception:
        pass    
    try:
        print("Dimission:",matrix.ndim)
    except Exception:
        pass  
    #print(matrix)
    print('-'*40)
def mic(x,y):
    '''Calculate the Maximal Information Coffecient'''
    mine = MINE(alpha=0.6, c=15, est="mic_approx")
    mine.compute_score(x,y)
    return(mine.mic())

#Super Parameters & Load Data
r = 0.33                                              #r is threshold
data_set_path = 'Data\ALL3.txt'
df = fLoadDataMatrix(data_set_path)
#print('\nVisualize Original DataSet:\n',df)
#print('-'*100)

#Step 0: Calculate MIC for each gene
f_array = df.iloc[:,1:].to_numpy()                    #Feature array, converted to numpy array
c_array = df.columns.str.contains('POS*')==True       
c_mic = c_array[1:]                                   #Class array, boolen 
MIC = []                                              #MIC array stores the MIC value for each row(gene)
f_sel_index = []                                      #The index for selected feature
for i in range(len(df)):
    f_i = f_array[i,:]
    mic_i= mic(f_i,c_mic)
    MIC.append(mic_i)
    if mic_i > r:
        f_sel_index.append(i)
print('\nCompute MICs finished.')
print('\nThe Maixmal MIC =',max(MIC))

'''------Rank Selectred Features using MIC------'''
#Step 1: Sort f_sel_index in descending order
df_new = df[['mdr']].copy()
df_new['MIC_VALUE'] = MIC                           #Create a new DataFrame with only Gene names and its MIC values respectivly.                 
df_new_selected = df_new.iloc[f_sel_index,:]        #Use the f_sel_index to choose the features from Step 1
df_new_selected_sorted = df_new_selected.sort_values(by='MIC_VALUE',ascending = False)  #Sort it into Descending Order
f_sel_index_sorted = df_new_selected_sorted.index   #New index for selected&sorted features

#Step 2: Delete Redundant Features
iter_num = len(f_sel_index_sorted)
dropped_features_index = []
for e in range(iter_num):
    q = e +1
    while q < iter_num:
        mic_eq = mic(f_array[f_sel_index_sorted[e]],f_array[f_sel_index_sorted[q]])
        mic_q = MIC[f_sel_index[q]]
        if mic_eq >= mic_q:
            df_new_selected_sorted = df_new_selected_sorted.drop(df_new_selected_sorted.index[q]) 
            #Drop feature indexed as [q in f_sel_index], which is indexed as [f_sel_index[q] in df_new_selected_sorted]
            iter_num -= 1
            dropped_features_index.append(f_sel_index_sorted[q])
            #print(f'Drop feature with index [{f_sel_index_sorted[q]}].')
            f_sel_index_sorted = f_sel_index_sorted.delete(q)
            #print(f_sel_index_sorted)
            #print('New selected feature arrya after droping:\n',df_new_selected_sorted)
        else:
            q += 1
print(f'\nFeatures selected by McOne under threshold r = {r} are:\n',df_new_selected_sorted['mdr'])
print('\nMcOne Finshed')
print('-'*100)
'''----END of McOne----'''

'''----McTwo----'''
print('-'*100)
print('McTwo Started')
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import balanced_accuracy_score

def calculate_BAcc(features_array):
    '''For given feature(s), use NN to calculate its BAcc_values with respect to their true lables. Leave-one-out'''
    features_array = features_array.T                      #Features subset used to calculate BAcc
    true_labels = 1*c_mic.T                                #true_labels =  Data set cloumns converted to boolean then 0/1.
    loo = LeaveOneOut()
    loo.get_n_splits(features_array)                       #The number of iterations.
    classifier = KNeighborsClassifier(n_neighbors = 1)     #KNN, k = 1
    y_pred = []
    for train_index, test_index in loo.split(features_array):
        X_train, X_test = features_array.iloc[train_index], features_array.iloc[test_index]
        y_train, y_test = true_labels[train_index], true_labels[test_index]
        classifier.fit(X_train, y_train)
        y_pred_i = classifier.predict(X_test)
        y_pred.append(y_pred_i)
    BAcc = balanced_accuracy_score(true_labels, y_pred)
    return(BAcc)
def return_selected_features(selected_feature_index_array):
    '''Function name should be different from selected_feature(pd.Dataframe)'''
    try:
        selected_features = df.iloc[selected_feature_index_array,1:].to_frame().T
        #If only one row is selected, pandas will convert it to a series column.
        #Convert series into dataframe then transpose.
    except Exception:
        selected_features = df.iloc[selected_feature_index_array,1:]
    return(selected_features)

'''Best-first-search'''
all_features = df.iloc[f_sel_index_sorted,1:]
selected_features_index = []                                                        #Only Manipulate the index
selected_features = return_selected_features(selected_features_index)

#Step 0: 1st iteration
BAcc_array = []
for i in range(len(all_features)):
    test_feature_i = all_features.iloc[i].to_frame().T
    BAcc_array.append(calculate_BAcc(test_feature_i))
Max_BAcc = max(BAcc_array)
selected_features_index.append(all_features.index[np.argmax(BAcc_array)])

#Step 1: The rest iterations
flag = True
while flag == True:
    left_features = all_features.drop(selected_features.index)
    BAcc_array = []
    for j in range(len(left_features)):
        #print('ITERATION:',j)
        test_feature_i = left_features.iloc[j].to_frame().T
        selected_features = return_selected_features(selected_features_index)
        temp_combination = selected_features.append(test_feature_i)
        test_included_bacc = calculate_BAcc(temp_combination)
        BAcc_array.append(test_included_bacc)
        if test_included_bacc > Max_BAcc:
            Max_BAcc = test_included_bacc
            selected_features_index = temp_combination.index
            print('\nNew Improvement')
            print('Temporary maximal BAcc: ',Max_BAcc)
            #print('Indices for selected features',selected_features_index)
        #print('-'*40)
    if Max_BAcc > max(BAcc_array):
        '''Algorithm stops if Max_BAcc stops increasing.'''
        flag = False
print('\nFeatures selected by McTwo are: \n',df.iloc[return_selected_features(selected_features_index).index]['mdr'])
print('\nThe best BAcc: ',Max_BAcc)
print('\nMcTwo Finshed')
print('-'*100)
'''----END of McTwo----'''