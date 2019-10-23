def minfo(matrix):
    '''Display Matrix Information, it receives two parameters(matrix,matrix_name)'''
    print('-'*19+'Start'+'-'*19)
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
    print(matrix)
    print('-'*20+'End'+'-'*20)

#Load Dataset
from scipy.io import loadmat
data_all = loadmat('N1.mat')['dataStruct']
#minfo(data_all)
channels_all = data_all['data'][0][0]
channels_index = data_all['channelIndices'][0][0][0]
sampling_rate = data_all['iEEGsamplingRate'][0][0][0][0]
samples_num = data_all['nSamplesSegment'][0][0][0][0]
sequence = data_all['sequence'][0][0][0][0]

#Visualization Data
import numpy as np
import matplotlib.pyplot as plt
time_start = 0
time_stop = samples_num//sampling_rate
time_step = 1/sampling_rate
x = np.arange(time_start,time_stop,time_step) #np.arange(), np.range(), np.linspace()
y = channels_all
fig, ax = plt.subplots()
ax.plot(x,y)
plt.show()