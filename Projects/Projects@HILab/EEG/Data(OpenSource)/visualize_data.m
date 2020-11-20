N1 = load('N1.mat');
data = N1.dataStruct.data;
channels = N1.dataStruct.channelIndices;
figure, plot(data(:,channels(5)))