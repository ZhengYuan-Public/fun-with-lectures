#-----Import library and System information-----#
import os,sys
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
#import keras
print ("System info:")
print ("Python version =", sys.version)
print ("Tensorflow version =",tf.__version__)
#print ("keras version =",keras.__version__)
if tf.test.gpu_device_name() == '/device:GPU:0':
    print ("Hardware = GPU")
else:
    print ("Hardware = CPU")
print ('-------------------------------')
#--------------------------------------------#

#Import Minist Data Sets
import input_data
print("Extrac data:")
mnist = input_data.read_data_sets('data/',one_hot = True)
train_img = mnist.train.images
train_label = mnist.train.labels
test_img = mnist.test.images
test_label = mnist.test.labels
print ('-------------------------------')

#Visualize Data Sets
#print("Visualize Data Sets:")
#n_sample = 2
#rand_index = np.random.randint(train_img.shape[0],size = n_sample)
#for i in rand_index:
#  curr_img = np.reshape(train_img[i,:],(28,28))
#  curr_label = np.argmax(train_label[i,:])
#  plt.matshow(curr_img,cmap=plt.get_cmap('gray'))
#  plt.title("The "+str(i)+"th Training Data"+"label is "+str(curr_label))
#  print("The "+str(i)+"th Training Data"+"label is "+str(curr_label))
#  plt.show()
#print ('-------------------------------')

#Logistic Regression

#print("train_img shape =",train_img.shape)
#print("train_label shape =",train_label.shape)
#print("test_img shape =",test_img.shape)
#print("test_label shape =",test_label.shape)

x = tf.placeholder("float",[None, 784])
y = tf.placeholder("float",[None, 10])
W = tf.Variable(tf.zeros([784,10]))
b = tf.Variable(tf.zeros([10]))

#There was a mistake for using tf.zeros(784,10), 
#which cause it's value to be false.

#print("x =",x)
#print("y =",y)
#print("W =",W)
#print("b =",b)

#Activation Function
#print(tf.matmul(x,W).shape)
activ_func = tf.nn.softmax(tf.matmul(x,W)+b)
#print(activation.shape)
cost_func = tf.reduce_mean(-tf.reduce_sum(y*tf.log(activ_func),reduction_indices=1))
#optimization = tf.train.GradientDescentOptmizer(learning_rate).minimize(cost_funct)
#optimization = tf.optmizers.SGD(learning_rate).minimize(cost_funct) #tf__version__=2.0
learning_rate = 0.01
optimization = tf.train.GradientDescentOptimizer(learning_rate).minimize(cost_func) #tf__version__=1.14.0
prediction = tf.equal(tf.argmax(activ_func,1),tf.argmax(y,1))
accuracy = tf.reduce_mean(tf.cast(prediction,"float")) #tf.cast(): Converto Boolean to float.
#--------------------------------------------#

#Training
training_epochs = 500
batch_size = 128
display_step = 5

#Start Session
print("Start Session:")

#with tf.Session() as sess:
#    sess.run(init_op)
init_op = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init_op)
for epoch in range(training_epochs): #training_epochs useage
    average_cost = 0
    num_batch = int(mnist.train.num_examples/batch_size)
    for i in range(num_batch):
        batch_xs, batch_ys = mnist.train.next_batch(batch_size) #learm this useage
        feeds = {x:batch_xs, y:batch_ys} #Load current batch to feeds
        sess.run(optimization, feed_dict = feeds) # {x:batch_xs, y:batch_ys})
        average_cost += sess.run(cost_func,feed_dict=feeds)/num_batch
    if epoch % display_step == 0:
        feeds_train = feeds
        feeds_test = {x:mnist.test.images, y:mnist.test.labels}
        train_acc = sess.run(accuracy, feed_dict = feeds_train)
        test_acc = sess.run(accuracy, feed_dict = feeds_test )
        print("Epoch: %03d/%03d Cost: %.9f Train_accuracy: %.3f Test_accuracy: %.3f" 
              % (epoch, training_epochs, average_cost, train_acc, test_acc))
print("Session Finished")