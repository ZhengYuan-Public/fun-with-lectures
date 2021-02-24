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
if tf.test.gpu_device_name() == "/device:GPU:0":
    print ("Hardware = GPU")
else:
    print ("Hardware = CPU")
print ('-------------------------------')

#Import Minist Data Sets
import input_data
print("Extrac data:")
mnist = input_data.read_data_sets("data/",one_hot = True)
train_img = mnist.train.images
train_label = mnist.train.labels
test_img = mnist.test.images
test_label = mnist.test.labels
print ('-------------------------------')

#Define network topologies
layer_0_input_layer = 784 #(?,784)
layer_1_hidden_layer = 256 #256 AN Units
layer_2_hidden_layer = 128 #128 AN Units
layer_3_output_layer = 10 #10 Classes
#(?*784)x(784*256)x(256*128)x(128*10) = (?*10)

#Define I/O
x = tf.placeholder("float",[None,layer_0_input_layer]) #(?*784)
y = tf.placeholder("float",[None,layer_3_output_layer]) #(?*10)
#pint("Input: ",x)
#pint("Output: ",y)

#Initialize Network Parameters
StandardDeviation = 0.1
weights = {
    "w01": tf.Variable(tf.random_normal([layer_0_input_layer,layer_1_hidden_layer], stddev = StandardDeviation)),
    "w12": tf.Variable(tf.random_normal([layer_1_hidden_layer,layer_2_hidden_layer], stddev = StandardDeviation)),
    "w23": tf.Variable(tf.random_normal([layer_2_hidden_layer,layer_3_output_layer], stddev = StandardDeviation)),
    #w01=(784*256); w12=(256*128); w23=(128*10)
}
biases = {
    "b1":tf.Variable(tf.random_normal([layer_1_hidden_layer], stddev = StandardDeviation)),
    "b2":tf.Variable(tf.random_normal([layer_2_hidden_layer], stddev = StandardDeviation)),
    "b3":tf.Variable(tf.random_normal([layer_3_output_layer], stddev = StandardDeviation)),
    #b1=(256); b2=(128); b3=(10)
}
#pint("Layer to layer weights: ",weights)
#pint("Biases: ",biases)
print("Network Parameters Initialized")
print('-------------------------------')

#Forward Propagation
def forward_propagation(x_i, weights_i, biases_i):
    layer_1_output = tf.nn.sigmoid(tf.add(tf.matmul(x_i,weights_i["w01"]),biases_i["b1"]))
    layer_2_output = tf.nn.sigmoid(tf.add(tf.matmul(layer_1_output,weights_i["w12"]),biases_i["b2"]))
    layer_3_output = tf.add(tf.matmul(layer_2_output,weights_i["w23"]),biases_i["b3"])
    return(layer_3_output)
    #Activation function = Sigmiod

#Cost Functions and Optimization
learning_rate = 0.001
prediction = forward_propagation(x,weights,biases)
#cost_function = tf.reduce_mean(-tf.reduce_sum(y*tf.log(prediction),reduction_indices=1))
cost_function = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(prediction,y)) #tf__version__=1.14.0
#cost_function = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(prediction,y)) #tf__version__=1.14.0
#optimization = tf.optmizers.SGD(learning_rate).minimize(cost_function) #tf__version__=2.0
optimization = tf.train.GradientDescentOptimizer(learning_rate).minimize(cost_function)  #tf__version__=1.14.0

#Calculate accuracy
correct_values = tf.equal(tf.argmax(prediction,1),tf.argmax(y,1))
accuracy = tf.reduce_mean(tf.cast(correct_values,"float"))

#--------------------------------------------#

#Network Super Parameter
training_epochs = 50
batch_size = 10
display_step = 5

print("Start Session:")
#with tf.Session() as sess:
#    sess.run(init_op)
init_op = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init_op)

for epoch in range(training_epochs):
    average_cost = 0
    num_batch = int(mnist.train.num_examples/batch_size)
    for i in range(num_batch):
        batch_xs, batch_ys = mnist.train.next_batch(batch_size)
        feeds = {x:batch_xs, y:batch_ys}
        sess.run(optimization, feed_dict = feeds)
        average_cost += sess.run(cost_function,feed_dict=feeds)/num_batch
    if epoch % display_step == 0:
        feeds_train = feeds
        feeds_test = {x:mnist.test.images, y:mnist.test.labels}
        train_acc = sess.run(accuracy, feed_dict = feeds_train)
        test_acc = sess.run(accuracy, feed_dict = feeds_test )
        print("Epoch: %03d/%03d Cost: %.9f Train_accuracy: %.3f Test_accuracy: %.3f" 
              % (epoch, training_epochs, average_cost, train_acc, test_acc))
print("Session Finished")