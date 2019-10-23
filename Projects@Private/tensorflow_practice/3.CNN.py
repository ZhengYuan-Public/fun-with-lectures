#-----Import library and System information-----#
import os,sys
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
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
print("Extract data:")
mnist = input_data.read_data_sets("data/",one_hot = True)
train_img = mnist.train.images
train_label = mnist.train.labels
test_img = mnist.test.images
test_label = mnist.test.labels
print ('-------------------------------')

#Define network topologies
#layer_0_input_layer = 784 #(?,784)
#layer_1_convolution_layer1
#layer_2_pooling_layer1 
#layer_3_convolution_layer2 
#layer_4_pooling_layer2
#layer_5_fc_layer1
#layer_6_fc_layer2_output = 10 #(?,10)
layer_0_input_layer = 784 #(?,784)
layer_6_fc_layer2_output = 10 #(?,10)

#Initialize Network Parameters
StandardDeviation = 0.1
weights = {
    #Weights for convolution layer
    "w01": tf.Variable(tf.random.normal([3,3,1,64], stddev = StandardDeviation)),
    #filter: [Height,Width,Depth,num_FeatureMaps]
    "w23": tf.Variable(tf.random.normal([3,3,64,128], stddev = StandardDeviation)),
    #filter: Depth = 64, num_FeatrureMaps = 128
    
    #Weights for fully connect layer
    "w45": tf.Variable(tf.random.normal([7*7*128,1024], stddev = StandardDeviation)),
    #Pooling_1 = 14*14; Pooling_2 = 7*7
    "w56": tf.Variable(tf.random.normal([1024,layer_6_fc_layer2_output], stddev = StandardDeviation)),
}
biases = {
    "b01":tf.Variable(tf.random.normal([64], stddev = StandardDeviation)),
    "b23":tf.Variable(tf.random.normal([128], stddev = StandardDeviation)),
    "b45":tf.Variable(tf.random.normal([1024], stddev = StandardDeviation)),
    "b56":tf.Variable(tf.random.normal([layer_6_fc_layer2_output],stddev = StandardDeviation)),
}
#print("Layer to layer weights: ",weights)
#print("Biases: ",biases)
print("Network Parameters Initialized")
print('-------------------------------')

def forward_propagation_conv(input_i, weights_i, biases_i, droupout_ratio):
    #input
    input_reshaped = tf.reshape(input_i, shape = [-1,28,28,1])
    #shape = [Batch_size,Width,height,Depth]
    #Convolution layer 1
    conv1 = tf.nn.conv2d(input_reshaped, weights_i['w01'], strides = [1,1,1,1], padding = 'SAME')
    conv1_output = tf.nn.relu(tf.nn.bias_add(conv1,biases_i['b01']))
    #Pooling layer 1
    pool1 = tf.nn.max_pool(conv1_output, ksize = [1,2,2,1], strides = [1,2,2,1], padding = 'SAME')
    pool1_dropout = tf.nn.dropout(pool1, droupout_ratio)
    #Convolution layer 2
    conv2 = tf.nn.conv2d(pool1_dropout, weights_i['w23'], strides = [1,1,1,1], padding = 'SAME')
    conv2_output = tf.nn.relu(tf.nn.bias_add(conv2,biases_i['b23']))
    #Pooling layer 2
    pool2 = tf.nn.max_pool(conv2_output, ksize = [1,2,2,1], strides = [1,2,2,1], padding = 'SAME')
    pool2_dropout = tf.nn.dropout(pool2, droupout_ratio)
    #Fully Connected Layer 1
    pool2_dropout_reshaped = tf.reshape(pool2_dropout,[-1,weights_i['w45'].get_shape().as_list()[0]])
    fc1 = tf.nn.relu(tf.add(tf.matmul(pool2_dropout_reshaped,weights_i['w45']),biases_i['b45']))
    fc1_dropout = tf.nn.dropout(fc1,droupout_ratio)
    #Fully Connected Layer 2
    fc2_output = tf.add(tf.matmul(fc1_dropout,weights_i['w56']),biases_i['b56'])
    return(fc2_output)
    #return(
    #   {
    #        input_reshaped,
    #        conv1,conv1_output,
    #        pool1,pool1_dropout,
    #        conv2,conv2_output,
    #        pool2,pool2_dropout,
    #        fc1,fc1_dropout,
    #        fc2_output
    #    }
    #      )

#Define I/O
x = tf.placeholder(tf.float32,[None,layer_0_input_layer])
y = tf.placeholder(tf.float32,[None,layer_6_fc_layer2_output])
droupout_ratio = tf.placeholder(tf.float32)
#print("Input: ",x)
#print("Output: ",y)
#print("Droupout_ratio: ",droupout_ratio)

#Prediction
prediction = forward_propagation_conv(x, weights, biases, droupout_ratio)
#print(prediction)

#Cost Functions and Optimization
learning_rate = 0.001
cost_function = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(prediction,y))
optimization = tf.train.AdamOptimizer(learning_rate).minimize(cost_function)

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