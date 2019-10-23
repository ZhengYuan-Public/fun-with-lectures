import os,sys
import tensorflow as tf
#import keras
print ("python version", sys.version)
print ("tensorflow version =",tf.__version__)
#print ("keras version =",keras.__version__)
if tf.test.gpu_device_name() == '/device:GPU:0':
    print ("hardware = GPU")
else:
    print ("hardware = CPU")

w = tf.Variable([[0.5,1.0]])
x = tf.Variable([[2.0],[1.0]])
y = tf.matmul(w,x)
init_rvs = tf.global_variables_initializer()
with tf.Session() as sess:
    sess.run(init_rvs)
    print("w =",w.eval())
    print("x =",x.eval())
    print("y =",y.eval())

w_layer1 = tf.zeros([3,4])
with tf.Session() as sess:
    sess.run(w_layer1)
    print("w =",w_layer1.eval())

#Linear Regression
import os,sys
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
print ("python version =",sys.version)
print ("tensorflow version =",tf.__version__)
num_points = 1000
vector = []
for i in range(num_points):
    x = np.random.normal(0,0.052)
    y = 5*x + np.random.normal(0,0.0201)
    vector.append([x,y])
x_data = [v[0] for v in vector]
y_data = [v[1] for v in vector]

"""
for i in vector_set:
   x_data.append(i[0])
for i in vector_set:
   Y_data.append(i[1])
"""
plt.scatter(x_data, y_data, c='g')



#Define Model
W = tf.Variable(tf.random_uniform([1],-1,1))
b = tf.Variable(tf.zeros([1]))
y = W*x_data + b

#Cost function
loss = tf.reduce_mean(tf.square(y-y_data))

#Optimizer
optimizer = tf.train.GradientDescentOptimizer(0.02)

#Train
train = optimizer.minimize(loss)

#Initialize variables
init_vs = tf.global_variables_initializer()

sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))

with tf.Session() as sess:
    sess.run(init_vs)
    print('W =',sess.run(W),'b =',sess.run(b),'loss =',sess.run(loss))
    for i in range(50000):
        sess.run(train)
        if i%500 == 0:
            print('W =',sess.run(W),'b =',sess.run(b),'loss =',sess.run(loss))
    plt.scatter(x_data, y_data, c='g')
    plt.plot(x_data, sess.run(W)*x_data+sess.run(b), c='r') 