array = [1,2,4,5,6,7,555,4698]
array.sort()
INPUT = array
sum = 0
for i in range(int(len(INPUT)/2)):
    num_i = INPUT[2*i]
    sum += num_i
print(sum)