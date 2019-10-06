array = [3,3,7,7,10,11,11]
i = 0
if len(array) == 1:
    print(array[i])
else:
    while i <= int(len(array)):
        if array[i] == array [i+1]:
            i += 2
        else:
            print(array[i])
            break