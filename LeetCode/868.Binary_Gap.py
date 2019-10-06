def binary_gap(num):
    print(format(num,"b"))
    N = str(format(num,"b"))
    gap_i = 0
    gap_max = 0
    for i in range(len(N)-1):
        if N[i] == N[i+1] == "1":
            gap_i += 1
            gap_max = max(gap_i,gap_max)
        else:
            gap_i = 0
    if gap_max == 0:
        print('There is no consequtive 1s')
        return("Nothing")
        
    else:
        print(f'The max consequtive number of "1" is: ',gap_max+1)
        return(gap_max+1)

binary_gap(111111111111111111111111111111111111)