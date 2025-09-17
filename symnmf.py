import numpy as np
import pandas as pd
import sys
import math
import symnmfmodule as sm

# Initialization of the matrix H for symNMF algorithm
np.random.seed(0)
def H_init(n, k, W):
    m = np.mean(W)
    upper_bound = 2 * math.sqrt(m/k)
    H = np.random.uniform(0,upper_bound,size=(n,k))
    return H

# Read the input file data into numpy array
def read_file(file_name):
    try:
        data = pd.read_csv(file_name, header=None).values
    except:
        return None
    data_array = np.array(data, dtype=float)
    return data_array

def main():
    try:
        K = int(float(sys.argv[1]))
        goal = sys.argv[2]
        file_name = sys.argv[3]
    except:
        print("An Error Has Occurred")
        return 1
    data = read_file(file_name)
    if(data is None):
        print("An Error Has Occurred")
        return 1
    N = len(data)   # number of points
    data = data.tolist()
    if(goal == "symnmf"):
        W = sm.norm(data)
        if(W is None):
            print("An Error Has Occurred")
            return 1
        H = H_init(N,K,W)
        W = np.array(W)
        H = np.array(H)
        W = W.flatten().tolist()
        H = H.flatten().tolist()
        sol = sm.symnmf(N,K,H,W)
        if(sol is None):
            print("An Error Has Occurred")
            return 1
        sol = np.array(sol)
        sol = sol.reshape(N,K)
    else:
        if (goal == "sym"):
            sol = sm.sym(data)
        elif (goal == "ddg"):
            sol = sm.ddg(data)
        else:
            sol = sm.norm(data)
        if(sol is None):
            print("An Error Has Occurred")
            return 1
        sol = np.array(sol)
        sol = sol.reshape(N,N)
    for vector in sol:
        vector = [("%.4f" % num) for num in vector]
        print(*vector, sep = ",")
    
 
if __name__ == "__main__":
    main()    