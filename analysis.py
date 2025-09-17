import sklearn.metrics as sk
import math
import numpy as np
import pandas as pd
import sys
import symnmfmodule as sm

epsilon = 10 ** -4

# Initialization of the matrix H for symNMF algorithm
def H_init(n, k, W):
    m = np.mean(W)
    np.random.seed(0)
    upper_bound = 2 * math.sqrt(m/k)
    H = np.random.uniform(0,upper_bound,size=(n,k))
    return H

# Calculate euclidean distance
def euclidean_distance(p, q, d):
    distance = 0.0
    for i in range(d):
        distance += math.pow((float(p[i]) - float(q[i])), 2)
    return math.sqrt(distance)

# Read the input file data into numpy array
def read_file(file_name):
    data = pd.read_csv(file_name, header=None).values
    data_array = np.array(data, dtype=float)
    return data_array

# Kmeans algorithm from HW1
def kmeans(k, data):
    d = len(data[0])  # points dimension
    N = len(data)  # number of points
    centroids = data[:k]
    labels = np.zeros(N)
    for i in range(300):
        save_vector = [[0.0 for r in range(d)] for j in range(k)]
        save_cnt = [0 for t in range(k)]
        for idx, xi in enumerate(data):
            dis = [euclidean_distance(xi, centroids[j], d) for j in range(k)]
            min_val = min(dis)
            min_index = dis.index(min_val)
            labels[idx] = min_index # Assign point to cluster
            save_cnt[min_index] += 1
            for l in range(d):
                save_vector[min_index][l] += xi[l]
        for j in range(k):
            if save_cnt[j] == 0:
                save_cnt[j] = 1
            for m in range(d):
                save_vector[j][m] /= save_cnt[j]
        flag = True
        for g in range(k):
            delta = euclidean_distance(save_vector[g], centroids[g], d)
            if delta >= epsilon:
                flag = False
                break
        centroids = save_vector
        if flag:
            break
    return labels

# SymNMF algorithm
def symnmf(N,K,data):
        data = data.tolist()
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
        labels = np.argmax(sol, axis=1) # Assign every point to cluster
        return labels

def main():
    K = int(float(sys.argv[1]))
    file_name = sys.argv[2]
    data = read_file(file_name)
    N = len(data)   # number of points
    symnmf_labels = symnmf(N,K,data)
    synmnmf_score = sk.silhouette_score(data, symnmf_labels)
    kmeans_score = sk.silhouette_score(data, kmeans(K,data))
    print("nmf:", "%.4f" % synmnmf_score)
    print("kmeans:", "%.4f" % kmeans_score)

if __name__ == "__main__":
    main()    