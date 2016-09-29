#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of DBSCAN.
# https://github.com/josl/DBSCAN

# Licensed under the MIT license:
# http://www.opensource.org/licenses/MIT-license
# Copyright (c) 2016, Jose L. Bellod Cisneros <bellod.cisneros@gmail.com>

from DBSCAN.version import __version__  # NOQA
import pickle
# from scipy.sparse.csr_matrix import csr_matrix
from scipy.spatial import cKDTree
import numpy as np
import math
from collections import defaultdict

class HashPermutation():
    def __init__(self, N, a=0, b=0, p=0):
        random = np.random.randint(100, size=2)
        self.a = a if a != 0 else random[0]
        self.b = b if b != 0 else random[1]
        self.p = p if p != 0 else self.random_prime(N)
        self.N = N
    def random_prime(self, N):
        for random_n in np.random.randint(N + 1, N*2, size=100):
            if self.__is_prime(random_n):
                return random_n
    def __is_prime(self, x):
        # TODO look at: https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
        if x <= 1:
            return False
        else:
            for i in range(2, int(np.sqrt(x)) + 1):
                if x % i == 0:
                    return False
        return True
    def hash(self, x):
        # TODO: Hashing a vector????
        return (((self.a * x) + self.b) % self.p) % self.N
class LSH():
    # Reference: http://www.mmds.org/mmds/v2.1/ch03-lsh.pdf
    def __init__(self, dist, sparse_matrix, k_permutations):
        self.dist = dist
        self.sparse_matrix = sparse_matrix
        # sparse_matrix = [
        #     [1, 0, 0, 1],
        #     [0, 0, 1, 0],
        #     [1, 1, 0, 1],
        #     [0, 0, 1, 1],
        #     [0, 0, 1, 0]
        # ]
        # row = np.array([0, 0, 1, 2, 2, 3, 3, 3, 4])
        # col = np.array([0, 3, 2, 1, 3, 0, 2, 3, 2])
        # data = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])
        # self.sparse_matrix = csr_matrix((data, (row, col)), shape=(5, 4)).transpose()
        self.point_set = {}
        for row_i, row in enumerate(self.sparse_matrix):
            self.point_set[row_i] = row.indices
        self.point_dict = {}
        self.k_permutations = k_permutations
        self.dimensions = self.sparse_matrix.shape[1]
        # self.dimensions = sparse_matrix.dim()
        self.n_points = self.sparse_matrix.shape[0]
        # Innit signature matrix: permutations x dimension
        self.signatures = [
            [math.inf for i in range(0, self.k_permutations)]
            for j in range(0, self.n_points)
        ]
        self.neigbors = defaultdict(set)
        self.hash_permutations = []
        # a = HashPermutation(self.dimensions, 1, 1, 5)
        # b = HashPermutation(self.dimensions, 3, 1, 5)
        # self.hash_permutations = [a.hash, b.hash]
        for k in range(0, self.k_permutations):
            perm = HashPermutation(self.dimensions)
            self.hash_permutations.append(perm.hash)
        self.permutations = {}
    def jaccard_distance(self, a, b):
        intersect = 0
        for index in range(0, len(a)):
            if a[index] == b[index]:
                intersect += 1
        return 1 - (intersect / len(a))
        # return 1 - (len(a & b) / len(a | b) )
    def createLHS(self):
        for col_j in range(0, self.n_points):
            for sign_i, hash_func in enumerate(self.hash_permutations):
                # We assume that the matrix has been transformed into an array of indexes
                for row_i in self.point_set[col_j]:
                    if hash_func(row_i) < self.signatures[col_j][sign_i]:
                        self.signatures[col_j][sign_i] = hash_func(row_i)
    def query_point_region(self):
        for index_a, point_a in enumerate(self.signatures):
            for index_b, point_b in enumerate(self.signatures):
                dist = self.jaccard_distance(point_a, point_b)
                print(point_a, point_b, dist)
                if dist <= self.dist:
                    self.neigbors[index_a].add(index_b)


class DBSCAN():
    def __init__(self, file_name, eps=0.1, minPTS=2, hash_k=100):
        self.sparse_matrix = pickle.load(open(file_name, 'rb'), encoding='latin1')
        # self.spatial_data = cKDTree(self.data.toarray())
        self.lsh = LSH(eps, self.sparse_matrix, hash_k)
        self.lsh.createLHS()
        self.lsh.query_point_region()
        self.neigbors = self.lsh.neigbors
        self.dist = eps
        self.minPTS = minPTS
        self.clusters = {}
        self.visited = {}
        for point in self.lsh.point_set:
            self.visited[point] = False
    def query_point_region(self, point, distance):
        pass
    def expand_cluster(self, point, point_neigh, cluster):
        for new_point in point_neigh:
            if not self.visited[point]:
                self.visited[point] = True
                neigh = self.neigbors[new_point]
                print(neigh)
                if len(neigh) >= self.minPTS:
                    point_neigh = self.neigbors[point].add(neigh)
            if new_point not in self.clusters:
                self.clusters[cluster].add(new_point)
    def start(self):
        cluster = 0
        for point in self.lsh.point_set:
            print(point)
            if not self.visited[point]:
                self.visited[point] = True
                neigh = self.neigbors[point]
                print(neigh)
                # New cluster
                cluster += 1
                self.clusters[cluster] = {}
                self.clusters[cluster] = set([point])
                if len(neigh) >= self.minPTS:
                    self.expand_cluster(point, neigh, cluster)
                print(neigh)
