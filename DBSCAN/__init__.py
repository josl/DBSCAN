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
from sklearn.metrics import jaccard_similarity_score
from scipy.sparse import csr_matrix


# class Coefficient():
#     coefficients = set()

coefficients = set()
class HashPermutation():
    global coefficients
    def __init__(self, N, p=None):
        self.p = p
        self.a, self.b = self.get_coefficients()
        # self.a = a if a is not None else 1
        # self.b = b if b is not None else 1
        self.N = N

    def get_coefficients(self):
        # print(self.p)
        ab = np.random.randint(1, self.p, size=1)[0], np.random.randint(0, self.p, size=1)[0]
        # print(ab)
        # ab = tuple(np.random.randint(1, self.p, size=2))
        return ab
        # print(ab, self.p)
        # if ab in coefficients:
        #     # print(ab)
        #     return self.get_coefficients()
        # else:
        #     coefficients.add(ab)
        #     return ab


    def primes(self, n):
        # https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
        # simple Sieve of Eratosthenes
        # http://stackoverflow.com/a/19498432
        odds = range(3, n+1, 2)
        sieve = set(sum([list(range(q*q, n+1, q+q)) for q in odds], []))
        return [2] + [p for p in odds if p not in sieve]

    def random_prime(n):
        # print('----------------------')
        for random_n in range(n, n * 10):
            # print(random_n)
            # if self.__is_prime(random_n):
            if random_n <= 1:
                continue
            else:
                for i in range(2, int(np.sqrt(random_n)) + 1, 2):
                    if random_n % i == 0:
                        break
                if random_n % i == 0:
                    continue
                return random_n
            # if random_n    in primes:
            #     return random_n

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
        return (((self.a * x) + self.b) % self.p) % self.N
        # return (((self.a * x) + self.b) % self.p)


class LSH():
    # Reference: http://www.mmds.org/mmds/v2.1/ch03-lsh.pdf

    def __init__(self, dist, sparse_matrix, k_permutations):
        global coefficients
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
            np.array([math.inf for i in range(0, self.k_permutations)])
            for j in range(0, self.n_points)
        ]
        self.neigbors = defaultdict(set)
        self.hash_permutations = []
        # a = HashPermutation(self.dimensions, 1, 1, 5)
        # b = HashPermutation(self.dimensions, 3, 1, 5)
        # self.hash_permutations = [a.hash, b.hash]
        # print(self.dimensions)
        # print(self.n_points)
        p = HashPermutation.random_prime(self.n_points)
        print(p)
        # coefficients = Coefficient()
        for k in range(0, self.k_permutations):
            perm = HashPermutation(self.n_points, p)
            # print('-------')
            # coefficients = set()
            self.hash_permutations.append(perm.hash)
        # print(coefficients.coefficients)
        self.permutations = {}

    def signature_distance(self, a, b):
        # intersect = np.intersect1d(a,b)
        intersect = (a == b).sum()
        return intersect / self.k_permutations
        # if intersect != 0:
        #     return intersect / self.k_permutations
        #     # return (self.k_permutations - np.intersect1d(a,b)[0])/self.k_permutations
        # else:
        #     return 0
        # intersect = 0
        # for index in zip(a, b):
        #     if a == b:
        #         intersect += 1
        # return intersect / self.k_permutations

    def jaccard_distance(self, a, b):
        intersect = 0
        union = len(a)
        intersect = (a != b).sum()
        for a_i, b_i in zip(a, b):
            if a_i == 0 and b_i == 0:
                union -= 1
        return (intersect / union)

    def createLHS(self):
        for col_j in range(0, self.n_points):
            for sign_i, hash_func in enumerate(self.hash_permutations):
                # We assume that the matrix has been transformed into an array
                # of indexes
                for row_i in self.point_set[col_j]:
                    # print(row_i)
                    # print(col_j)
                    # print(sign_i)
                    hash_row = hash_func(row_i)
                    # print(hash_row , self.signatures[col_j][sign_i], hash_row < self.signatures[col_j][sign_i])
                    if hash_row < self.signatures[col_j][sign_i]:
                        self.signatures[col_j][sign_i] = hash_row
                    # if hash_func(row_i)[0] < self.signatures[col_j][sign_i]:
                    #     self.signatures[col_j][sign_i] = hash_func(row_i)[0]
    def query_point_region(self):
        for index_a, point_a in enumerate(self.signatures):
            for index_b, point_b in enumerate(self.signatures):
                dist = 1 - self.signature_distance(point_a, point_b)
                # print(point_a, point_b, dist)
                if dist <= self.dist:
                    self.neigbors[index_a].add(index_b)

    def brute_force_distance(self):
        for index_a, point_a in enumerate(self.sparse_matrix.toarray()):
            for index_b, point_b in enumerate(self.sparse_matrix.toarray()):
                dist = self.jaccard_distance(point_a, point_b)
                if dist <= self.dist:
                    self.neigbors[index_a].add(index_b)


class DBSCAN():

    def __init__(self, file_name, eps=0.1, minPTS=2, hash_k=100):
        self.sparse_matrix = pickle.load(
            open(file_name, 'rb'), encoding='latin1')
        self.lsh = LSH(eps, self.sparse_matrix, hash_k)
        self.lsh.createLHS()
        self.lsh.query_point_region()
        self.neigbors = self.lsh.neigbors
        self.dist = eps
        self.minPTS = minPTS
        self.clusters = {}
        self.visited = set()
        self.clustered = set()
        self.noise = set()

    def expand_cluster(self, point, point_neigh, cluster):

        self.clusters[cluster] = set([point])
        # print('init len!', len(point_neigh))
        list_neigh = list(point_neigh)
        tem_neigh = point_neigh
        for new_point in list_neigh:
            # print('new new point', new_point)
            if new_point not in self.visited:
                self.visited.add(new_point)
                neighs = self.neigbors[new_point]
                # print('new \'Burb point', new_point, neighs)
                # print(new_point, neighs)
                if len(neighs) >= self.minPTS:
                    # print('we expand!!!', point_neigh, neighs)
                    # point_neigh = point_neigh.union(neighs)
                    # point_neigh += list(neighs)
                    # print(type(neighs),type(tem_neigh))
                    tem_neigh.update(neighs - tem_neigh)
                    list_neigh += list(neighs - tem_neigh)
                    # list_neigh += list(neighs - set(list_neigh))

                    # print(point_neigh)
            if new_point not in self.clustered:
                self.clustered.add(new_point)
                self.clusters[cluster].add(new_point)
        # print('finish len!', len(point_neigh))

    def start(self):
        cluster = 0
        for point in self.neigbors:
            # for point in range(0, self.lsh.n_points):
            if point not in self.visited:
                self.visited.add(point)
                neigh = self.neigbors[point]
                # print('new_point!', point, neigh)
                if len(neigh) >= self.minPTS:
                    self.expand_cluster(point, neigh, cluster)
                    cluster += 1
                else:
                    self.noise.add(point)
        if len(self.noise) > 0:
            self.clusters[len(self.clusters)] = self.noise
