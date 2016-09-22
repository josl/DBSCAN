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


class DBSCAN():
    def __init__(self, file_name, eps=0.1, minPTS=2):
        self.data = pickle.load(open(file_name, 'rb'), encoding='latin1')
        self.spatial_data = cKDTree(self.data.toarray())
        self.data = eps
        self.minPTS = minPTS
        self.clusters = {}
        self.visited = {}

    def query_point_region(self, point, distance):
        pass

    def expand_cluster(self, point, point_neigh, cluster):
        self.clusters[cluster] = [point]
        for new_point in point_neigh:
            visited = [True for i in self.visited if self.visited[i] == point]
            if True in visited:
                self.visited[len(self.visited)] = point
                # Change this to own jaccard distance
                # neigh = self.query_point_region(point, self.eps)
                neigh = self.spatial_data.query_ball_point(new_point, self.eps)
                print(neigh)
                if len(neigh) >= self.minPTS:
                    point_neigh = point_neigh + neigh
            else:
                self.clusters[cluster].append(new_point)

    def start(self):
        cluster = 0
        for point in self.spatial_data.data:
            print(point)
            # print(self.visited)
            visited = [True for i in self.visited if self.visited[i] == point]
            if True in visited:
                self.visited[len(self.visited)] = point
                # Change this to query_point_region
                # neigh = self.query_point_region(point, self.eps)
                neigh = self.spatial_data.query_ball_point(point, self.eps)
                print(neigh)
                if len(neigh) < self.minPTS:
                    continue
                else:
                    # New cluster
                    cluster += 1
                    self.clusters[cluster] = {}
                    self.expand_cluster(P, neigh, C)

                print(neig)
