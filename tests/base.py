#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of DBSCAN.
# https://github.com/josl/DBSCAN

# Licensed under the MIT license:
# http://www.opensource.org/licenses/MIT-license
# Copyright (c) 2016, Jose L. Bellod Cisneros <bellod.cisneros@gmail.com>

from unittest import TestCase as PythonTestCase
import pickle
from DBSCAN import DBSCAN


class TestCase(PythonTestCase):
    def test_data_10(self):
        file_name = './test_files/data_10points_10dims.dat'
        # data = pickle.load(open(file_name, 'rb'), encoding='latin1')
        # self.assertTrue(data.shape[0] == 10)
        dbscan = DBSCAN(file_name, 0.4, 2)
        dbscan.start()
        print(dbscan.clusters)
        self.assertTrue(len(dbscan.clusters) == 4)

    # def test_data_100(self):
    #     file_name = './test_files/data_100points_100dims.dat'
    #     data = pickle.load(open(file_name, 'rb'), encoding='latin1')
    #     self.assertTrue(data.shape[0] == 100)
    #     # dbscan = new DBSAN(0.4, 2)
    #     # dbscan.start()
    #     # self.assertTrue(dbscan.n_clusters == 6)
    #
    # def test_data_1000(self):
    #     file_name = './test_files/data_1000points_1000dims.dat'
    #     data1000 = pickle.load(open(file_name, 'rb'), encoding='latin1')
    #     self.assertTrue(data1000.shape[0] == 1000)
    #     # dbscan = new DBSAN(0.4, 2)
    #     # dbscan.start()
    #     # self.assertTrue(dbscan.n_clusters == 9)
    #
    # def test_data_10000(self):
    #     file_name = './test_files/data_10000points_10000dims.dat'
    #     data = pickle.load(open(file_name, 'rb'), encoding='latin1')
    #     self.assertTrue(data.shape[0] == 10000)
    #     # dbscan = new DBSAN(0.4, 2)
    #     # dbscan.start()
    #     # self.assertTrue(dbscan.n_clusters == 394)
    #
    # def test_data_100000(self):
    #     file_name = './test_files/data_100000points_100000dims.dat'
    #     data = pickle.load(open(file_name, 'rb'), encoding='latin1')
    #     self.assertTrue(data.shape[0] == 100000)
    #     # dbscan = new DBSAN(0.4, 2)
    #     # dbscan.start()
    #     # self.assertTrue(dbscan.n_clusters == 1692)
