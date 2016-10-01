# DBSCAN

Exercise 4.1 of 02807 Computational Tools for Big Data at DTU

Implement the DBSCAN clustering algorithm to work with Jaccard-distance as its metric. It should be able to handle sparse data.

We suggest that you implement the DBSCAN clustering algorithm as described in the Wikipedia article, linked to above. The algorithm should use Jaccard-distance (1 minus Jaccard index) when measuring distance between points. In this exercise two points that are exactly a distance epsilon apart are taken to be in the same cluster, thus greater-than-or-equal dist(a, b) <= epsilon and not dist(a, b) < epsilon should be used.

Your method should support a large number (~100,000) of sparse points in very high-dimensional space (~100,000 dimensions). This means that the points – even though they lie in very high dimensional space – only have few non-zero coordinates.

Debug and test your algorithm with the following:

Use this data. The files are “pickled”, and must first be loaded with pickle. Parse the loaded data to a compressed sparse row matrix using csr_matrix, where each row corresponds to a point and each column to a dimension.
For the five input files, use the following parameters

- File 1: eps=0.4 and M=2
- File 2: eps=0.3 and M=2
- File 3: eps=0.15 and M=2
- File 4: eps=0.15 and M=2
- File 5: eps=0.15 and M=2

You should get this many clusters when running your algorithm on the five test files (also counting the special cluster of points not assigned to a cluster):
4, 6, 9, 394, 1692
Remember to describe how you process the data in order for the algorithm to work with very high dimensional data.

You should also document how many points is contained in the largest cluster for each of the 5 test files. The correct answers for this will be released for giving peer feedback.

It is much more important to have a solution which can only solve the first few input files correctly than a fancy solution solving none of them. There will be partial credit for the exercise even if you can only make a solution for the smallest input file.

## Install

It's recommended to create a virtual environment (conda env preffered)
``` bash
make setup
```
## Usage

## Testing

# Unit test
``` bash
make unit
```


## License
MIT © [Jose Luis Bellod Cisneros](http://josl.github.io)
