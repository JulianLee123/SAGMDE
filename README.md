# Simulated Annealing Clustering

Code base for the paper ["A simulated annealing algorithm with a dual perturbation method for clustering"](https://www.sciencedirect.com/science/article/abs/pii/S0031320320305161#preview-section-cited-by) authored by Julian Lee and Professor David Perkins and published in *Pattern Recognition*, which currently has 24 citations.

Clustering is a powerful form of unsupervised learning in exploratory data analysis that partitions a set of objects into clusters with the goal of maximizing the similarity of objects within each cluster. Use cases include fraud detection and recommender systems. We propose a new SA-based clustering algorithm that uses two perturbation methods to allow for both large and small shifts in solutions, leading to improved cluster quality. 

## Description

This code base provides an implementation of our proposed simulated annealing-based clustering algorithm SAGMDE (Simulated Annealing with Gaussian Mutation and Distortion Equalization). Existing high-accuracy clustering algorithms were implemented for comparison, including K-means, [SAKM](https://www.worldscientific.com/doi/abs/10.1142/S0218001401000927), and [SAGM](https://www.semanticscholar.org/paper/A-Simulated-Annealing-Clustering-Algorithm-Based-On-Merendino-Celebi/e8ef335803b287ac5c36d1dd2b5afa035dafdf43). In the paper, algorithms are evaluated using common clustering metrics (Normalized Mutual Information and Sum of Squared Errors) across various real and synthetic data sets.

## Usage

The project requires the Java SE Runtime Environment (v. 1.8.0_371 and 1.8.0_161 have been tested). The public static void runner(int dataset, int algorithm, int num_it) method of AllKM.java runs the selected algorithm on the selected dataset and repeats the process num_it times. Results are outputted on the terminal, which can then be copied into a file and processed by the ToCSV class to produce a CSV results file.  After compiling all .java files, the following demo can be run to generate results.

```
java AllKM
java ToCSV sample_output 10
```


