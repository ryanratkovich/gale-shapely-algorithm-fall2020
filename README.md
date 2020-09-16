# COP4531_Project1_Ratkovich

GENERAL:
This is the Gitlab repository for all files related to my solution to project 1 for COP4531 Fall 2020.

NOTE ON COMPILATION AND EXECUTION:
This program can be compiled using the makefile provided within this repository. Type "make" at the commandline.
To run this program: enter "project1 < inputfile.txt" at the commandline.

PROJECT DESCRIPTION:
This program aims to construct a closed, smooth curve from a set of points in a two-dimensional plane. This program accomplishes this by computing each points' nearest neighboring point and by then storing the results of this calculation in a vector. The program then computes which points are "single and looking" (i.e. have only one incident edge on them). Once the sl-vertices are computed, the program calculates each sl-vertices' Euclidean distance to each other sl-vertex and determines if they are "admissable" to have an edge between them. Each sl-vertex is then assigned a "weight" which is based on their Euclidean distance and admissability to be paired with other sl-vertices. Finally, this program performs a modified version of the Gale-Shapely algorithm using the weights of the sl-vertices to calculate the remaining, stable edges that need to be made on the set of points so that a smooth, closed curve can be constructed. The program then outputs this result in the form of a sequence consisting of each point one after another in the counter-clockwise direction.

NOTE ON FUNCTIONS:
This program uses the following functions to accomplish this goal:

- compute_euclidean_distance : computes the distance between any two points on a 2D plane.

- compute_nn_edges : computes the nearest neighboring point of each point.

- compute_sl_vertices : computes which points are "single and looking" (one incident edge on the point)

- is_admissable : using the Law of Cosines, computes of two points are admissable to be paired based on the angle that is made between them on the 2D plane.

- compute_weights: computes each points preference list (a list containing the weights of every other point in ascending order).

- compute_matching: uses date from all other functions to construct a vector containing the points necessary to construct a smooth, closed curve on the set of points given as input. Peforms the Gale-Shapely algorithm on all points.





