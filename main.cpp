/*
* Class: COP4531
* Assignment: 1
* Author: Ryan Ratkovich
* Date: 09/15/2020
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <climits>

using namespace std;

/* 	HELPER FUNCTIONS */

int compute_euclidean_distance(const pair<int, int> & p1, const pair<int, int> & p2){ // computes Euclidean distance of 2 points
	return sqrt(pow(p2.first - p1.first,2) + pow(p2.second - p1.second,2));
}

vector<int> compute_nn_edges(vector<pair<int, int>> points){	// computes nearest neighbor of each point
	vector<int> G1(points.size());
	int distance = INT_MAX;
	int nn = -1;	// holds position of point's nearest neighbor
	for (int i = 0; i < points.size(); ++i){
		for (int j = 0; j < points.size(); ++j){
			if (i != j){
				if (compute_euclidean_distance(points[i], points[j]) <= distance){
					distance = compute_euclidean_distance(points[i], points[j]);
					nn = j;
				}
			}
		}
		G1[i] = nn;	// assign nearest neighbor in G1
		nn = -1;	// reset nn
	}
	return G1;
}

vector<int> compute_sl_vertices(vector<int> G1){	// computes the sl-vertices

	vector<int> deg_v(G1.size());	// vector that holds degree of each vertex
	int degree = 0;
	for (int i = 0; i < deg_v.size(); ++i){
		if (G1[i] != -1)	//	if point "i" is connected to any other point, increment deg(i)
			++degree;
		for (int j = 0; j < G1.size(); ++j){
			if (G1[j] == i)	// if point "j" is connected to "i", increment deg(i)
				++degree;
		}		
		deg_v[i] = degree;
		degree = 0;
	}
	for (int i = 0; i < deg_v.size(); ++i){
		if (deg_v[i] == 2){
			for (int j = 0; j < deg_v.size(); ++j){
				if (i == G1[j] && j == G1[i])
					--deg_v[i];
			}
		}
	}

	vector<int> sl_vertices;	// vector that holds sl-vertices
	for (int i = 0; i < deg_v.size(); ++i){
		if (deg_v[i] == 1)
			sl_vertices.push_back(i);
	}
	return sl_vertices;
}

vector<int> compute_matching(vector<pair<int, int>> points, vector<int> G1){

	vector<int> G2(G1.size());	// create G2 vector which will contain the edges that complete the reconstruction
	for (int & i : G2)	// initialize all G2 values to -1
		i = -1;

	/* COMPUTE SL-VERTICES */

	vector<int> sl_vertices = compute_sl_vertices(G1);

	vector<pair<int,int>> sl_vertices_pairs;	// vector that holds std::pair's of sl-vertices
	for (const int & i : sl_vertices)
		sl_vertices_pairs.push_back(points[i]);

	/* COMPUTE WEIGHTS OF SL-EDGES */

	return sl_vertices;
}

bool is_admissable(pair<int, int> p1, pair<int, int> p2){

}

/* MAIN */

int main() {

	ifstream file;	// open input file
	file.open("inputfile.txt");

	int numLines = 0;	
	cin >> numLines;	// store number of lines

	string s_arr[numLines];				// store string representation of points from file
	for (int i = 0; i < numLines; ++i)
		cin >> s_arr[i];
	file.close();

	vector<pair<int, int> > points;		// store the points as std::pair in a vector
	for (const string & s : s_arr){
		int first = s[0] - 48;			// convert ascii char's to correct int
		int second = s[2] - 48;
		points.push_back(pair<int, int>(first, second));
	}

	cout << "Points:" << endl;
	for (const auto & p : points)
		cout << p.first << "," << p.second << endl;
	cout << endl;

	vector<int> G1 = compute_nn_edges(points);	// store nearest neighbors

	cout << "G1 (nearest neighbors):" << endl;
	for (int i = 0; i < G1.size(); ++i)
		cout << "[" << i << "]: " << G1[i] << endl;
	cout << endl;

	vector<int> sl_vertices = compute_matching(points, G1); // store G-S algorithm edges

	cout << "sl_vertices:" << endl;
	for (int i = 0; i < sl_vertices.size(); ++i)
		cout << "[" << i << "]: " << sl_vertices[i] << endl;
	cout << endl;

	return 0;
}