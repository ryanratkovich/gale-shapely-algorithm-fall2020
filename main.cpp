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
#include <algorithm>

using namespace std;

/* 	HELPER FUNCTIONS */

double compute_euclidean_distance(const pair<int, int> & p1, const pair<int, int> & p2){ // computes Euclidean distance of 2 points
	return sqrt(pow(p2.first - p1.first,2) + pow(p2.second - p1.second,2));
}

vector<int> compute_nn_edges(vector<pair<int, int>> points){	// computes nearest neighbor of each point
	vector<int> G1(points.size());
	double distance = INT_MAX;
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

bool is_admissable(pair<int, int> p1, pair<int, int> p2, pair<int, int> p3){	// calculates if an sl-edge is admissable using Law of Cosines
	double p1_p2 = compute_euclidean_distance(p1, p2);
	double p3_p2 = compute_euclidean_distance(p3, p2);
	double p1_p3 = compute_euclidean_distance(p1, p3);

	//cout << "Distances for: " << "(" << p1.first << "," << p1.second << ")" << "(" << p2.first << "," << p2.second << ")" << "(" << p3.first << "," << p3.second << ")" << endl;
	//cout << "Length p1->p2: " << p1_p2 << endl;
	//cout << "Length p3->p2: " << p3_p2 << endl;
	//cout << "Length p1->p3: " << p1_p3 << endl;

	double deg = acos(((pow(p3_p2,2) + pow(p1_p2,2) - pow(p1_p3,2)) / (2*p3_p2*p1_p2))) * (180 / M_PI);

	//cout << "deg: " << deg << endl;
	//cout << endl;
	if (deg <= 90 || deg >= 270)
		return false;
	return true;
}

vector<vector<pair<int, int>>> compute_weights(vector<pair<int, int>> sl_vertices_pairs, vector<int> G1, vector<pair<int, int>> points){	// computes weights for each sl-edge
	vector<vector<pair<int ,int>>> preference_list(sl_vertices_pairs.size());	// list of weights for each sl-edge
	for (int i = 0; i < preference_list.size(); ++i){
		double weight = INT_MAX;
		pair<int, int> incident_point_of_i; // get incident point to point "i"
		for (int j = 0; j < G1.size(); ++j){
				if (i != j){
					if (G1[i] == j){
						incident_point_of_i = points[j];
				}
			}
		}
		for (int j = 0; j < sl_vertices_pairs.size(); ++j){
			if (i != j){
				if (sl_vertices_pairs[j] == incident_point_of_i){
					weight = compute_euclidean_distance(sl_vertices_pairs[i], sl_vertices_pairs[j]);
					preference_list[i].push_back(pair<int, int>(weight, G1[j]));
					continue;
				}
				if (is_admissable(incident_point_of_i, sl_vertices_pairs[i], sl_vertices_pairs[j])){
					weight = compute_euclidean_distance(sl_vertices_pairs[i], sl_vertices_pairs[j]);
				} else {
					weight = INT_MAX;
				}
				preference_list[i].push_back(pair<int, int>(weight, G1[j]));
			}
		}
	}

	for (int i = 0; i < preference_list.size(); ++i)
		sort(preference_list[i].begin(), preference_list[i].end());

	return preference_list;
}

vector<int> compute_matching(vector<pair<int, int>> points, vector<int> G1){
	vector<int> G2(G1.size());	// create G2 vector which will contain the edges that complete the reconstruction
	for (int & i : G2)	// initialize all G2 values to -1
		i = -1;

	/* COMPUTE SL-VERTICES */
	vector<int> sl_vertices = compute_sl_vertices(G1);	// stores sl-vertices in a vector

	// cout << "sl-vertices:" << endl;
	// for (int i = 0; i < sl_vertices.size(); ++i){
	// 	cout << "[" << i << "]: " << sl_vertices[i] << endl;
	// }
	// cout << endl;

	vector<pair<int,int>> sl_vertices_pairs;	// vector that holds std::pair's of sl-vertices
	for (const int & i : sl_vertices)
		sl_vertices_pairs.push_back(points[i]);

	// cout << "sl-vertices (pairs):" << endl;
	// for (const auto & p : sl_vertices_pairs)
	// 	cout << p.first << ", " << p.second << endl;
	// cout << endl;

	/* COMPUTE WEIGHTS OF SL-EDGES */
	vector<vector<pair<int, int>>> preference_list = compute_weights(sl_vertices_pairs, G1, points);	// stores weights in a matrix

	cout << "preference_list:" << endl;
	for (int i = 0; i < preference_list.size(); ++i){
		cout << "sl[" << i << "]-(" << sl_vertices_pairs[i].first << "," << sl_vertices_pairs[i].second << "): ";
		for (int j = 0; j < preference_list[i].size(); ++j){
			cout << "(w:" << preference_list[i][j].first << ",G1[" << preference_list[i][j].second << "]) ";
		}
		cout << endl;
	}
	cout << endl;

	/* PERFORM MODIFIED GALE-SHAPELY ALGORITHM */

	vector<vector<int>> proposal_table(G1.size());	// keeps track of which points have proposed to other points
	int i, j = 0;
	int point = G1[0];	// let this be the first point to begin proposing to other points
	bool free_point_available = true;

	while (free_point_available){
		int proposed_point = preference_list[i][j].first;	// point to propose to
		++proposal_table[proposed_point];	// mark that this point has been proposed to by the proposing point
		G2[i] = proposed_point;
		G2[proposed_point] = i;


		for(int k = 0; k < G1.size(); ++k) {
            if(manCurrentMatch[x] == 0) {
                m = x;
                free_point_available = true;
                break;
         	}
       	}
	}

	return G2;
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

	// cout << "Points:" << endl;
	// for (const auto & p : points)
	// 	cout << p.first << "," << p.second << endl;
	// cout << endl;

	vector<int> G1 = compute_nn_edges(points);	// store nearest neighbors

	// cout << "G1 (nearest neighbors):" << endl;
	// for (int i = 0; i < G1.size(); ++i)
	// 	cout << "[" << i << "]: " << G1[i] << endl;
	// cout << endl;

	vector<int> G2 = compute_matching(points, G1); // store G-S algorithm edges

	// cout << "G2 (G-S points):" << endl;
	// for (int i = 0; i < G2.size(); ++i)
	// 	cout << "[" << i << "]: " << G2[i] << endl;
	// cout << endl;

	return 0;
}