/*
* Class: COP4531
* Assignment: 1
* Author: Ryan Ratkovich
* Date: 09/16/2020
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

bool is_admissable(pair<int, int> p1, pair<int, int> p2, pair<int, int> p3, pair<int, int> p4){	// calculates if an sl-edge is admissable using Law of Cosines
	double p1_p2 = compute_euclidean_distance(p1, p2);
	double p3_p2 = compute_euclidean_distance(p3, p2);
	double p1_p3 = compute_euclidean_distance(p1, p3);

	double deg1 = acos(((pow(p3_p2,2) + pow(p1_p2,2) - pow(p1_p3,2)) / (2*p3_p2*p1_p2))) * (180 / M_PI);

	double p2_p3 = compute_euclidean_distance(p2, p3);
	double p4_p3 = compute_euclidean_distance(p4, p3);
	double p2_p4 = compute_euclidean_distance(p2, p4);

	double deg2 = acos(((pow(p4_p3,2) + pow(p2_p3,2) - pow(p2_p4,2)) / (2*p4_p3*p2_p3))) * (180 / M_PI);

	if (deg1 < 90 && deg2 < 90)	// if the angles between points in a potential edge are small, reject
		return false;
	return true;
}

vector<vector<pair<int, int>>> compute_weights(vector<pair<int, int>> sl_vertices_pairs, vector<int> G1, vector<pair<int, int>> points, vector<int> sl_vertices){	// computes weights for each sl-edge
	vector<vector<pair<int ,int>>> preference_list(sl_vertices_pairs.size());	// list of weights for each sl-edge
	for (int i = 0; i < preference_list.size(); ++i){
		double weight = INT_MAX;	// weight = INFINITY 
		pair<int, int> incident_point_of_i; // get incident point to point "i"
		for (int j = 0; j < G1.size(); ++j){
			if (i != j){
				if (G1[i] == j){
					incident_point_of_i = points[j];
					break;
				}
			}
		}
		for (int j = 0; j < sl_vertices_pairs.size(); ++j){
			if (i != j){
				pair<int, int> incident_point_of_j; // get incident point to point "j"
				for (int k = 0; k < G1.size(); ++k){
					if (j != k){
						if (G1[j] == k){
							incident_point_of_j = points[k];
							break;
						}
					}
				}
				if (sl_vertices_pairs[j] == incident_point_of_i){
					// weight = compute_euclidean_distance(sl_vertices_pairs[i], sl_vertices_pairs[j]);
					// preference_list[i].push_back(pair<int, int>(weight, G1[j]));
					continue;
				}
				if (is_admissable(incident_point_of_i, sl_vertices_pairs[i], sl_vertices_pairs[j], incident_point_of_j)){
					weight = compute_euclidean_distance(sl_vertices_pairs[i], sl_vertices_pairs[j]);
				} else {
					weight = INT_MAX;
				}
				preference_list[i].push_back(pair<int, int>(weight, j));
			}
		}
	}

	for (int i = 0; i < preference_list.size(); ++i)
		sort(preference_list[i].begin(), preference_list[i].end());	//	sort the preference lists

	return preference_list;
}

vector<int> compute_matching(vector<pair<int, int>> points, vector<int> G1){
	vector<int> G2(G1.size());	// create G2 vector which will contain the edges that complete the reconstruction
	for (int & i : G2)	// initialize all G2 values to -1
		i = -1;

	/* COMPUTE SL-VERTICES */
	vector<int> sl_vertices = compute_sl_vertices(G1);	// stores sl-vertices in a vector

	vector<pair<int,int>> sl_vertices_pairs;	// vector that holds std::pair's of sl-vertices
	for (const int & i : sl_vertices)
		sl_vertices_pairs.push_back(points[i]);

	/* COMPUTE WEIGHTS OF SL-EDGES */
	vector<vector<pair<int, int>>> preference_list = compute_weights(sl_vertices_pairs, G1, points, sl_vertices);	// stores weights in a matrix

	/* PERFORM MODIFIED GALE-SHAPELY ALGORITHM */

	vector<vector<int>> proposal_table(preference_list.size(), vector<int>(preference_list.size()));;	// keeps track of which points have proposed to other points

	for (int i = 0; i < proposal_table.size(); ++i){	// initialize proposal table
		for (int j = 0; j < proposal_table[i].size(); ++j){
			if (j == proposal_table[i].size() - 1)
				proposal_table[i][j] = 1;
			if (i == j)
				proposal_table[i][j] = 1;
		}
	}

	int point = 0;	// let this be the first point to begin proposing to other points
	bool free_point_exists = true;
	while (free_point_exists){
		free_point_exists = false;
		int preferred_point_index = 0;
		int preferred_point = 0;	// placeholder for the point to be proposed to
		for (int i = 0; i < preference_list[point].size(); ++i){
			if (proposal_table[point][i] == 0){	// if point hasnt proposed to point "i"
				preferred_point_index = preference_list[point][i].second;	// index of preferred point in sl_vertices
				preferred_point = sl_vertices[preference_list[point][i].second];	// value of preferred point
				proposal_table[point][i] = 1;
				break;
			}
		}
		if (G2[preferred_point] == -1){	// if G2[preferred_point] is free, the points are matched in G2
			G2[point] = preferred_point;
			G2[preferred_point] = point;
		} else {	// else if preferred_point prefers point to their current match
			int preferred_point_current_partner = G2[preferred_point];
			int weight_of_point = 0;
			for (int i = 0; i < preference_list[preferred_point_index].size(); ++i){	// get point's weight in proposed points preference list
				if (preference_list[preferred_point_index][i].second == point){
					weight_of_point = preference_list[preferred_point_index][i].first;
					break;
				}
			}
			for (int i = 0; i < preference_list[preferred_point_index].size(); ++i){
				if (weight_of_point < preference_list[preferred_point_index][i].first){
					G2[point] = preferred_point;
					G2[preferred_point] = point;
					G2[preferred_point_current_partner] = -1;
					break;
				}
			}
		}

		for(int i = point; i < proposal_table.size(); ++i) {	// find another free point, if one exists
			bool has_proposed_to_every_point = true;
			for (int j = 0; j < proposal_table[i].size(); ++j){	// check if point has proposed to every other point
				if (proposal_table[i][j] == 0){
					has_proposed_to_every_point = false;
					break;
				}
			}
            if(G2[i] == -1 && !has_proposed_to_every_point) { // if free and hasnt proposed to all points
                point = i;	// i becomes next point to propose
                free_point_exists = true;
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

	vector<int> G1 = compute_nn_edges(points);	// store nearest neighbors

	vector<int> G2 = compute_matching(points, G1); // store G-S algorithm edges

	int point = 0;	// print out the bipartite matching of G1 and G2
	for (int i = 0; i < G2.size() + 1; ++i){
		cout << point << endl;
		if (i % 2 == 0)
			point = G2[point];
		else
			point = G1[point];
	}

	return 0;
}