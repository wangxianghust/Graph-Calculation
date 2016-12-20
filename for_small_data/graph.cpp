// graph.cpp : kernel start point.
//

#include "stdafx.h"
#include "graph.h"

#include <string>
#include <fstream>
#include <algorithm>
#include <set>
#include <algorithm>
#include <map>
#include <utility>

using std::string;
using std::ifstream;
using std::ofstream;
using std::set;
using std::map;
using std::pair;

// customed split function for string process
vector<string> split(const string &s, const string &delim) {
	vector<string> res;
	string::size_type front = 0;
	string::size_type last = s.find_first_of(delim, front);
	while (last != string::npos) {
		if (last > front) {
			string tmp = s.substr(front, last - front);
			res.push_back(tmp);
		}
		front = last + 1;
		last = s.find_first_of(delim, front);
	}
	if (last > front) {
		res.push_back(s.substr(front, last - front));
	}
	return res;
}

/*
	Input file required: 
	from(int) \t	to(int)
	from..	  \t    to..
*/
auto formatData(string filePath) {
	vector<int> edges; 
	//edges for saved the edges information, 
    //edge0: 0,1 edge1:2,3 etc..the adjacent element forms a edge
	ifstream input(filePath);
	if (input) {
		string line;
		while ( getline(input, line) ) {
			auto from = split(line, "\t")[0];
			edges.push_back(stoi(from));
			auto to = split(line, "\t")[1];
			edges.push_back(stoi(to));
		}
	}
	return edges;
}

Graph constuctGraph(vector<int> edges) {
	Graph G;
	auto max = max_element(edges.begin(), edges.end());
	auto vertexNum = *max + 1;
	vector<VNode> v(vertexNum);
	// add graph data
	for (int i = 0; i < vertexNum; ++i) {
		v[i].data = i;
		G.InsertVertex(v[i]);
	}
	for (int i = 0, j = 1; j < edges.size(); i += 2, j += 2) {
		G.AddEdge(edges[i], edges[j]);
		G.AddEdge(edges[j], edges[i]);  // for undirected-graph
	}
	return G;
}

/*
	degree centrality
*/
vector<double> degreeCentrality(Graph G) {
	vector<double> degree;
	for (int i = 0; i < G.vexNum; ++i) {
		auto d = G.Vertex[i].e.size();
		degree.push_back(d);
		cout << "Node " << i << " degree is " << d << endl;
	}
	auto max = *max_element(degree.begin(), degree.end());
	double centrality = 0;
	for (int i : degree) {
		centrality += max - i;
	}
	centrality = centrality / (G.vexNum * (G.vexNum - 1));
	cout << "centrality is : " << centrality << endl;
	return degree;
}

/*//////////////////////////////////////////
betweenness centerality compute codes strat
*///////////////////////////////////////////
// print the paths of selected two nodes
void printShortestPath(vector<deque<int>> &paths) {
	for (auto &path : paths) {
		cout << endl;
		for (auto &node : path) {
			cout << node << "-->";
		}
	}
}

// print the all paths of the graph(all two non-equal nodes)
void printShortestPath(vector<vector<deque<int>>> &allpath) {
	for (auto &paths : allpath) {
		cout << "all path of two new nodes" << endl;
		for (auto &path : paths) {
			for (auto &node : path) {
				cout << node << "-->";
			}
			cout << endl;
		}
	}
}

/*	
	param: Graph graph
	return:	get all shortest paths of the graph.
		vector<paths>: paths: 0--1, 0--2,...0--5, 1--0,1--2,1--3...1--5...5--0...5--4
	note£ºs-->t and t-->s are all computed.
*/
vector<vector<deque<int>>> getAllShortestPaths(Graph graph) {
	//get all nodes
	auto &nodes = graph.Vertex;
	vector<vector<deque<int>>> allShortestPaths; // store all the shortes path of any two nodes
	for (auto nodeStart = nodes.begin(); nodeStart != nodes.end(); ++nodeStart) {
		for (auto nodeTarget = nodes.begin(); nodeTarget != nodes.end(); ++nodeTarget) {
			if (nodeStart != nodeTarget) {   // ignor the same start and target
				auto paths = graph.shortestPaths(*nodeStart, *nodeTarget);
				allShortestPaths.push_back(paths);
			}
		}
	}
	return allShortestPaths;
}

/*
TEST: print the parents of every node.
*/
void testParents(Graph graph) {
	for (auto &v : graph.Vertex) {
		cout << v.data;
		auto &p = v.getParents();
		for (auto i : p) {
			cout << " ---> " << i.data ;
		}
		cout << endl;
	}
}
/*
 param: nodeNum, the number of vertex
 param: allpath, the all shortest path of the graph
 return: void
 notice: bc_tmp, calculate all internal nodes's bc of two nodes
		 bc_result, accumulate all nodes pair bc_tmp
		 why ? Because using the index of the vector as the id of the node, it is just like a map but cost less than it.
*/
vector<double> betweennessCentrality(int nodeNum, vector<vector<deque<int>>> &allpath){
	vector<double> bc_result(nodeNum, 0);
	vector<int> bc_tmp(nodeNum, 0);
	for (auto &paths : allpath) {	  // calculate the intermediate-bc(for the bc_result) of the two nodes
		float g_st = paths.size();   // float type is for the next calculation precision
		for (auto &i : bc_tmp) {    // the new two nodes, initilize the bc_tmp.
			i = 0;
		}
		for (auto &path : paths) {     
			// the first and the last node on the path should not be considered, just ignore it.
			for (auto begin = ++path.begin(); begin != --path.end(); ++begin) {
				++bc_tmp[*begin];
			}
		}
		for (auto &path : paths) {
			for (auto begin = ++path.begin(); begin != --path.end(); ++begin) {
				//cout << "id: " << *begin << "\t g_st: " << g_st << "\t bc_tmp: " << bc_tmp[*begin];
				bc_result[*begin] += bc_tmp[*begin] / g_st;
				//cout << "\t acc_result: " << bc_result[*begin] << endl;
			}
		}
	}

	for (auto ret : bc_result) {
		cout << ret << endl;
	}
	return bc_result;
}

/*//////////////////////////////////////////
betweenness centerality compute codes end
*///////////////////////////////////////////

/*//////////////////////////////////////////
closeness centerality compute codes start
*///////////////////////////////////////////
/*
	param: nodeNum
	allpath: all shortest paths of the graph
	return: void
	notice: because the ordered saving of the path, every continuous numbers(nodeNum-1) elements own the same start node.
*/
vector<double> closenessCentrality(int nodeNum, vector<vector<deque<int>>> &allpath) {
	int allpathNum = allpath.size();
	int step = nodeNum - 1;
	vector<double> allpathLength(nodeNum, 0);   // as the bc, use the vector index to indexed the node id. 
	for (int i = 0, nodeID = 0; i < allpathNum; i += step, ++nodeID) {
		for (int j = i; j < i + step && j < allpathNum; ++j) {
			auto path = allpath[j].begin();
			int pathLength = (*path).size() - 1; // the path is include the start node , ignore it.
			allpathLength[nodeID] += pathLength;
		}
	}
	for (auto &i : allpathLength) {
		i /= (nodeNum - 1);    // use the number(n-1) to normalized the value.
		i = 1 / i;  // use the iverse to denote it, the average distance is smaller, the more important.
		cout << i << endl;
	}
	return allpathLength;
}


/*//////////////////////////////////////////
closeness centerality compute codes start
*///////////////////////////////////////////



/*
	a interface to output the index value.
*/
void formatOutput(vector<double>& index_value, const string& filePath) {
	ofstream output(filePath);
	auto size = index_value.size();
	for (int i = 0; i < size; ++i) {
		output << i << "\t" << index_value[i] << endl;
	}
}

void formatOutput(map<int, double>& index_value, const string& filePath) {
	ofstream output(filePath);
	for (auto ele : index_value) {
		output << ele.first << "\t" << ele.second << endl;
	}
}

void formatOutput(vector<pair<int, double>>& sorted_index, const string& filePath) {
	ofstream output(filePath);
	for (auto ele : sorted_index) {
		output << ele.first << "\t" << ele.second << endl;
	}
}


/*//////////////////////////////////////////
semi local centerality compute codes start
*///////////////////////////////////////////

void getLinkNodeId(Graph& graph, set<int>& s_id) {  // all node in the old set walk one more step,and get the new id set.
	set<int> s(s_id);								// just modify the parameter use reference.
	for (auto n : s) {
		auto edges = graph.Vertex[n].e;
		for (auto linkNode : edges) {
			s_id.insert(linkNode); 
		}
	}
}

void SCL_walk(Graph& graph, set<int>& nodeIdSet, int stepLength) {
	for (int i = 0; i < stepLength; ++i) {
		getLinkNodeId(graph, nodeIdSet);
	}
}

double semiLocalCentralityNode(Graph& graph,int i) {
	int stepLength = 3;
	set<int> nodeIdSet;
	auto edges = graph.Vertex[i].e;
	for (auto n : edges) {          // Insert first layer of the node i, because the graph.e structure is different from the set<int>
		nodeIdSet.insert(n);        // so the first step we should walk by ourself without call getLinkNodeId.
	}
	SCL_walk(graph, nodeIdSet, stepLength);		// walk three more steps 
	double result = nodeIdSet.size() - 1 ; // ignore the self-node(because it is a undirected graph,the self-node must be in the set)
	return result;
}

vector<double> semiLocalCentrality(Graph& graph) {
	vector<double> sc;
	for (int i = 0; i < graph.vexNum; ++i) {
		auto scn = semiLocalCentralityNode(graph,i);
		sc.push_back(scn);
		cout << scn << endl;
	}
	return sc;
}

/*//////////////////////////////////////////
semi local centerality compute codes end
*///////////////////////////////////////////


/*//////////////////////////////////////////
Eccentricity centerality compute codes start
*///////////////////////////////////////////
double eccentricityNode(int nodeNum, int nodeId, vector<vector<deque<int>>> & allpaths) {
	int allpathNum = allpaths.size();
	int step = nodeNum - 1;
	double m;
	vector<double> tmp(nodeNum, 0);
	for (int i = nodeId * step; i < (nodeId + 1)* step; ++i) {
		auto path = allpaths[i].begin();
		double distance = (*path).size();
		tmp.push_back(distance);
	}
	auto ecc = *max_element(tmp.begin(), tmp.end()) - 1; // path lenth = nodeNumInPath - 1
	ecc = 1 / ecc; // ecc denotes the max_shortestpath of all other nodes, we use inverse to express it's importance.
	return ecc;
}

vector<double> Eccentricity(int nodeNum,vector<vector<deque<int>>> & allpaths) {
	vector <double> result;
	for (int i = 0; i < nodeNum; ++i) {
		auto ret = eccentricityNode(nodeNum, i,allpaths);
		result.push_back(ret);
		cout << ret << endl;
	}
	return result;
}

/*//////////////////////////////////////////
Eccentricity centerality compute codes end
*///////////////////////////////////////////


/*//////////////////////////////////////////
//A function to rank the vector number.
*///////////////////////////////////////////
// a compare function for sort, then we can use it the sort the map by value 
struct CmpByValue {
	bool operator()(const pair<int, double> & lhs, const pair<int, double> & rhs)
	{
		return lhs.second > rhs.second;
	}
};

void importanceRank(vector<double>& index_value, vector<pair<int, double>> &sorted_index) {
	map<int, double> tmp;
	int nodeID = 0;
	for (auto value : index_value) {
		tmp.insert({ nodeID++, value });
		//tmp.insert({ ++nodeID, value });
	}
	for (auto i : tmp) {
		sorted_index.push_back(i);
	}
	sort(sorted_index.begin(), sorted_index.end(), CmpByValue());
}
/////////////////////////////////////////////////////////

void testInfo(const string& info) {
	cout << "-----" << info << "-----" << endl;
}

int main() {
	string filePath("data/test.txt");
	auto data = formatData(filePath);
	Graph graph = constuctGraph(data);
	graph.BreadthFirstSearch();

	testInfo("DEGREE TEST");
	auto degree = degreeCentrality(graph);

	testInfo("Running,please waite...");
	auto allpaths = getAllShortestPaths(graph);

	testInfo("ALL GRAPH PRINT");
	printShortestPath(allpaths);

	testInfo("BC PRINT");
	auto bc = betweennessCentrality(graph.vexNum, allpaths);

	testInfo("CC PTRINT");
	auto cc = closenessCentrality(graph.vexNum, allpaths);

	testInfo("SLC PRINT");
	auto slc = semiLocalCentrality(graph);

	testInfo("ECC PRINT");
	auto ecc = Eccentricity(graph.vexNum, allpaths);

	formatOutput(degree, "./result/result_degree.txt");
	formatOutput(bc, "./result/result_bc.txt");
	formatOutput(cc, "./result/result_cc.txt");
	formatOutput(slc, "./result/result_slc.txt");
	formatOutput(ecc, "./result/result_ecc.txt");

	/*
	vector<pair<int, double>> degree_rank;
	importanceRank(degree, degree_rank);
	formatOutput(degree_rank, "./result/result_degree.txt");

	vector<pair<int, double>> bc_rank;
	importanceRank(bc, bc_rank);
	formatOutput(bc_rank, "./result/result_bc.txt");

	vector<pair<int, double>> cc_rank;
	importanceRank(cc, cc_rank);
	formatOutput(cc_rank, "./result/result_cc.txt");

	vector<pair<int, double>> slc_rank;
	importanceRank(slc, slc_rank);
	formatOutput(slc_rank, "./result/result_slc.txt");

	vector<pair<int, double>> ecc_rank;
	importanceRank(ecc, ecc_rank);
	formatOutput(ecc_rank, "./result/result_ecc.txt");

	*/
}


