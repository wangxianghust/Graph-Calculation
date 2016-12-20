#pragma once
#include <iostream>
#include <vector>
#include <queue>
#include <deque>

using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::queue;
using std::deque;


typedef int DataType;
const int INIT_DATA = -1;
const int NO_NODE = -1;
const int NO_EDGE = -2;

class VNode {
public:
	bool visited;
	int layer;
	int numPaths;
	vector<VNode> parents;
	DataType data;
	vector<int> e; // edge lists of this node

				   //member function
	VNode(DataType val = INIT_DATA, bool flag = false, int layer = 0, int num = 1) : data(val), visited(flag), layer(layer), numPaths(num) {}
	void Visit() {
		cout << data << ", layer is : " << layer;
		cout << " parent is : ";
		for (auto &node : parents) {
			cout << node.data << "---";
		}
		cout << endl;
	}

	bool operator== (VNode &y) {
		if (data == y.data && e == y.e)
			return true;
		else
			return false;
	}

	bool operator!= (VNode &y) {
		if (data != y.data  || e != y.e)
			return true;
		else
			return false;
	}

	vector<VNode> getParents();
};

vector<VNode> VNode::getParents() {
	return parents;
}

class Graph {
public:
	int vexNum;
	int edgeNum;
	vector<VNode> Vertex;

	Graph() {
		vexNum = 0;
		edgeNum = 0;
	}

	Graph(vector<VNode> &v);

	int GetVexNum();
	int GetEdgeNum();

	void InsertVertex(VNode &v);

	bool AddEdge(int x, int y);

	void BFS(VNode &v);
	void BreadthFirstSearch();

	void DFS(VNode &x);     //recursion funcion express

	auto shortestPaths(VNode &v, VNode &u);   // get all the shortest of node u and v.

	~Graph() {};



};


Graph::Graph(vector<VNode>& v)
{
	Vertex.assign(v.begin(), v.end());
	vexNum = Vertex.size();
	edgeNum = GetEdgeNum();
}

int Graph::GetVexNum()
{
	return Vertex.size();
}

int Graph::GetEdgeNum()
{
	edgeNum = 0;
	for (auto iter = Vertex.begin(); iter != Vertex.end(); ++iter) {
		edgeNum += (iter->e).size();
	}
	return edgeNum;
}

void Graph::InsertVertex(VNode & v)
{
	Vertex.push_back(v);
	++vexNum;
}

// add the edge x->y to the graph, no multiple-dege allowed. FOR UNDIRECTED GRAPH, use AddEdge(x,y) and AddEdge(y,x)
bool Graph::AddEdge(int x, int y)
{
	if ((x >= Vertex.size() || x < 0) || (y >= Vertex.size()) || y < 0) {
		cout << "AddEdge problem : the vertex id is out of range !" << endl;
		return false;
	}
	int size = Vertex[x].e.size();
	for (auto iter = Vertex[x].e.begin(); iter != Vertex[x].e.end(); ++iter) {
		if (Vertex[*iter] == Vertex[y]) {
			cout << "The edge " << x << " ---> " << y << "is already existed !" << endl;
			return false;
		}
	}
	Vertex[x].e.push_back(y);
	++edgeNum;
	return true;
}

void Graph::BFS(VNode & v)
{
	queue<VNode> q;
	q.push(v);

	while (!q.empty())
	{
		VNode _v = q.front();
		if (!_v.visited) {
			//_v.Visit();
			_v.visited = true;
			int size = _v.e.size();
			for (int i = 0; i < size; ++i) {
				VNode &w = Vertex[_v.e[i]];
				if (!w.visited) {           // the visited element is not added
					//w.Visit();
					w.visited = true;
					q.push(w);
				}
			}
		}
		q.pop();
	}
}

void Graph::BreadthFirstSearch()
{
	for (int i = 0; i != vexNum; ++i)
	{
		if (!Vertex[i].visited)
			BFS(Vertex[i]);     // use Vertex[0], Vertex[1]...as the parameter of the BFS
	}
	for (auto k = Vertex.begin(); k != Vertex.end(); ++k)
		k->visited = false;
}

// for backup
void Graph::DFS(VNode &x)
{

	x.Visit();
	x.visited = true;
	int size = x.e.size();
	for (int i = 0; i < size; ++i)
	{
		VNode &w = Vertex[x.e[i]];
		if (!w.visited)
			DFS(w);
	}
}

/* 
	使用BFS获得无向图的多个最短路径，适用每条边的权重是1，
	parents 是最短路径上的回溯节点
	dfsPath 是用dfs根据parent回溯path
*/
void dfsPath(VNode &v, deque<int> path, vector<deque<int>> &result) { // path is saved by deque，result is the set of all path
	path.push_back(v.data);
	if ((v.getParents()).size() == 0) {
		result.push_back(path);
	}
	for (auto &node : v.getParents()) {	
		//cout << node.data << "++";
		dfsPath(node, path, result);
	}
	//path.pop_front();
	path.pop_back();
}

auto Graph::shortestPaths(VNode &v, VNode &u)
{
	// Iitial the graph before the next shortest path find
	for (int i = 0; i < vexNum; ++i) {
		VNode &initial_node = Vertex[i];
		initial_node.visited = false;
		initial_node.numPaths = 1;
		initial_node.layer = 0;
		initial_node.parents.clear();
	}
	queue<VNode*> q;
	// initialize the start point
	v.visited = true;
	q.push(&v);

	while (!q.empty())
	{
		VNode* _v = q.front();
		if (*_v == u) {
			break; // get the target node
		}
		q.pop();
		int size = (*_v).e.size();
		for (int i = 0; i < size; ++i) {
			VNode &w = Vertex[(*_v).e[i]];
			if (!w.visited) {           // the visited element is not added
				w.visited = true;
				w.layer = (*_v).layer + 1; //update layer
				w.numPaths = (*_v).numPaths;
				w.parents.push_back(*_v); // _v is the parent of the w
				q.push(&w);
				//w.Visit();
			}
			else {
				if (w.layer == (*_v).layer + 1) {
					w.parents.push_back(*_v);
					w.numPaths += (*_v).numPaths;
					//w.Visit();
				}
			}
		}
	} 
	vector<deque<int>> result;
 	deque<int> path;
	dfsPath(u, path, result);
	return result;
}
