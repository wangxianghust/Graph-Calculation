#include "stdafx.h"
#include "Snap.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::string;
using std::ofstream;

void formatOutput(TIntFltH& index_value, const string& filePath) {
	index_value.SortByDat(false);   // rank from the max to min
	ofstream output(filePath);
	for (TIntFltH::TIter iter = index_value.BegI(); iter != index_value.EndI(); iter++) {
		output << iter->Key << "\t" << iter->Dat << endl;
	}
}

int main(int argc, char* argv[]) {
  //// what type of graph do you want to use?
  //typedef PUNGraph PGraph; // undirected graph
  //typedef PNGraph PGraph;  //   directed graph
  //typedef PNEGraph PGraph;  //   directed multigraph
  //typedef TPt<TNodeNet<TInt> > PGraph;
  //typedef TPt<TNodeEdgeNet<TInt, TInt> > PGraph;

  // this code is independent of what particular graph implementation/type we use
	printf("Creating graph:\n");

	//PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./data/test.txt"); 
	PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("./data/Email-Enron.txt");

	TIntFltH CC;
	for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
		double cc = TSnap::GetClosenessCentr(G, NI.GetId(), true);
		CC.AddDat(NI.GetId()) = cc;
		//cout << cc << endl;
	}
	formatOutput(CC, "./result/CC.txt");

	TIntFltH PageRank;
	TSnap::GetPageRank(G, PageRank);
	formatOutput(PageRank, "./result/PageRank.txt");

	TIntFltH Hub;
	TIntFltH Auth;
	TSnap::GetHits(G, Hub, Auth);
	formatOutput(Hub, "./result/HITS_Hub.txt");
	formatOutput(Auth, "./result/HITS_Auth.txt");

	TIntFltH BC;
	TSnap::GetBetweennessCentr(G, BC);
	//for (TIntFltH::TIter iter = BC.BegI(); iter != BC.EndI(); iter++) {
	//	printf("%d \t %f\n", iter->Key, iter->Dat);
	//}
	formatOutput(BC, "./result/BC.txt");

	TIntFltH Degree;
	for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
		//printf("%d %d\n", NI.GetId(), NI.GetInDeg());
		Degree.AddDat(NI.GetId()) = NI.GetDeg();
	}
	formatOutput(Degree, "./result/Degree.txt");

	printf("END:\n");
	system("PAUSE");

	//{ TFOut FOut("./result/as20graph.graph"); G->Save(FOut); } // 保存成为二进制文件
	//{ TFIn FIn("./result/as20graph.graph"); PNGraph G2 = TNGraph::Load(FIn); } // 读取这个文件
	return 0;
}

void TestEigSvd() {
  PNGraph G = TSnap::GenRndGnm<PNGraph>(100,1000, true);
  PUNGraph UG = TSnap::ConvertGraph<PUNGraph>(G);

  TSnap::SaveMatlabSparseMtx(G, "test1.mtx");
  TSnap::SaveMatlabSparseMtx(UG, "test2.mtx");

  TFltV SngValV; TVec<TFltV> LeftV, RightV;
  TSnap::GetSngVec(G, 20, SngValV, LeftV, RightV);
  printf("Singular Values:\n");
  for (int i =0; i < SngValV.Len(); i++) {
    printf("%d\t%f\n", i, SngValV[i]()); }
  printf("LEFT Singular Vectors:\n");
  for (int i=0; i < LeftV[0].Len(); i++) {
    printf("%d\t%f\t%f\t%f\t%f\t%f\n", i, LeftV[0][i](), LeftV[1][i](), LeftV[2][i](), LeftV[3][i](), LeftV[4][i]());
  }
  printf("RIGHT Singular Vectors:\n");
  for (int i=0; i < RightV[0].Len(); i++) {
    printf("%d\t%f\t%f\t%f\t%f\t%f\n", i, RightV[0][i](), RightV[1][i](), RightV[2][i](), RightV[3][i](), RightV[4][i]());
  }
  TFltV EigValV;
  TVec<TFltV> EigV;
  TSnap::GetEigVec(UG, 20, EigValV, EigV);
  printf("Eigen Values:\n");
  for (int i =0; i < EigValV.Len(); i++) {
    printf("%d\t%f\n", i, EigValV[i]()); }
  printf("Eigen Vectors %d:\n", EigV.Len());
  for (int i =0; i < EigV[0].Len(); i++) {
    printf("%d\t%f\t%f\t%f\t%f\t%f\n", i, EigV[0][i](), EigV[1][i](), EigV[2][i](), EigV[3][i](), EigV[4][i]());
  }

}
