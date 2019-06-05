#pragma once
#include<vector>
class Graph {
private:
	static const int MAX_SIZE = 100;
	int vertixesCount_;
	int edgesMatrix_[MAX_SIZE][MAX_SIZE];
public:
	Graph(int n);
	Graph(int n, int minDist, int maxDist);
	Graph(int vertixes, std::vector<std::vector<int>> edges);
	bool  isConnected(int i, int j) const;
	int  distance(int i, int j) const;
	void  connect(int i, int j, int distance);
	int calculateLenght(const int* path, int pathLenght) const;
	const int* antPath(double alpha, double beta, double rho, double tau0, double Q, int antCount, int iterations, int bestPath[MAX_SIZE]) const;
	const int* antPathOpenMp(double alpha, double beta, double rho, double tau0, double Q, int antCount, int iterations, int bestPath[MAX_SIZE]) const;
};