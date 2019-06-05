#include"Graph.h"
#include"Ant.h"
#include<iostream>
#include <random>
#include<fstream>
Graph::Graph(int n) {
	vertixesCount_ = n;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			edgesMatrix_[i][j]=-1;
		}
	}
}
Graph::Graph(int n, int minDist, int maxDist) {
	vertixesCount_ = n;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (i != j) {
				edgesMatrix_[i][j] = minDist + rand() % (maxDist - minDist + 1);
				edgesMatrix_[j][i] = edgesMatrix_[i][j];
			}
			else
				edgesMatrix_[i][j] = -1;
		}
	}
	std::ofstream fout;
	fout.open("D:\\graph.txt");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << edgesMatrix_[i][j] << (j==n-1?"\n":" ");
		}
	}
	fout.close();
}
bool  Graph::isConnected(int i, int j) const {
	return edgesMatrix_[i][j] != -1;
}
int  Graph::distance(int i, int j) const {
	return edgesMatrix_[i][j];
}
void  Graph::connect(int i, int j, int distance) {
	edgesMatrix_[i][j] = distance;
	edgesMatrix_[j][i] = distance;
}
int Graph::calculateLenght(const int* path, int pathLenght) const {
	int lenght = 0;
	for (int i = 1; i < pathLenght; i++) {
		lenght += edgesMatrix_[path[i - 1]][path[i]];
	}
	lenght += edgesMatrix_[path[0]][path[pathLenght- 1]];
	return lenght;
}
const int* Graph::antPath(double alpha, double beta, double rho, double tau0, double Q, int antCount, int iterations, int bestPath[MAX_SIZE]) const {
	double pheromones[MAX_SIZE][MAX_SIZE];
	double pheromonesChange[MAX_SIZE][MAX_SIZE];
	for (int i = 0; i < vertixesCount_; i++) {
		for (int j = 0; j < vertixesCount_; j++) {
			pheromones[i][j] = tau0;
			pheromonesChange[i][j] = 0;
		}
	}
	int bestLenght = -1;
	std::uniform_real_distribution<double> unif(0, 1);
	std::default_random_engine re;
	Ant ant = Ant(vertixesCount_);
	for (int i = 0; i < iterations; i++) {
		for (int j = 0; j < vertixesCount_; j++) {
			for (int k = 0; k < vertixesCount_; k++) {
				pheromonesChange[j][k] = 0;
			}
		}
		for (int j = 0; j < antCount; j++) {
			int lenght = 0;
			ant.Init(rand() % vertixesCount_, vertixesCount_);
			while (ant.getCount() < vertixesCount_) {
				double currentRnd = unif(re);
				double probSum = 0.0;
				double fullProb = 0.0;
				for (int k = 0; k < vertixesCount_; k++) {
					if (!ant.isVisited(k))
						probSum += pow(pheromones[ant.getPosition()][k], alpha)*pow(1.0 / distance(ant.getPosition(), k), beta);
				}
				bool isFound = false;
				for (int k = 0; k < vertixesCount_ && !isFound; k++) {
					if (!ant.isVisited(k)) {
						fullProb += pow(pheromones[ant.getPosition()][k], alpha)*
							pow(1.0 / distance(ant.getPosition(), k), beta) / probSum;
						if (currentRnd < fullProb) {
							lenght += distance(ant.getPosition(), k);
							ant.MoveTo(k);
							isFound = true;
						}
					}
				}
			}
			lenght += distance(ant.getPosition(), ant.getStart());
			for (int k = 0; k < vertixesCount_; k++) {
				pheromonesChange[ant.getPathNode(k)][ant.getPathNode((k + 1) % vertixesCount_)] += Q / lenght;
				pheromonesChange[ant.getPathNode((k + 1) % vertixesCount_)][ant.getPathNode(k)] =
					pheromonesChange[ant.getPathNode(k)][ant.getPathNode((k + 1) % vertixesCount_)];
			}
			if (bestLenght == -1 || lenght < bestLenght) {
				for (int k = 0; k < vertixesCount_; k++) {
					bestPath[k] = ant.getPathNode(k);
				}
				bestLenght = lenght;
			}
		}
		for (int j = 0; j < vertixesCount_; j++) {
			for (int k = j; k < vertixesCount_; k++) {
				pheromones[j][k] *= (1 - rho);
				if (pheromones[j][k] < tau0)
					pheromones[j][k] = tau0;
				pheromones[k][j] = pheromones[j][k];
			}
		}
		for (int j = 0; j < vertixesCount_; j++) {
			for (int k = j + 1; k < vertixesCount_; k++) {
				pheromones[j][k] += pheromonesChange[j][k];
				pheromones[k][j] = pheromones[j][k];
			}
		}
	}
	return bestPath;
}
const int* Graph::antPathOpenMp(double alpha, double beta, double rho, double tau0, double Q, int antCount, int iterations, int bestPath[MAX_SIZE]) const {
	double pheromones[MAX_SIZE][MAX_SIZE];
	double pheromonesChange[MAX_SIZE][MAX_SIZE];
	for (int i = 0; i < vertixesCount_; i++) {
		for (int j = 0; j < vertixesCount_; j++) {
			pheromones[i][j] = tau0;
			pheromonesChange[i][j] = 0;
		}
	}
	int bestLenght = -1;
	std::uniform_real_distribution<double> unif(0, 1);
	std::default_random_engine re;
	Ant ant = Ant(vertixesCount_);
	for (int i = 0; i < iterations; i++) {
		for (int j = 0; j < vertixesCount_; j++) {
			for (int k = 0; k < vertixesCount_; k++) {
				pheromonesChange[j][k] = 0;
			}
		}
		#pragma omp parallel for firstprivate(ant)
		for (int j = 0; j < antCount; j++) {
			int lenght = 0;
			ant.Init(rand() % vertixesCount_, vertixesCount_);
			while (ant.getCount() < vertixesCount_) {
				double currentRnd = unif(re);
				double probSum = 0.0;
				double fullProb = 0.0;
				for (int k = 0; k < vertixesCount_; k++) {
					if (!ant.isVisited(k))
						probSum += pow(pheromones[ant.getPosition()][k], alpha)*pow(1.0 / distance(ant.getPosition(), k), beta);
				}
				bool isFound = false;
				for (int k = 0; k < vertixesCount_ && !isFound; k++) {
					if (!ant.isVisited(k)) {
						fullProb += pow(pheromones[ant.getPosition()][k], alpha)*
							pow(1.0 / distance(ant.getPosition(), k), beta) / probSum;
						if (currentRnd < fullProb) {
							lenght += distance(ant.getPosition(), k);
							ant.MoveTo(k);
							isFound = true;
						}
					}
				}
			}
			lenght += distance(ant.getPosition(), ant.getStart());
			for (int k = 0; k < vertixesCount_; k++) {
				pheromonesChange[ant.getPathNode(k)][ant.getPathNode((k + 1) % vertixesCount_)] += Q / lenght;
				pheromonesChange[ant.getPathNode((k + 1) % vertixesCount_)][ant.getPathNode(k)] =
					pheromonesChange[ant.getPathNode(k)][ant.getPathNode((k + 1) % vertixesCount_)];
			}
			if (bestLenght == -1 || lenght < bestLenght) {
				for (int k = 0; k < vertixesCount_; k++) {
					bestPath[k] = ant.getPathNode(k);
				}
				bestLenght = lenght;
			}
		}
		#pragma omp parallel for
		for (int j = 0; j < vertixesCount_; j++) {
			for (int k = j; k < vertixesCount_; k++) {
				pheromones[j][k] *= (1 - rho);
				if (pheromones[j][k] < tau0)
					pheromones[j][k] = tau0;
				pheromones[k][j] = pheromones[j][k];
			}
		}
		#pragma omp parallel for
		for (int j = 0; j < vertixesCount_; j++) {
			for (int k = j + 1; k < vertixesCount_; k++) {
				pheromones[j][k] += pheromonesChange[j][k];
				pheromones[k][j] = pheromones[j][k];
			}
		}
	}
	return bestPath;
}
//std::vector<int> Graph::antPathOpenMp(double alpha, double beta, double rho, double tau0, double Q, int antCount, int iterations) const {
//	std::vector<std::vector<double>> pheromones = std::vector<std::vector<double>>(vertixesCount_, std::vector<double>(vertixesCount_, tau0));
//	std::vector<std::vector<double>> pheromonesChange = std::vector<std::vector<double>>(vertixesCount_, std::vector<double>(vertixesCount_, 0));
//	std::vector<int> bestPath= std::vector<int>(vertixesCount_);;
//	int bestLenght = -1;
//	std::uniform_real_distribution<double> unif(0, 1);
//	std::default_random_engine re;
//	Ant ant = Ant(vertixesCount_);
//	for (int i = 0; i < iterations; i++) {
//		int lenght = 0;
//		for (int j = 0; j < vertixesCount_; j++) {
//			for (int k = 0; k < vertixesCount_; k++) {
//				pheromonesChange[j][k] = 0;
//			}
//		}
//		//#pragma omp parallel for firstprivate(ant)
//		for (int j = 0; j < antCount; j++) {
//			ant.Init(rand()%vertixesCount_, vertixesCount_);
//			while (ant.getCount() < vertixesCount_) {
//				double currentRnd = unif(re);
//				double probSum = 0.0;
//				double fullProb = 0.0;
//				for (int k = 0; k < vertixesCount_; k++) {
//					if (!ant.isVisited(k))
//						probSum += pow(pheromones[ant.getPosition()][k], alpha)*pow(1.0 / distance(ant.getPosition(), k), beta);
//				}
//				bool isFound = false;
//				for (int k = 0; k < vertixesCount_ && !isFound; k++) {
//					if (!ant.isVisited(k)) {
//						fullProb += pow(pheromones[ant.getPosition()][k], alpha)*
//							pow(1.0 / distance(ant.getPosition(), k), beta) / probSum;
//						if (currentRnd < fullProb) {
//							lenght += distance(ant.getPosition(), k);
//							ant.MoveTo(k);
//							isFound = true;
//						}
//					}
//				}
//			}
//			lenght += distance(ant.getPosition(), ant.getStart());
//			for (int k = 0; k < vertixesCount_; k++) {
//				pheromonesChange[ant.getPathNode(k)][ant.getPathNode((k + 1) % vertixesCount_)] += Q / lenght;
//				pheromonesChange[ant.getPathNode((k + 1) % vertixesCount_)][ant.getPathNode(k)] =
//					pheromonesChange[ant.getPathNode(k)][ant.getPathNode((k + 1) % vertixesCount_)];
//			}
//			if (bestLenght == -1 || lenght < bestLenght) {
//				for (int k = 0; k < vertixesCount_; k++) {
//					bestPath[k] = ant.getPathNode(k);
//				}
//				bestLenght = lenght;
//			}
//		}
//		//#pragma omp parallel for
//		for (int j = 0; j < vertixesCount_; j++) {
//			for (int k = j; k < vertixesCount_; k++) {
//				pheromones[j][k] *= (1 - rho);
//				if (pheromones[j][k] < tau0)
//					pheromones[j][k] = tau0;
//				pheromones[k][j] = pheromones[j][k];
//			}
//		}
//		//#pragma omp parallel for
//		for (int j = 0; j < vertixesCount_; j++) {
//			for (int k = j + 1; k < vertixesCount_; k++) {
//				pheromones[j][k] += pheromonesChange[j][k];
//				pheromones[k][j] = pheromones[j][k];
//			}
//		}
//	}
//	return bestPath;
//}