#include<iostream>
#include<vector>
#include<ctime>
#include<chrono>
#include<omp.h>
#include"Graph.h"
using namespace std;
int main(int *arg, char** argv) {
#ifdef _OPENMP
	cout << "Parallel\n";
#endif
	int n, min, max;
	srand(time(0));
	cout << "Enter n, min distance and max distance: ";
	cin >> n>>min>>max;

	Graph g = Graph(n, min, max);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << g.distance(i, j)<<" ";
		}
		cout << endl;
	}
	int res[100];
	time_t currTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	g.antPath(1, 1, 0.5, 10, 200, 50, 50, res);
	currTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - currTime;
	int lenght = g.calculateLenght(res, n);
	cout << lenght << ":\\";
	for (int i = 0; i < n; i++) {
		cout << res[i] << " -> ";
	}
	cout << currTime << endl; 
	cout << "OpenMp:\n";
	currTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	g.antPathOpenMp(1, 1, 0.5, 10, 200, 50, 50, res);
	currTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - currTime;
	lenght = g.calculateLenght(res, n);
	cout << lenght << ":\\";
	for (int i = 0; i < n; i++) {
		cout << res[i] << " -> ";
	}
	cout << currTime << endl;
	system("PAUSE");
	return 0;
}