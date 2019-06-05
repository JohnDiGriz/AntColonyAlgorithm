#pragma once
#include<vector>
class Ant {
private:
	static const int MAX_SIZE = 100;
	int Position;
	int Path[MAX_SIZE];
	bool Visited[MAX_SIZE];
	int Count;
	int VertexCount;
public:
	Ant();
	Ant(int vertexCount);
	void MoveTo(int position);
	int getPosition() const;
	bool isVisited(int vertex) const;
	const int* getPath() const;
	int getCount() const;
	int getStart() const;
	int getPathNode(int index) const;
	void Init(int position, int vertexCount);
};