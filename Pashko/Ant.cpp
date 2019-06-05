#include"Ant.h"
#include"Graph.h"
Ant::Ant() {
	Position = -1;
}
Ant::Ant(int vertexCount) {
	VertexCount = vertexCount;
}
void Ant::MoveTo(int position) {
	Position = position;
	Path[Count] = position;
	Visited[position] = true;
	Count++;
}
int Ant::getPosition() const {
	return Position;
}
bool Ant::isVisited(int vertex) const {
	return Visited[vertex];
}
const int* Ant::getPath() const {
	return Path;
}
int Ant::getCount() const {
	return Count;
}
int Ant::getStart() const {
	return Path[0];
}
int Ant::getPathNode(int index) const {
	return Path[index];
}
void Ant::Init(int position, int vertexCount) {
	Position = position;
	Count = 1;
	Path[0] = position;
	for (int i = 0; i < vertexCount; i++) {
		Visited[i] = false;
	}
	Visited[position] = true;
}