#include <sch-creator/sch-creator-2d.h>
#include <Eigen/Core>
#include <iostream>
#include <list>
#include <math.h>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

namespace SCH
{
	/* Strictly Convex Hull */
	SchCreator2D::SchCreator2D(vector<Vector2d> points, double alpha)
	{
		_alpha = alpha;
		_points = points;
	}

	void SchCreator2D::swap(Radius &a, Radius &b) 
	{
		Radius temp = a;
		a = b;
		b = temp;
	}

	/* Heap insertion while maintaining heap */
	void SchCreator2D::insertElementToHeap(vector<Radius> &heap, Radius a) {
		heap.push_back(a);
		if(heap.size() > 1){
			int childIndex = heap.size();
			int parentIndex = int(childIndex / 2);
			
			while(heap[parentIndex - 1].radius < a.radius){
				swap(heap[parentIndex - 1], heap[childIndex - 1]);
				if(parentIndex - 2 < 0) break;
				childIndex = parentIndex;
				parentIndex = int(childIndex / 2);
			}
		}
	}

	/* Listing triangles and inserting their circumcircle radius to the heap */
	list<Triangle> SchCreator2D::listTriangles(vector<Point> points, vector<Radius> &heap) {
		list<Triangle> triangles;
		int n = points.size();
		
		for (int i = 0; i < n; i++) {
			Triangle t = Triangle(points[i % n].point, points[(i + 1) % n].point, points[(i + 2) % n].point);
			Radius r = Radius(i % n, (i + 1) % n, (i + 2) % n, i, t.d);
			triangles.push_back(t);
			insertElementToHeap(heap, r);
		}

		return triangles;
	}

	/* Changes the inHull Boleean in the respective point */
	void SchCreator2D::removePointFromHull(vector<Point> &points, list<Triangle> &triangles, const Radius &heap) {
		list<Triangle>::iterator it;
		int trianglesSize = triangles.size();
		int trianglesIndex =  trianglesSize + heap.triangleIndex;
		for(int i = -1; i <= 1; i++){
			it = triangles.begin();
			if((trianglesIndex + i) >= 2 * points.size() && heap.triangleIndex < points.size()) {
				advance(it, (points.size() + heap.triangleIndex + i) % points.size());
			} else{
				advance(it, (trianglesIndex + i) % trianglesSize);
			}	

			(*it).removeFromHeap();
		}
		points[heap.midpointIndex].removeFromHull();
	}

	/* Heap operation, comparing node to left and right child */
	void SchCreator2D::compareToChild(vector<Radius> &heap) {
		int parentIndex = 1;
		int childIndex = 2 * parentIndex;
		while(childIndex < heap.size() && heap[parentIndex - 1].radius < heap[childIndex - 1].radius || heap[parentIndex - 1].radius < heap[childIndex].radius){
			if(heap[childIndex - 1].radius < heap[childIndex].radius) childIndex++;
			while(heap[parentIndex - 1].radius < heap[childIndex - 1].radius){
				swap(heap[parentIndex - 1], heap[childIndex - 1]);
				parentIndex = childIndex;
				childIndex = 2 * parentIndex;
			}
		}
	}

	/* Remove max heap while mantaining heap */
	void SchCreator2D::removeHeap(vector<Radius> &heap){
		swap(heap.front(), heap.back());
		heap.pop_back();
		int parentIndex = 1;
		int childIndex = 2 * parentIndex;
		while(childIndex + 1 < heap.size() && (heap[parentIndex - 1].radius < heap[childIndex - 1].radius || heap[parentIndex - 1].radius < heap[childIndex].radius)){
			if(heap[childIndex - 1].radius < heap[childIndex].radius) childIndex++;
			if(heap[parentIndex - 1].radius < heap[childIndex - 1].radius){
				swap(heap[parentIndex - 1], heap[childIndex - 1]);
				parentIndex = childIndex;
				childIndex = 2 * parentIndex;
			}
		}
	}

	/* Verifies all points corresponding to the max heap exist in the sch */
	bool SchCreator2D::checkIfMaxHeapIsInHull(vector<Radius> &heap, vector<Point> points, list<Triangle> &triangles){
		list<Triangle>::iterator it = triangles.begin();
		advance(it, heap[0].triangleIndex);
		bool trianglePointsExist = points[heap[0].frontpointIndex].inHull && (points[heap[0].midpointIndex].inHull && points[heap[0].endpointIndex].inHull);
		return (*it).inHeap && trianglePointsExist;
	}

	/* Finds the closest previous point in the Hull*/
	int SchCreator2D::findPreviousPoint(const vector<Point> &points, int pointIndex) {
		while(!points[pointIndex % points.size()].inHull) {
			pointIndex--;
		}
		return pointIndex;
	}

	/* Finds the closest next point in the Hull */
	int SchCreator2D::findNextPoint(const vector<Point> &points, int pointIndex) {
		while(!points[pointIndex % points.size()].inHull) {
			pointIndex++;
		}
		return pointIndex;
	}

	/* Makes two new triangles based on the point that must be deleted from the sch */
	void SchCreator2D::makeTriangles(const vector<Point> &points, list<Triangle> &triangles, vector<Radius> &heap, int previousMidpoint) {
		int n = points.size();

		int newMidpoint = findPreviousPoint(points, n + previousMidpoint - 1);

		int newFrontpoint = findPreviousPoint(points, newMidpoint - 1);

		int newEndpoint = findNextPoint(points, n + previousMidpoint + 1); 

		Triangle newTriangle1 = Triangle(points[newFrontpoint % n].point, points[newMidpoint % n].point, points[newEndpoint % n].point);
		Radius newRadius1 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, triangles.size(), newTriangle1.d);

		newFrontpoint = newMidpoint;
		newMidpoint = newEndpoint;
		newEndpoint = findNextPoint(points, newMidpoint + 1);

		Triangle newTriangle2 = Triangle(points[newFrontpoint % n].point, points[newMidpoint % n].point, points[newEndpoint % n].point);
		Radius newRadius2 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, triangles.size() + 1, newTriangle2.d);

		triangles.push_back(newTriangle1);
		triangles.push_back(newTriangle2);
		insertElementToHeap(heap, newRadius1);
		insertElementToHeap(heap, newRadius2);
	}

	/* Stores the index of the max heap midpoint, removes it from the heap and gets the two new triangles*/
	void SchCreator2D::updateTriangles(const vector<Point> &points, list<Triangle> &triangles, vector<Radius> &heap ){
		int eliminatedPointIndex = heap.front().midpointIndex;
		removeHeap(heap);
		makeTriangles(points, triangles, heap, eliminatedPointIndex);
	}

	

	vector<Vector2d> SchCreator2D::FindSch2D(){
		vector<Point> points;
		for(auto i = _points.begin(); i != _points.end(); i++) points.push_back(Point((*i)));

		if(points.size() < 3) {
			cout << "You need at least 3 points.\n" << endl;
			return {};
		} 
		vector<Vector2d> strictlyConvexHull;
		int pointsInSCH = points.size();

		vector<Radius> radiusHeap;
		list<Triangle> triangles = listTriangles(points, radiusHeap);
		
		while(radiusHeap[0].radius > _alpha && checkIfMaxHeapIsInHull(radiusHeap, points, triangles)) {
			removePointFromHull(points, triangles, radiusHeap[0]);
			updateTriangles(points, triangles, radiusHeap);
			pointsInSCH--;

			if(pointsInSCH < 3) {
				cout << "\nAlpha is too small.\n" << endl;
				return {};
			}

			while(!checkIfMaxHeapIsInHull(radiusHeap, points, triangles)) {
				removeHeap(radiusHeap);
			}
		}

		for (auto i = points.begin(); i != points.end(); i++) {
			if((*i).inHull){
				strictlyConvexHull.push_back((*i).point);
			}
		}

		return strictlyConvexHull;
	}

	bool SchCreator2D::checkHull(const vector<Vector2d> &points) {
		int n = points.size();
		Vector2d circleCenter;
		
		for(int i = 0; i < n; i++) {
			Vector2d p1 = points[i % n];
			Vector2d p2 = points[(i + 1) % n];
			Vector2d pa = (p2-p1)/2;
			Vector2d rhombusCenter = p1 + pa;
			double a = pa.norm();
			double b = sqrt(pow(_alpha,2) - pow(a,2));

			Vector2d distanceToRhombusCenter = Vector2d(b * (p2[1] - rhombusCenter[1]) / a, -b * (p2[0] - rhombusCenter[0]) / a); 
			circleCenter = rhombusCenter + distanceToRhombusCenter;

			double m = _points.size();
			for(int j = 0; j < m; j++) {
				double distancePointToCenter = (_points[j] - circleCenter).norm();

				if(distancePointToCenter > _alpha + 1e-5) return false;
			}
		}

		return true;
	}

	
}

int main() {
	string pointsInfo;
	ifstream readFile;
	double x, y;
	vector<Vector2d> points;

	readFile.open("points.txt");
	if(!readFile) {
		cout << "File not found." << endl;
		exit(1);
	}

	while(getline(readFile, pointsInfo)) {
		// get substring and convert to double
		x = atof(pointsInfo.substr(0, pointsInfo.find(' ')).c_str());
		y = atof(pointsInfo.substr(pointsInfo.find(' ') + 1, pointsInfo.length()).c_str());
		// add point to vector
		points.push_back(Vector2d(x, y));
	}
	SCH::SchCreator2D sch(points,4.5);
	vector<Vector2d> schPoints = sch.FindSch2D();
	cout << "Is hull strictly convex? " << sch.checkHull(schPoints) <<endl;

	cout << "\nPoints in strictly convex hull" << endl;
	for(auto i = schPoints.begin(); i != schPoints.end(); i++) cout << *i << '\n' << endl;

	return 0;
}