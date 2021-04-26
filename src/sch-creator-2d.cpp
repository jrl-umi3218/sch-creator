#include <sch-creator/sch-creator-2d.h>

namespace SCH
{
	/* Strictly Convex Hull */
	SchCreator2D::SchCreator2D(const std::string &points)
	{
		_pointsPath = points;
		readPointsFromFile();

	}

	void SchCreator2D::swap(Radius &a, Radius &b) 
	{
		Radius temp = a;
		a = b;
		b = temp;
	}

	void SchCreator2D::readPointsFromFile()
	{
		std::string fileString;
		std::ifstream file;
		double x, y;
		std::vector<Eigen::Vector2d> newPoints;

		file.open(_pointsPath);
		if(!file) {
			std::cout << "File not found." << std::endl;
			exit(1);
		}

		while(std::getline(file, fileString)) {
			// get substring and convert to double
			x = atof(fileString.substr(0, fileString.find(' ')).c_str());
			y = atof(fileString.substr(fileString.find(' ') + 1, fileString.length()).c_str());
			// add point to vector
			newPoints.push_back(Eigen::Vector2d(x,y));
		}
		
		file.close();

		_points = newPoints;
	}

	/* Heap insertion while maintaining heap */
	void SchCreator2D::insertElementToHeap(std::vector<Radius> & heap, Radius & a) {
		heap.push_back(a);
		if(heap.size() > 1){
			size_t childIndex = heap.size();
			size_t parentIndex = (childIndex / 2);
			
			while(heap[parentIndex - 1].radius < a.radius){
				swap(heap[parentIndex - 1], heap[childIndex - 1]);
				if(parentIndex - 2 < 0) break;
				childIndex = parentIndex;
				parentIndex = int(childIndex / 2);
			}
		}
	}

	/* Listing triangles and inserting their circumcircle radius to the heap */
	std::list<SchCreator2D::Triangle> SchCreator2D::listTriangles(std::vector<Point> & points, std::vector<Radius> & heap) 
	{
		std::list<SchCreator2D::Triangle> triangles;
		size_t n = points.size();
		
		for (int i = 0; i < n; i++) {
			Triangle t = Triangle(points[i % n].point, points[(i + 1) % n].point, points[(i + 2) % n].point);
			Radius r = Radius(i % n, (i + 1) % n, (i + 2) % n, i, t.d);
			triangles.push_back(t);
			insertElementToHeap(heap, r);
		}

		return triangles;
	}

	/* Changes the inHull Boleean in the respective point */
	void SchCreator2D::removePointFromHull(std::vector<Point> & points, std::list<Triangle> & triangles, const Radius & heap) {
		std::list<Triangle>::iterator it;
		size_t trianglesSize = triangles.size();
		size_t trianglesIndex =  trianglesSize + heap.triangleIndex;
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
	void SchCreator2D::compareToChild(std::vector<Radius> & heap) {
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
	void SchCreator2D::removeHeap(std::vector<Radius> & heap){
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
	bool SchCreator2D::checkIfMaxHeapIsInHull(std::vector<Radius> & heap, std::vector<Point> & points, std::list<Triangle> & triangles){
		std::list<Triangle>::iterator it = triangles.begin();
		advance(it, heap[0].triangleIndex);
		bool trianglePointsExist = points[heap[0].frontpointIndex].inHull && (points[heap[0].midpointIndex].inHull && points[heap[0].endpointIndex].inHull);
		return (*it).inHeap && trianglePointsExist;
	}

	/* Finds the closest previous point in the Hull*/
	size_t SchCreator2D::findPreviousPoint(const std::vector<Point> & points, size_t pointIndex) {
		while(!points[pointIndex % points.size()].inHull) {
			pointIndex--;
		}
		return pointIndex;
	}

	/* Finds the closest next point in the Hull */
	size_t SchCreator2D::findNextPoint(const std::vector<Point> & points, size_t pointIndex) {
		while(!points[pointIndex % points.size()].inHull) {
			pointIndex++;
		}
		return pointIndex;
	}

	/* Makes two new triangles based on the point that must be deleted from the sch */
	void SchCreator2D::makeTriangles(const std::vector<Point> & points, std::list<Triangle> & triangles, std::vector<Radius> & heap, size_t previousMidpoint) {
		size_t n = points.size();

		size_t newMidpoint = findPreviousPoint(points, n + previousMidpoint - 1);

		size_t newFrontpoint = findPreviousPoint(points, newMidpoint - 1);

		size_t newEndpoint = findNextPoint(points, n + previousMidpoint + 1); 

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
	void SchCreator2D::updateTriangles(const std::vector<Point> & points, std::list<Triangle> & triangles, std::vector<Radius> & heap ){
		int eliminatedPointIndex = heap.front().midpointIndex;
		removeHeap(heap);
		makeTriangles(points, triangles, heap, eliminatedPointIndex);
	}

	

	std::vector<Eigen::Vector2d> SchCreator2D::FindSch2D(double alpha)
	{
		_alpha = alpha;
		std::vector<Point> points;
		for(auto i = _points.begin(); i != _points.end(); i++) points.push_back(Point((*i)));

		if(points.size() < 3) {
			std::cout << "You need at least 3 points.\n" << std::endl;
			return {};
		} 

		std::vector<Eigen::Vector2d> strictlyConvexHull;
		size_t pointsInSCH = points.size();

		std::vector<Radius> radiusHeap;
		std::list<Triangle> triangles = listTriangles(points, radiusHeap);
		
		while(radiusHeap[0].radius > alpha && checkIfMaxHeapIsInHull(radiusHeap, points, triangles)) {
			removePointFromHull(points, triangles, radiusHeap[0]);
			updateTriangles(points, triangles, radiusHeap);
			pointsInSCH--;

			if(pointsInSCH < 3) {
				std::cout << "\nAlpha is too small.\n" << std::endl;
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

	bool SchCreator2D::checkHull(const std::vector<Eigen::Vector2d> &points) {
		size_t n = points.size();
		Eigen::Vector2d circleCenter;
		
		for(size_t i = 0; i < n; i++) {
			Eigen::Vector2d p1 = points[i % n];
			Eigen::Vector2d p2 = points[(i + 1) % n];
			Eigen::Vector2d pa = (p2-p1)/2;
			Eigen::Vector2d rhombusCenter = p1 + pa;
			double a = pa.norm();
			double b = sqrt(pow(_alpha,2) - pow(a,2));

			Eigen::Vector2d distanceToRhombusCenter = Eigen::Vector2d(b * (p2[1] - rhombusCenter[1]) / a, -b * (p2[0] - rhombusCenter[0]) / a); 
			circleCenter = rhombusCenter + distanceToRhombusCenter;

			size_t m = _points.size();
			for(size_t j = 0; j < m; j++) {
				double distancePointToCenter = (_points[j] - circleCenter).norm();

				if(distancePointToCenter > _alpha + 1e-5) return false;
			}
		}

		return true;
	}
}

int main() {
	SCH::SchCreator2D sch("C:/Users/Home/Documents/UDLAP/2021/japon/convexhull/sch/points.txt");
	std::vector<Eigen::Vector2d> schPoints = sch.FindSch2D(4.5);
	std::cout << "Is hull strictly convex? " << sch.checkHull(schPoints) << std::endl;

	std::cout << "\nPoints in strictly convex hull" << std::endl;
	for(auto i = schPoints.begin(); i != schPoints.end(); i++) std::cout << *i << '\n' << std::endl;

	return 0;
}