#include <sch-creator/sch-creator-2d.h>
// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>


// PYBIND11_MODULE(Sch2D, handle) {
// 	py::class_<SchCreator2D>(
// 		handle, "PySchCreator2D"
// 	)

// 	.def(py::init<std::string CONS&>());
// 	.def("FindSch2D", [](SchCreator2D &self, double alpha) {
// 		return py::return_value_policy::reference_internal
// 	});
// }

namespace SCH
{
	/* Strictly Convex Hull */
	SchCreator2D::SchCreator2D(const std::string &points)
	{
		_pointsPath = points;
		readPointsFromFile();

	}

	void SchCreator2D::readPointsFromFile()
	{
		std::string fileString;
		std::ifstream file;
		double x, y;

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
			_points.push_back(Eigen::Vector2d(x,y));
			// create vector of SCH::Point
			_pointsStructure.push_back(SCH::SchCreator2D::Point(Eigen::Vector2d(x,y)));
		}
		
		file.close();
	}

	/* Heap insertion while maintaining heap */
	void SchCreator2D::insertElementToHeap(Radius & a) {
		_heap.push(a);
	}

	/* Listing triangles and inserting their circumcircle radius to the heap */
	std::list<SchCreator2D::Triangle> SchCreator2D::listTriangles(std::vector<Point> & points) 
	{
		std::list<SchCreator2D::Triangle> triangles;
		size_t n = points.size();

		std::vector<SchCreator2D::Radius> radii(n);

		for (int i = 0; i < n; i++) {
			Triangle t = Triangle(points[i % n].point, points[(i + 1) % n].point, points[(i + 2) % n].point);
			Radius r = Radius(i % n, (i + 1) % n, (i + 2) % n, i, t.d);
			triangles.push_back(t);
			radii[i] = r;
		}

		_heap = std::priority_queue<SchCreator2D::Radius>(radii.begin(), radii.end());

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

	/* Remove max heap while mantaining heap */
	void SchCreator2D::removeHeap(std::priority_queue<Radius> & heap){
		heap.pop();
	}

	/* Verifies all points corresponding to the max heap exist in the sch */
	bool SchCreator2D::checkIfMaxHeapIsInHull(std::list<Triangle> & triangles){
		std::list<Triangle>::iterator it = triangles.begin();
		advance(it, (_heap.top()).triangleIndex);
		bool trianglePointsExist = _pointsStructure[(_heap.top()).frontpointIndex].inHull && (_pointsStructure[(_heap.top()).midpointIndex].inHull && _pointsStructure[(_heap.top()).endpointIndex].inHull);
		return (*it).inHeap && trianglePointsExist;
	}

	/* Finds the closest previous point in the Hull*/
	size_t SchCreator2D::findPreviousPoint(size_t pointIndex) 
	{
		while(!_pointsStructure[pointIndex % _pointsStructure.size()].inHull) {
			pointIndex--;
		}
		return pointIndex;
	}

	/* Finds the closest next point in the Hull */
	size_t SchCreator2D::findNextPoint(size_t pointIndex) 
	{
		while(!_pointsStructure[pointIndex % _pointsStructure.size()].inHull) {
			pointIndex++;
		}
		return pointIndex;
	}

	/* Makes two new triangles based on the point that must be deleted from the sch */
	void SchCreator2D::makeTriangles(std::list<Triangle> & triangles, size_t previousMidpoint) {
		size_t n = _pointsStructure.size();

		size_t newMidpoint = findPreviousPoint(n + previousMidpoint - 1);

		size_t newFrontpoint = findPreviousPoint(newMidpoint - 1);

		size_t newEndpoint = findNextPoint(n + previousMidpoint + 1); 

		Triangle newTriangle1 = Triangle(_pointsStructure[newFrontpoint % n].point, _pointsStructure[newMidpoint % n].point, _pointsStructure[newEndpoint % n].point);
		Radius newRadius1 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, triangles.size(), newTriangle1.d);

		newFrontpoint = newMidpoint;
		newMidpoint = newEndpoint;
		newEndpoint = findNextPoint(newMidpoint + 1);

		Triangle newTriangle2 = Triangle(_pointsStructure[newFrontpoint % n].point, _pointsStructure[newMidpoint % n].point, _pointsStructure[newEndpoint % n].point);
		Radius newRadius2 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, triangles.size() + 1, newTriangle2.d);

		triangles.push_back(newTriangle1);
		triangles.push_back(newTriangle2);
		insertElementToHeap(newRadius1);
		insertElementToHeap(newRadius2);
	}

	/* Stores the index of the max heap midpoint, removes it from the heap and gets the two new triangles*/
	void SchCreator2D::updateTriangles(std::list<Triangle> & triangles)
	{
		size_t eliminatedPointIndex = (_heap.top()).midpointIndex;
		removeHeap(_heap);
		makeTriangles(triangles, eliminatedPointIndex);
	}

	

	std::vector<Eigen::Vector2d> SchCreator2D::FindSch2D(double alpha)
	{
		_alpha = alpha;

		std::vector<Eigen::Vector2d> strictlyConvexHull;
		size_t pointsInSCH = _pointsStructure.size();

		if(_pointsStructure.size() < 3) {
			std::cout << "You need at least 3 points.\n" << std::endl;
			return {};
		} 

		

		std::list<Triangle> triangles = listTriangles(_pointsStructure);
		
		while((_heap.top()).radius > alpha && checkIfMaxHeapIsInHull(triangles)) {
			removePointFromHull(_pointsStructure, triangles, _heap.top());
			updateTriangles(triangles);
			pointsInSCH--;

			if(pointsInSCH < 3) {
				std::cout << "\nAlpha is too small.\n" << std::endl;
				return {};
			}

			while(!checkIfMaxHeapIsInHull(triangles)) {
				removeHeap(_heap);
			}
		}

		for (auto i = _pointsStructure.begin(); i != _pointsStructure.end(); i++) {
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
	SCH::SchCreator2D::Radius rad(2,3,4,6,2.5);
	SCH::SchCreator2D::Radius radb(4,5,6,2,6.45);
	std::cout << (rad < radb) << std::endl;
	exit(0);

	SCH::SchCreator2D sch("C:/Users/Home/Documents/UDLAP/2021/japon/convexhull/sch/points.txt");
	std::vector<Eigen::Vector2d> schPoints = sch.FindSch2D(4.5);
	std::cout << "Is hull strictly convex? " << sch.checkHull(schPoints) << std::endl;

	std::cout << "\nPoints in strictly convex hull" << std::endl;
	for(auto i = schPoints.begin(); i != schPoints.end(); i++) std::cout << *i << '\n' << std::endl;


	return 0;
}