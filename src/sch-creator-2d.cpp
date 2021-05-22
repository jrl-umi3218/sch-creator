#include <sch-creator/sch-creator-2d.h>

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
			// get substring before the space and convert to double
			x = atof(fileString.substr(0, fileString.find(' ')).c_str());
			// get substring after the space and convert to double
			y = atof(fileString.substr(fileString.find(' ') + 1, fileString.length()).c_str());
			// add Vector2d to vector
			_points.push_back(Eigen::Vector2d(x,y));
			// add SCH::Point to vector
			_pointsStructure.push_back(SCH::SchCreator2D::Point(Eigen::Vector2d(x,y)));
		}
		
		file.close();
	}

	/* Listing triangles and inserting their circumcircle radius to the heap */
	void SchCreator2D::listTriangles() 
	{
		size_t n = _pointsStructure.size();

		std::vector<SchCreator2D::Radius> radii(n); // unnordered vector of Radius

		for (int i = 0; i < n; i++) {
			// make triangle and corresponding Radius
			Triangle t = Triangle(_pointsStructure[i % n].point, _pointsStructure[(i + 1) % n].point, _pointsStructure[(i + 2) % n].point);
			Radius r = Radius(i % n, (i + 1) % n, (i + 2) % n, i, t.d);
			// add triangle to list
			_triangles.push_back(t);
			// add radius to unordered vector
			radii[i] = r;
		}

		// create max heap of Radius 
		_heap = std::priority_queue<SchCreator2D::Radius>(radii.begin(), radii.end());
	}

	/* Changes the inHull Boleean in the respective point */
	void SchCreator2D::removePointFromHull(const Radius & heap) 
	{
		std::list<Triangle>::iterator it;
		size_t trianglesSize = _triangles.size();
		size_t trianglesIndex =  trianglesSize + heap.triangleIndex;
		size_t n = _pointsStructure.size();

		for(int i = -1; i <= 1; i++){
			it = _triangles.begin();
			// check if the triangle is part of the original n triangles or
			// if it's a triangle generated after eliminating a point from the sch
			if((trianglesIndex + i) >= 2 * n && heap.triangleIndex < n) {
				// if it's an original triangle, advance the iterator with respect to n
				advance(it, (n + heap.triangleIndex + i) % n);
			} else{
				// if it's a "new" triangle, advance the iterator with respect to trianglesSize
				advance(it, (trianglesIndex + i) % trianglesSize);
			}	
			
			// set the inHeap property of the triangle to false
			(*it).removeFromHeap();
		}

		// set the inHull property of the respecting Radius to false
		_pointsStructure[heap.midpointIndex].removeFromHull();
	}

	/* Verifies all points corresponding to the max heap exist in the sch */
	bool SchCreator2D::checkIfMaxHeapIsInHull()
	{
		std::list<Triangle>::iterator it = _triangles.begin();
		// move the iterator to point at the triangle corresponding to the max heap
		advance(it, (_heap.top()).triangleIndex);
		// check if all vertices of the triangle are in the Hull (inHull == true)
		bool trianglePointsExist = _pointsStructure[(_heap.top()).frontpointIndex].inHull && (_pointsStructure[(_heap.top()).midpointIndex].inHull && _pointsStructure[(_heap.top()).endpointIndex].inHull);
		// true only when all vertices are in the Hull and the triangle is in the heap
		return (*it).inHeap && trianglePointsExist;
	}

	/* Finds the closest previous point in the Hull*/
	size_t SchCreator2D::findPreviousPoint(size_t pointIndex) 
	{
		// as long as the points are in not the hull, keep looking for a point anti-clockwise
		while(!_pointsStructure[pointIndex % _pointsStructure.size()].inHull) {
			pointIndex--;
		}
		return pointIndex;
	}

	/* Finds the closest next point in the Hull */
	size_t SchCreator2D::findNextPoint(size_t pointIndex) 
	{
		// as long as the points are in not the hull, keep looking for a point clockwise
		while(!_pointsStructure[pointIndex % _pointsStructure.size()].inHull) {
			pointIndex++;
		}
		return pointIndex;
	}

	/* Makes two new triangles based on the point that must be deleted from the sch */
	void SchCreator2D::makeTriangles(size_t previousMidpoint) {
		size_t n = _pointsStructure.size();

		size_t newMidpoint = findPreviousPoint(n + previousMidpoint - 1);

		size_t newFrontpoint = findPreviousPoint(newMidpoint - 1);

		size_t newEndpoint = findNextPoint(n + previousMidpoint + 1); 

		Triangle newTriangle1 = Triangle(_pointsStructure[newFrontpoint % n].point, _pointsStructure[newMidpoint % n].point, _pointsStructure[newEndpoint % n].point);
		Radius newRadius1 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, _triangles.size(), newTriangle1.d);

		newFrontpoint = newMidpoint;
		newMidpoint = newEndpoint;
		newEndpoint = findNextPoint(newMidpoint + 1);

		Triangle newTriangle2 = Triangle(_pointsStructure[newFrontpoint % n].point, _pointsStructure[newMidpoint % n].point, _pointsStructure[newEndpoint % n].point);
		Radius newRadius2 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, _triangles.size() + 1, newTriangle2.d);

		_triangles.push_back(newTriangle1);
		_triangles.push_back(newTriangle2);
		_heap.push(newRadius1);
		_heap.push(newRadius2);
	}

	/* Stores the index of the max heap midpoint, removes it from the heap and gets the two new triangles*/
	void SchCreator2D::updateTriangles(std::list<Triangle> & triangles)
	{
		// get the index of the point to eliminate
		size_t eliminatedPointIndex = (_heap.top()).midpointIndex;
		// remove max heap
		_eliminatedVertex.push_back(_heap.top());
		_heap.pop();
		// make the mew triangles
		makeTriangles(eliminatedPointIndex);
	}

	void SchCreator2D::FindSch2D(double alpha)
	{
		_alpha = alpha;

		std::vector<Eigen::Vector2d> strictlyConvexHull;
		size_t pointsInSCH = _pointsStructure.size();

		if(pointsInSCH < 3) {
			std::cout << "You need at least 3 points.\n" << std::endl;
			exit(0);
		} 

		// make list of all initial triangles
		listTriangles();

		// while the max Heap is larger than alpha and the max Heap points are still in the Hull
		while((_heap.top()).radius > alpha && checkIfMaxHeapIsInHull()) {
			// Remove middle point and make the new triangles
			removePointFromHull(_heap.top());
			updateTriangles(_triangles);
			pointsInSCH--;

			if(pointsInSCH < 3) {
				std::cout << "\nAlpha is too small.\n" << std::endl;
				exit(0);
			}

			while(!checkIfMaxHeapIsInHull()) {
				_heap.pop();
			}
		}

		for (auto i = _pointsStructure.begin(); i != _pointsStructure.end(); i++) {
			if((*i).inHull){
				strictlyConvexHull.push_back((*i).point);
			}
		}
		
		_eliminatedPoints = _points.size() - pointsInSCH;

		makeYAML(strictlyConvexHull);
	}

	bool SchCreator2D::checkHull(const std::vector<Eigen::Vector2d> &points) 
	{
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

	// Overwriting the << operation for Eigen's Vector2d and the YAML emitter
	YAML::Emitter& operator << (YAML::Emitter& out, const Eigen::Vector2d& v) {
		out << YAML::Flow;
		out << YAML::BeginSeq << v[0] << v[1] << YAML::EndSeq;
		// Vector2d's elements added to a YAML flow
		return out;
	}

	// Overwriting the << operation for SCH's Radius and the YAML emitter
	YAML::Emitter& operator << (YAML::Emitter& out, const SCH::SchCreator2D::Radius& r) {
		out << YAML::Flow;
		out << YAML::BeginSeq << r.radius << r.midpointIndex << YAML::EndSeq;
		// Radius' radius value and eliminated point index added to a YAML flow
		return out;
	}

	void SchCreator2D::makeYAML(const std::vector<Eigen::Vector2d> &schPoints) {
		YAML::Emitter out; // make YAML Emmitter
		out << YAML::Comment("Output file containing the information about the computed strictly convex hull.");

		out << YAML::BeginMap; 
		//Begin sequence for original ch points
		out << YAML::Key << "convexHull_points";
		out << YAML::Value << YAML::BeginSeq;
		for(size_t i = 0; i <  _points.size(); i++) {
			out << _points[i];
		}
		out << YAML::EndSeq;

		//Number of eliminated points
		out << YAML::Key << "eliminated_points";
		out << YAML::Value << _eliminatedPoints;
		
		//Begin sequence for removed ch points
		out << YAML::Key << "removed_points_radius_and_index";
		out << YAML::Value << YAML::BeginSeq;
		out << YAML::Comment("[Radius, index]");
		for(size_t i = 0; i <  _eliminatedVertex.size(); i++) {
			out << _eliminatedVertex[i];
		}
		out << YAML::EndSeq;

		//Begin sequence for sch points
		out << YAML::Key << "strictlyConvexHull_points";
		out << YAML::Value << YAML::BeginSeq;
		for(size_t i = 0; i <  schPoints.size(); i++) {
			out << schPoints[i];
		}
		out << YAML::EndSeq;

		out << YAML::EndMap;

		std::ofstream fout;
		fout.open("output.yaml");
		fout << out.c_str();
		fout.close();
	}

}


int main() {
	SCH::SchCreator2D sch("C:/Users/Home/Documents/UDLAP/2021/japon/convexhull/sch/points.txt");
	sch.FindSch2D(4.5);
	//std::cout << "Is hull strictly convex? " << sch.checkHull(schPoints) << std::endl;

	return 0;
}