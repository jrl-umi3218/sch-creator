#include <sch-creator/sch-creator-2d.h>

namespace SCH
{
	/* Strictly Convex Hull */
	SchCreator2D::SchCreator2D(const std::string &points)
	{
		_pointsPath = points;
		readPointsFromFile();
	}

	/*	reads the files from the given path and builds the
		arrays to be procesed throughout the algorithm */
	void SchCreator2D::readPointsFromFile()
	{
		std::string fileString;
		std::ifstream file;
		double x, y;

		// open the file
		file.open(_pointsPath);
		//verify file is correctly opened, exit otherwise
		if(!file) {
			std::cout << "File not found." << std::endl;
			exit(1);
		}

		// build arrays of points
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
		
		// close the file
		file.close();
	}

	/* Listing triangles and inserting their circumcircle radius to the heap */
	void SchCreator2D::listTriangles() 
	{
		// get the amount of points
		size_t n = _pointsStructure.size();

		// declare the unnordered vector of Radius
		std::vector<SchCreator2D::Radius> radii(n); 

		for (int i = 0; i < n; i++) {
			// make triangle and corresponding Radius
			Triangle t = Triangle(_pointsStructure[i % n].point, 
						_pointsStructure[(i + 1) % n].point, 
						_pointsStructure[(i + 2) % n].point);
			Radius r = Radius(i % n, (i + 1) % n, (i + 2) % n, i, t.d);
			// add triangle to list
			_triangles.push_back(t);
			// add radius to unordered vector
			radii[i] = r;
		}

		// create max heap of Radius 
		_heap = std::priority_queue<SchCreator2D::Radius>(radii.begin(), radii.end());
	}

	/* 	changes the in Hull state of the removed point and changes its corresponding
		triangle's inHeap state to false */
	void SchCreator2D::removePointFromHull(const Radius & maxHeap) 
	{
		// create iterator for the list of triangles
		std::list<Triangle>::iterator it;
		// get triangles size, index and amount of points
		size_t trianglesSize = _triangles.size();
		size_t trianglesIndex =  trianglesSize + maxHeap.triangleIndex;
		size_t n = _pointsStructure.size();

		// initialize iterator to the begining of the list
		it = _triangles.begin();
		// advance iterator to point at triangle
		advance(it, (trianglesIndex + trianglesSize) % trianglesSize);
		// set the inHeap property of the triangle to false
		(*it).removeFromHeap();

		// set the inHull property of the corresponding midpoint to false
		_pointsStructure[maxHeap.midpointIndex].removeFromHull();
	}

	/* 	Verifies all points corresponding to the max heap exist in the sch */
	bool SchCreator2D::checkIfMaxHeapIsInHull()
	{
		// initialize iterator at the begining of the list
		std::list<Triangle>::iterator it = _triangles.begin();

		// move the iterator to point at the triangle corresponding to the max heap
		advance(it, (_heap.top()).triangleIndex);

		// check if all vertices of the triangle are in the Hull (inHull == true)
		bool trianglePointsExist = _pointsStructure[(_heap.top()).startpointIndex].inHull 
									&& (_pointsStructure[(_heap.top()).midpointIndex].inHull 
									&& _pointsStructure[(_heap.top()).endpointIndex].inHull);

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

		// get the three points of the first new triangle
		size_t newMidpoint = findPreviousPoint(n + previousMidpoint - 1) % n;
		size_t newStartpoint = findPreviousPoint(n + newMidpoint - 1) % n;
		size_t newEndpoint = findNextPoint(previousMidpoint + 1) % n; 

		// build the first new triangle and radius
		Triangle newTriangle1 = Triangle(_pointsStructure[newStartpoint].point,
										_pointsStructure[newMidpoint].point, 
										_pointsStructure[newEndpoint].point);
		Radius newRadius1 = Radius(newStartpoint, newMidpoint, newEndpoint, 
										_triangles.size(), newTriangle1.d);

		// get the three points of the second new triangle
		newStartpoint = newMidpoint % n;
		newMidpoint = newEndpoint % n;
		newEndpoint = findNextPoint(newMidpoint + 1) % n;

		// build the second triangle and radius
		Triangle newTriangle2 = Triangle(_pointsStructure[newStartpoint].point, 
										_pointsStructure[newMidpoint].point, 
										_pointsStructure[newEndpoint].point);
		Radius newRadius2 = Radius(newStartpoint, newMidpoint, newEndpoint, 
										_triangles.size() + 1, newTriangle2.d);
		
		// add the triangles to the list 
		_triangles.push_back(newTriangle1);
		_triangles.push_back(newTriangle2);
		// add radius to the heap
		_heap.push(newRadius1);
		_heap.push(newRadius2);
	}

	/* 	Stores the index of the max heap midpoint, removes it from the heap 
		and gets the two new triangles*/
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

	/*	processes the points and obtains the ones belonging to the sch,
		generates the required data for the output file and generates it */
	void SchCreator2D::FindSch2D(double alpha)
	{
		_alpha = alpha;
		
		// get the initial number of points for the sch
		size_t pointsInSCH = _pointsStructure.size();

		// ensure there's at least 3 points, exit program otherwise
		if(pointsInSCH < 3) {
			std::cout << "You need at least 3 points.\n" << std::endl;
			exit(1);
		} 

		// make list of all initial triangles
		listTriangles();

		// while the max Heap radius is larger than alpha
		while((_heap.top()).radius > alpha) {
			// Remove middle point and make the new triangles
			removePointFromHull(_heap.top());
			updateTriangles(_triangles);
			// decrease the amount of points in the sch
			pointsInSCH--;
			// ensure that, after elimination, there's still at least 3 points
			if(pointsInSCH < 3) {
				std::cout << "\nAlpha is too small.\n" << std::endl;
				exit(1);
			}

			// verify that all points in the new max heap belong to the hull,
			// eliminate the max heap otherwise
			while(!checkIfMaxHeapIsInHull()) {
				_heap.pop();
			}
		}

		// add all points that are still part of the hull to the schPoints array
		for (auto i = _pointsStructure.begin(); i != _pointsStructure.end(); i++) {
			if((*i).inHull){
				_schPoints.push_back((*i).point);
			}
		}

		// get the amount of eliminated points
		_eliminatedPoints = _points.size() - pointsInSCH;

		// generate the output YAML file
		makeYAML();
	}

	/*	check to ensure the sch is actually strictly convex */
	bool SchCreator2D::checkHull() 
	{
		std::cout << "\n\n********** CHECK HULL **********\n" << std::endl;
		size_t n = _schPoints.size();
		Eigen::Vector2d circleCenter;
		double distancePointToCenter;
		
		// for all points in the sch
		for(size_t i = 0; i < n; i++) {
			std::cout << "------------------------------" << std::endl;
			std::cout << "Measuring from point " << i << std::endl;
			// get two contiguous points
			Eigen::Vector2d p1 = _schPoints[i % n];
			Eigen::Vector2d p2 = _schPoints[(i + 1) % n];
			// get distance between the points and rhombus center
			Eigen::Vector2d pa = (p2-p1)/2;
			Eigen::Vector2d rhombusCenter = p1 + pa;
			// get a and b
			double a = pa.norm();
			double b = sqrt(pow(_alpha,2) - pow(a,2));

			// get the distance to the center of the rhombus
			Eigen::Vector2d distanceToRhombusCenter = Eigen::Vector2d(b * (p2[1] - rhombusCenter[1]) / a,
																	-b * (p2[0] - rhombusCenter[0]) / a); 
			// get the center of the circle going through pq and p2
			circleCenter = rhombusCenter + distanceToRhombusCenter;

			// get the amount of points in the original ch
			size_t m = _points.size();
			// for all points in the ch
			for(size_t j = 0; j < m; j++) {
				// calculate distance from the center of the circle to the point
				distancePointToCenter = (_points[j] - circleCenter).norm();
				std::cout << "Point " << i << " distance to center: " << distancePointToCenter << std::endl;

				// check if distance is larger than alpha, if so, return false
				if(distancePointToCenter > _alpha + 1e-5) return false;
			}
		}

		std::cout << "------------------------------" << std::endl;
		
		// if all points are at a distance smaller than alpha, return true
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

	void SchCreator2D::makeYAML() {
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
		
		//Alpha value
		out << YAML::Key << "alpha";
		out << YAML::Value << _alpha;

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
		for(size_t i = 0; i <  _schPoints.size(); i++) {
			out << _schPoints[i];
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
	sch.FindSch2D(10.5);
	std::cout << "Is hull strictly convex? " << sch.checkHull() << std::endl;

	return 0;
}