#include <sch-creator/sch-creator-2d.h>

namespace sch
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
		std::cout << "Reading points... ";
		std::string fileString;
		std::ifstream file;
		double x, y;

		// open the file
		file.open(_pointsPath);
		//verify file is correctly opened, exit otherwise
		if(!file) {
			std::cout << "\nFile not found. Exiting..." << std::endl << std::endl;
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
			// add sch::Point to vector
			_pointsStructure.push_back(sch::SchCreator2D::Point(Eigen::Vector2d(x,y)));
		}

		//std::cout << "***** Initial points *****" << std::endl;
		//printPoints();
		
		// close the file
		file.close();
		std::cout << "Done." << std::endl;
	}

	// print triangles
	void SchCreator2D::printPoints()
	{
		std::cout << "----- Points -----" << std::endl;
		int n = 0;

		for (auto it = _pointsStructure.begin(); it != _pointsStructure.end(); it++)
		{
			std::cout << "Point " << n << std::endl;
			std::cout << "inHull: " << (*it).inHull << std::endl;
			std:: cout << (*it).point << '\n' << std::endl;
			n++;
		}
	}

	void SchCreator2D::printTriangles()
	{
		std::cout << "----- Triangles -----" << std::endl;
		int n = 0;

		for (auto it = _triangles.begin(); it != _triangles.end(); it++)
		{
			std::cout << "Triangle " << n << std::endl;
			std::cout << "inHeap: " << (*it).inHeap << std::endl;
			std::cout << "Start Point:" << (*it).a[0] << ' ' << (*it).a[1] << std::endl;
			std::cout << "Midpoint Point:" << (*it).b[0] << ' ' << (*it).b[1] << std::endl;
			std::cout << "Endpoint Point:" << (*it).c[0] << ' ' << (*it).c[1] << std::endl;
			std:: cout << "Circumcircle radius: " << (*it).r << '\n' << std::endl;
			n++;
		}
	}

	/* Listing triangles and inserting their circumcircle radius to the heap */
	void SchCreator2D::listTriangles(size_t n) 
	{
		std::cout << "Finding initial triangles... ";
		Triangle t;
		Heap h;

		// declare the unnordered vector of Heap
		std::vector<SchCreator2D::Heap> radii(n); 

		for (size_t i = 0; i < n; i++) {
			// make triangle and corresponding Heap
			t = Triangle(_pointsStructure[i % n].point, 
						 _pointsStructure[(i + 1) % n].point, 
						 _pointsStructure[(i + 2) % n].point);
			h = Heap(i % n, (i + 1) % n, (i + 2) % n,
				    (i + n - 1) % n, (i + n) % n, (i + 1 + n) % n, t.r);

			// std::cout << "Midpoint: " << (i + 1) % n << std::endl;
			// std::cout << "Prev Triangle: " << (i + n - 1) % n << std::endl;
			// std::cout << "Curr Triangle: " << (i + n) % n << std::endl;
			// std::cout << "Next Triangle: " << (i + 1 + n) % n << std::endl;

			// add triangle to list
			_triangles.push_back(t);
			// add radius to unordered vector
			radii[i] = h;
		}

		// std::cout << "***** Initial triangles *****" << std::endl;
		//printTriangles();

		// create max heap of Heap 
		_heap = std::priority_queue<SchCreator2D::Heap>(radii.begin(), radii.end());

		std::cout << "Done." << std::endl;
	}

	/* 	changes the in Hull state of the removed point and changes its corresponding
		triangle's inHeap state to false */
	void SchCreator2D::removePointFromHull(const Heap & maxHeap) 
	{
		// create iterator for the list of triangles
		std::list<Triangle>::iterator it;
		// get triangles size, index and amount of points
		size_t trianglesSize = _triangles.size();
		size_t trianglesIndex[3] = { maxHeap.prevTriangleIndex,
									 maxHeap.triangleIndex,
									 maxHeap.nextTriangleIndex };
		for(int i = 0; i < 3; i++)
		{
			// initialize iterator to the begining of the list
			it = _triangles.begin();
			// advance iterator to point at triangle
			advance(it, (trianglesIndex[i] + trianglesSize) % trianglesSize);
			// set the inHeap property of the triangle to false
			(*it).removeFromHeap();
		}
		
		

		// set the inHull property of the corresponding midpoint to false
		_pointsStructure[maxHeap.midpointIndex].removeFromHull();

		
		// std::cout << "Point removed" << std::endl;
		// std::cout << "Circumcircle: " << maxHeap.radius << std::endl;
		// std::cout << "Prev Triangle: " << (trianglesIndex[0] + trianglesSize) % trianglesSize << std::endl;
		// std::cout << "Curr Triangle: " << (trianglesIndex[1] + trianglesSize) % trianglesSize << std::endl;
		// std::cout << "Next Triangle: " << (trianglesIndex[2] + trianglesSize) % trianglesSize << std::endl;
		
		
	}

	/* 	Verifies all points corresponding to the max heap exist in the sch */
	bool SchCreator2D::checkIfMaxHeapIsInHull()
	{
		// initialize iterator at the begining of the list
		std::list<Triangle>::iterator it = _triangles.begin();

		// move the iterator to point at the triangle corresponding to the max heap
		// advance(it, (_heap.top()).prevTriangleIndex);
		// bool prevTriangleInHeap = (*it).inHeap;

		it = _triangles.begin();
		advance(it, (_heap.top()).triangleIndex);
		bool currTriangleInHeap = (*it).inHeap;

		// it = _triangles.begin();
		// advance(it, (_heap.top()).endpointIndex);
		// bool nextTraingleInHeap = (*it).inHeap;
		
		// check if all vertices of the triangle are in the Hull (inHull == true)
		bool trianlgePointsInHull = _pointsStructure[(_heap.top()).startpointIndex].inHull 
									&& (_pointsStructure[(_heap.top()).midpointIndex].inHull 
									&& _pointsStructure[(_heap.top()).endpointIndex].inHull);

		// true only when all vertices are in the Hull and the triangle is in the heap
		return currTriangleInHeap && trianlgePointsInHull;
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

	/* Finds the closest previous triangle */
	size_t SchCreator2D::findPreviousTriangle(size_t triangleIndex)
	{
		std::list<Triangle>::iterator it;
		size_t i, n = _triangles.size();
		for(i = 1; i < n-1; i++)
		{
			it = _triangles.begin();
			advance(it, (triangleIndex - i + n) % n);
			if ((*it).inHeap) break;
		}

		return (triangleIndex - i + n) % n;
	}

	// find next triangle in heap
	size_t SchCreator2D::findNextTriangle(size_t triangleIndex)
	{
		std::list<Triangle>::iterator it;
		size_t i, n = _triangles.size();
		for(i = 1; i < n; i++)
		{
			it = _triangles.begin();
			advance(it, (triangleIndex + i + n) % n);
			if ((*it).inHeap) break;
		}

		return (triangleIndex + i + n) % n;
	}

	/* Makes two new triangles based on the point that must be deleted from the sch */
	void SchCreator2D::makeTriangles(size_t previousMidpoint) {
		size_t n = _pointsStructure.size();

		// get the three points of the first new triangle
		size_t newMidpoint = findPreviousPoint(n + previousMidpoint - 1) % n;
		size_t newStartpoint = findPreviousPoint(n + newMidpoint - 1) % n;
		size_t newEndpoint = findNextPoint(previousMidpoint + 1) % n; 

		// get the previous triangle in the heap
		size_t prevTriangle = findPreviousTriangle(_heap.top().prevTriangleIndex);

		// build the first new triangle and radius
		Triangle newTriangle1 = Triangle(_pointsStructure[newStartpoint].point,
										_pointsStructure[newMidpoint].point, 
										_pointsStructure[newEndpoint].point);
		Heap newRadius1 = Heap(newStartpoint, newMidpoint, newEndpoint, 
										prevTriangle, _triangles.size(), 
										_triangles.size()+1, newTriangle1.r);

		// get the three points of the second new triangle
		newStartpoint = newMidpoint % n;
		newMidpoint = newEndpoint % n;
		newEndpoint = findNextPoint(newMidpoint + 1) % n;

		// get the new prev and next triangle
		prevTriangle = _triangles.size();
		size_t nextTriangle = findNextTriangle(_heap.top().nextTriangleIndex);

		// build the second triangle and radius
		Triangle newTriangle2 = Triangle(_pointsStructure[newStartpoint].point, 
										_pointsStructure[newMidpoint].point, 
										_pointsStructure[newEndpoint].point);
		Heap newRadius2 = Heap(newStartpoint, newMidpoint, newEndpoint, 
										prevTriangle, prevTriangle + 1,
										nextTriangle, newTriangle2.r);
		
		// add the triangles to the list 
		_triangles.push_back(newTriangle1);
		_triangles.push_back(newTriangle2);
		// add radius to the heap
		_heap.push(newRadius1);
		_heap.push(newRadius2);
	}

	/* 	Stores the index of the max heap midpoint, removes it from the heap 
		and gets the two new triangles*/
	void SchCreator2D::updateTriangles()
	{
		// get the index of the point to eliminate
		size_t eliminatedPointIndex = (_heap.top()).midpointIndex;
		// remove max heap
		_eliminatedVertex.push_back(_heap.top());
		// make the mew triangles
		makeTriangles(eliminatedPointIndex);

		_heap.pop();
		
	}

	void SchCreator2D::findNewAlpha()
	{
		std::cout << "Computing new value for alpha... ";
		_initialAlpha = _alpha;
		double topAlpha = (_heap.top()).radius + 0.05;

		while(true)
		{
			if((_heap.top()).radius < _alpha) break;
			else
			{
				_alpha += (topAlpha - _alpha) /2;
				
			}
			
		}
		_alpha -= 0.05;
		std::cout << "Done. New Alpha: " << _alpha << std::endl;	
	}

	void SchCreator2D::findClosestPoint()
	{
		std::priority_queue<Heap> newHeap = _heap;
		Heap r1 = newHeap.top();
		size_t point1Index = r1.midpointIndex;
		newHeap.pop();
		Heap r2 = newHeap.top();
		size_t point2Index = r2.midpointIndex;
		newHeap.pop();
		Heap r3 = newHeap.top();
		size_t point3Index = r3.midpointIndex;

		Eigen::Vector2d p1 = _pointsStructure[point1Index].point;
		Eigen::Vector2d p2 = _pointsStructure[point2Index].point;
		Eigen::Vector2d p3 = _pointsStructure[point3Index].point;

		double d12 = (p1-p2).norm();
		double d31 = (p1-p3).norm();
		double d32 = (p2-p3).norm();

		if (d12 > d31)
		{
			if (d31 > d32) 
			{
				_pointsStructure[point2Index].removeFromHull();
				_eliminatedVertex.push_back(r2);
			}
			else 
			{
				_pointsStructure[point3Index].removeFromHull();
				_eliminatedVertex.push_back(r3);
			}
		} else
		{
			if (d12 > d32)
			{
				_pointsStructure[point3Index].removeFromHull();
				_eliminatedVertex.push_back(r3);
			}
			else 
			{
				_pointsStructure[point1Index].removeFromHull();
				_eliminatedVertex.push_back(r1);
			}
		}
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
			std::cout << "The convex hull must have at least 3 points. Exiting..." << std::endl;
			exit(1);
		} 

		// make list of all initial triangles
		listTriangles(pointsInSCH);

		std::cout << "Finding SCH..." << std::endl;
		// while the max Heap radius is larger than alpha
		while((_heap.top()).radius > _alpha) {
			// Remove middle point and make the new triangles
			removePointFromHull(_heap.top());
			
			updateTriangles();
			// decrease the amount of points in the sch
			pointsInSCH--;
			// ensure that, after elimination, there's still at least 3 points
			if(pointsInSCH == 3) {
				std::cout << "Alpha is too small." << std::endl;
				findNewAlpha();
				findClosestPoint();
				pointsInSCH--;
				break;
				// exit(1);
			}

			// verify that all points in the new max heap belong to the hull,
			// eliminate the max heap otherwise
			while(!checkIfMaxHeapIsInHull()) {
				_heap.pop();
			}

			// std::cout << '\n' << "New heap: " << _heap.top().radius << '\n' << std::endl;
			// std::cout << "-------------------------" << std::endl;
		}

		//printTriangles();
		std::cout << "SCH successfully found. " << std::endl;
		// printPoints();

		// add all points that are still part of the hull to the schPoints array
		for (auto i = _pointsStructure.begin(); i != _pointsStructure.end(); i++) {
			if((*i).inHull){
				_schPoints.push_back((*i).point);
			}
		}

		// get the amount of eliminated points
		_eliminatedPoints = _points.size()-_schPoints.size();

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
				std::cout << "Point " << j << " distance to center: " << distancePointToCenter << std::endl;

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

	// Overwriting the << operation for sch's Heap and the YAML emitter
	YAML::Emitter& operator << (YAML::Emitter& out, const sch::SchCreator2D::Heap& r) {
		out << YAML::Flow;
		out << YAML::BeginSeq << r.radius << r.midpointIndex << YAML::EndSeq;
		// Heap' radius value and eliminated point index added to a YAML flow
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
		out << YAML::Comment("[Heap, index]");
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

int main(int argc, char** argv) {
	std::cout << "\n SCH parameters: a = " <<  argv[2] << std::endl << std::endl;
    std::cout << "Opening " <<  argv[1] << std::endl;

	sch::SchCreator2D sch(argv[1]);
	sch.FindSch2D(std::stod(argv[2]));

	if(argc == 4)
		std::cout << "Is hull strictly convex? " << sch.checkHull() << std::endl;

	return 0;
}