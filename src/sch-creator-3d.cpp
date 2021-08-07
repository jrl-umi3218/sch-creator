// #include <Eigen/Dense>
// #include <iostream>
// #include <vector>
#include <sch-creator/sch-creator-3d.h>

namespace sch
{
SchCreator3D::SchCreator3D(double r, double R) 
{
  _r = r;
  _R = R;
  _alpha = R-r;
  _epsilon = 1e-5;
}

std::ostream & operator<<(std::ostream & os, const SchCreator3D::Sphere & s)
{
  os << "Center: \t" << s.center[0] << ' ' << s.center[1] << ' ' << s.center[2];
  os << "\nRadius: " << s.radius << '\n';

  return os;
}

/*
*   Prints every vertex along with its neighboring cones
*/
void SchCreator3D::printVertexNeighbours()
{
  std::cout << std::endl;
  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    int index = 0;
    std::cout << "***** VERTEX " << i << " *****"<< std::endl;
    std::cout << "INDEX" << '\t' << "COSANGLE" << "\t" << "AXIS COORDINATES" << std::endl;
    for(auto j = _vertexNeighbours[i].begin(); j != _vertexNeighbours[i].end(); j++)
    {
      std::cout << (*j) << '\t' << _cones[*j].cosangle << '\t';
      std::cout << _cones[*j].axis[0] << ' ' << _cones[*j].axis[1];
      std::cout << ' ' << _cones[*j].axis[2] << std::endl;
      index++;
    }
    std::cout << std::endl;
  }
}

/*
*   Prints every big sphere along with its plane normals
*/
void SchCreator3D::printBigSpherePlanes()
{
  for(size_t i = 0; i < _bigSpheres.size(); i++)
  {
    std::cout << "***** SPHERE " << i << " *****"<< std::endl;
    std::cout << _bigSpheres[i].s;
    std::cout << "----- Normals -----" << std::endl;
    for(auto j = _bigSphereNormals[i].begin(); j != _bigSphereNormals[i].end(); j++)
    {
      std::cout << (*j)[0] << ' ' << (*j)[1] << ' ' << (*j)[2] << std::endl;
    }
    std::cout << std::endl;
  }

}

/*
*   Prints the index of every vertex on the big sphere along with its edges
*/
void SchCreator3D::printSphereEdges()
{
  for(size_t i = 0; i < _bigSpheres.size(); i++)
  {
    std::cout << "Sphere " << i << std::endl;
    std::cout << _bigSpheres[i].p1 << ' ' << _bigSpheres[i].p2 << ' ' << _bigSpheres[i].p3 << std::endl;
    for(auto j = _bigSphereEdgess[i].begin(); j != _bigSphereEdgess[i].end(); j++)
    {
      std::cout << "Edge: ";
      std::cout << poly.edges_[(*j)].a << ' ' << poly.edges_[(*j)].b << std::endl;
    }
  }
}

/*
 *   Given 3 points, the function returns the circumsphere
 */
SchCreator3D::Sphere SchCreator3D::findCircumSphere3(size_t a, size_t b, size_t c)
{
  // find vectors a to b and a to c
  Eigen::Vector3d ab, ac;
  ab = _vertexes[b] - _vertexes[a];
  ac = _vertexes[c] - _vertexes[a];

  // compute dot products
  double v11 = ab.dot(ab);
  double v22 = ac.dot(ac);
  double v12 = ab.dot(ac);

  double num = v11 * v22 - pow(v12, 2);
  double k1 = 0.5 * v22 * (v11 - v12) / num;
  double k2 = 0.5 * v11 * (v22 - v12) / num;

  Eigen::Vector3d center = _vertexes[a] + k1 * ab + k2 * ac;
  double radius = (k1 * ab + k2 * ac).norm();

  return Sphere(center,radius);
}

/*
 *   Given the index of 3 points, the function finds the plane that touches them
 *   and returns a matrix containing the plane's base vectors.
 */
SchCreator3D::Plane SchCreator3D::findPlaneBase(size_t a,
                                                size_t b,
                                                size_t c)
{
  // find the normal to the plane
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], bc = _vertexes[c] - _vertexes[b];
  Eigen::Vector3d n = ab.cross(bc).normalized();

  // find base vector ex and ey
  Eigen::Vector3d ex = ab.cross(n).normalized();
  Eigen::Vector3d ey = n.cross(ex);
  Eigen::MatrixXd base(2, 3);
  base << ex.transpose(), ey.transpose();

  return Plane(base, n, a, b, c);
}

/*
 *   Given the index of 2 points and a point, the function finds the plane
 *   that touches them and returns a matrix containing the plane's base vectors.
 */
SchCreator3D::Plane SchCreator3D::findPlaneBase(size_t a,
                                                size_t b,
                                                const Eigen::Vector3d &c)
{
  // find the normal to the plane
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], bc = c - _vertexes[b];
  Eigen::Vector3d n = ab.cross(bc).normalized();

  // find base vector ex and ey
  Eigen::Vector3d ex = ab.cross(n).normalized();
  Eigen::Vector3d ey = n.cross(ex);
  Eigen::MatrixXd base(2, 3);
  base << ex.transpose(), ey.transpose();

  return Plane(base, n);
}

/*
 *   Given 3 points, the function returns the center of the circumcircle.
 *
 *   NOTE: Points must be 2D.
 */
Eigen::Vector2d SchCreator3D::findCircleThroughPoints(const Eigen::Vector2d & a,
                                                      const Eigen::Vector2d & b,
                                                      const Eigen::Vector2d & c)
{
  // Using crammer's rule, find the base matrix
  Eigen::MatrixXd A(3, 4);
  A << pow(a.norm(), 2), a[0], a[1], 1, pow(b.norm(), 2), b[0], b[1], 1, pow(c.norm(), 2), c[0], c[1], 1;

  // From the base matrix get M1n
  Eigen::Matrix3d M11 = A.block(0, 1, 3, 3), M12, M13, M14;
  M12 << A.leftCols(1), A.rightCols(2);
  M13 << A.leftCols(2), A.rightCols(1);
  M14 << A.leftCols(3);

  // compute the center's x and y coordinates
  Eigen::Vector2d center;
  center(0) = 0.5 * (M12.determinant() / M11.determinant());
  center(1) = -0.5 * (M13.determinant() / M11.determinant());

  return center;
}

/*
 *   Given 3 points, the function returns the Sphere of radius R
 *   on whose surface lay the points.
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a,
                                                           size_t b,
                                                           size_t c)
{
  // get plane base and normal
  SchCreator3D::Plane p = findPlaneBase(a, b, c);
  Eigen::MatrixXd base = p.base;
  Eigen::Vector3d n = p.normal;

  // find the 2D projection of the points using the base
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], ac = _vertexes[c] - _vertexes[a];
  Eigen::Vector2d a2d(0, 0), b2d = base * ab, c2d = base * ac;

  // get the center of the circle (in 2D)
  Eigen::Vector2d circleCenter2D = findCircleThroughPoints(a2d, b2d, c2d);

  // convert coordinates to 3D
  Eigen::Vector3d circleCenter3D;
  circleCenter3D(0) = _vertexes[a](0) + circleCenter2D.dot(base.col(0));
  circleCenter3D(1) = _vertexes[a](1) + circleCenter2D.dot(base.col(1));
  circleCenter3D(2) = _vertexes[a](2) + circleCenter2D.dot(base.col(2));

  // get sphere and circle radius
  double circleRadius = circleCenter2D.norm();
  double sphereRadius = _R;

  // check that sphere's radius is larger than the circle's,
  // else, increase the sphere radius
  while(sphereRadius < circleRadius) sphereRadius += 0.125;

  // get the distance from the center of the sphere
  double distanceFromCenter = sqrt(pow(sphereRadius, 2) - pow(circleRadius, 2));

  // find the center of the sphere
  Eigen::Vector3d sphereCenter = circleCenter3D + distanceFromCenter * n;

  // add center values to vector
  _sphereCenters.push_back(SphereCenter(circleRadius, circleCenter3D, n));

  return Sphere(sphereCenter, sphereRadius);
}

/*
 *   Given 3 points and a radius, the function returns the Sphere
 *   on whose surface lay the points.
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a,
                                                           size_t b,
                                                           size_t c,
                                                           double radius)
{
  // get plane base and normal
  SchCreator3D::Plane p = findPlaneBase(a, b, c);
  Eigen::MatrixXd base = p.base;
  Eigen::Vector3d n = p.normal;

  // find the 2D projection of the points using the base
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], ac = _vertexes[c] - _vertexes[a];
  Eigen::Vector2d a2d(0, 0), b2d = base * ab, c2d = base * ac;

  // get the center of the circle (in 2D)
  Eigen::Vector2d circleCenter2D = findCircleThroughPoints(a2d, b2d, c2d);

  // convert coordinates to 3D
  Eigen::Vector3d circleCenter3D;
  circleCenter3D(0) = _vertexes[a](0) + circleCenter2D.dot(base.col(0));
  circleCenter3D(1) = _vertexes[a](1) + circleCenter2D.dot(base.col(1));
  circleCenter3D(2) = _vertexes[a](2) + circleCenter2D.dot(base.col(2));

  // get sphere and circle radius
  double circleRadius = circleCenter2D.norm();
  double sphereRadius = radius;

  // check that sphere's radius is larger than the circle's,
  // else, increase the sphere radius
  while(sphereRadius < circleRadius) sphereRadius += 0.125;

  // get the distance from the center of the sphere
  double distanceFromCenter = sqrt(pow(sphereRadius, 2) - pow(circleRadius, 2));

  // find the center of the sphere
  Eigen::Vector3d sphereCenter = circleCenter3D + distanceFromCenter * n;

  // add center values to vector
  _sphereCenters.push_back(SphereCenter(circleRadius, circleCenter3D, n));

  return Sphere(sphereCenter, sphereRadius);
}

/*
 *   Given 4 points, the function returns the Sphere
 *   on whose surface lay the points.
 */
SchCreator3D::Sphere SchCreator3D::findCircumSphere4(size_t a,
                                                     size_t b,
                                                     size_t c,
                                                     size_t d)
{
  // ussing crammer's rule find base matrix
  Eigen::Matrix4d T;
  T << _vertexes[a][0], _vertexes[a][1], _vertexes[a][2], 1,
       _vertexes[b][0], _vertexes[b][1], _vertexes[b][2], 1, 
       _vertexes[c][0], _vertexes[c][1], _vertexes[c][2], 1, 
       _vertexes[d][0], _vertexes[d][1], _vertexes[d][2], 1;

  // get squared norm values
  Eigen::Vector4d t;
  t << -pow(_vertexes[a].norm(), 2), -pow(_vertexes[b].norm(), 2),
       -pow(_vertexes[c].norm(), 2), -pow(_vertexes[d].norm(), 2);

  // from base matrix get Mn
  Eigen::Matrix4d M1, M2, M3, M4;
  M1 << t, T.rightCols(3);
  M2 << T.col(0), t, T.rightCols(2);
  M3 << T.leftCols(2), t, T.col(3);
  M4 << T.leftCols(3), t;

  // get the determinants
  double Tdet = T.determinant(), D = M1.determinant() / Tdet, E = M2.determinant() / Tdet, F = M3.determinant() / Tdet,
         G = M4.determinant() / Tdet;

  return Sphere(Eigen::Vector3d(-D / 2, -E / 2, -F / 2), 0.5 * sqrt(pow(D, 2) + pow(E, 2) + pow(F, 2) - 4 * G));
}

bool SchCreator3D::getDerivative(const SchCreator3D::SphereCenter & currentSphereCenter,
                                 const Eigen::Vector3d B,
                                 double radius)
{
  double circleRadius = currentSphereCenter.circleRadius;
  double distanceToCenter = sqrt(pow(radius, 2) - pow(circleRadius, 2));

  Eigen::Vector3d n = currentSphereCenter.planeNormal;
  Eigen::Vector3d c = currentSphereCenter.circleCenter;

  double diff = n[0] * ((c[0] + distanceToCenter * n[0]) - B[0]) 
                + n[1] * ((c[1] + distanceToCenter * n[1]) - B[1])
                + n[2] * ((c[2] + distanceToCenter * n[2]) - B[2]);

  return ((2 * radius / distanceToCenter) * diff - 2 * radius) > 0;
}

/*
*   Checks edges to find the largest distance between vertexes.
*   Returns true if it's larger than two times alpha.
*/
bool SchCreator3D::findMaxDistance()
{
  double maxD = 0, d;

  for(size_t i = 0; i < poly.edges_.size(); i++)
  {
    // get the length of current vertex
    d = poly.edges_[i].edge.norm();
    // update lardest distance
    if(maxD < d) maxD = d;
  }

  std::cout << "Maximum body distance: " << maxD << std::endl;

  if(maxD >= 2*_alpha) return true;
  return false;
}

/*
*   Builds small spheres from the vertexes obtained with
*   Polyhedron_algorithm. Radius is equal to _r
*/
void SchCreator3D::getSmallSpheres()
{
  Vector3 temp;
  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    temp = poly.vertexes_[i]->getCoordinates();
    // add vertex to vector
    _vertexes.push_back(Eigen::Vector3d(temp.m_x, temp.m_y, temp.m_z));
    // get small spheres
    _smallSpheres.push_back(Sphere(_vertexes[i],_r));
  }
}

/*
*   Builds the big spheres from the triangles obtained with
*   Polyhedron_algorithms and adds them to a vector.
*   Adds the circum radius and the triangle index to the heap.
*/
void SchCreator3D::getBigSpheres()
{
  std::cout << "Finding big spheres...";
  Sphere s;
  BigSphere bs;
  size_t currTriangle = 0;

  for(auto i = poly.triangles_.begin(); i != poly.triangles_.end(); i++)
  {
    // get circum sphere
    s = findCircumSphere3((*i).a,(*i).b,(*i).c);
    // add circum radius to heap
    heap_.insert(std::make_pair(s.radius,currTriangle));
    // get big sphere
    s = findSphereThroughPoints((*i).a,(*i).b,(*i).c);
    bs = BigSphere(s,(*i).c,(*i).b,(*i).a);
    // add to big sphere vector
    _bigSpheres.push_back(bs);
    
    // increase the current triangle index
    currTriangle++;
  }
  std::cout << " Found " << _bigSpheres.size() << " big spheres." << std::endl;
}

void SchCreator3D::getBigSpherePlanes()
{
  std::cout << "Finding planes...";

  Plane p;
  _bigSphereNormals.resize(_bigSpheres.size());

  for(size_t i = 0; i < _bigSpheres.size(); i++)
  {
    p = findPlaneBase(_bigSpheres[i].p1,_bigSpheres[i].p2,_bigSpheres[i].s.center);
    _bigSphereNormals[i].push_back(p.normal);

    p = findPlaneBase(_bigSpheres[i].p1,_bigSpheres[i].p3,_bigSpheres[i].s.center);
    _bigSphereNormals[i].push_back(p.normal);

    p = findPlaneBase(_bigSpheres[i].p2,_bigSpheres[i].p3,_bigSpheres[i].s.center);
    _bigSphereNormals[i].push_back(p.normal);
  }

  std::cout << " Done." << std::endl;
}

/*
*   Finds a unique key/ID given two numbers.
*/
size_t SchCreator3D::getEdgeKey(size_t a, size_t b)
{
  return ((a<b)?(a*_numberOfVertexes+b):(b*_numberOfVertexes+a));
}


/*
*     Finds the three neighboring edges to a big sphere.
*/
void SchCreator3D::getBigSphereEdges()
{
  size_t index = 0;
  std::map<size_t,size_t> orderedEdges;
  std::map<size_t,size_t>::iterator it;
  _bigSphereEdgess.resize(_bigSpheres.size());
  
  for(size_t i = 0; i < poly.edges_.size(); i++)
  {
    orderedEdges.insert(std::make_pair(getEdgeKey(poly.edges_[i].a,
                                       poly.edges_[i].b),i));
  }

  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    it = orderedEdges.find(getEdgeKey((*i).p1, (*i).p2));
    _bigSphereEdgess[index].push_back((*it).second);

    it = orderedEdges.find(getEdgeKey((*i).p2, (*i).p3));
    _bigSphereEdgess[index].push_back((*it).second);

    it = orderedEdges.find(getEdgeKey((*i).p1, (*i).p3));
    _bigSphereEdgess[index].push_back((*it).second);

    index++;
  }
}


/*
*   Finds the neighbpours of an edge (2 vertexes and 2 triangles), finds the circumsphere
*   made up by the four neighbours and adds the radius to the heap.
*/
void SchCreator3D::getEdgeNeighbours()
{
  size_t index = 0,a;
  std::pair<size_t,size_t> vertexes, triangles;
  std::vector<std::pair<std::pair<size_t,size_t>,std::pair<size_t,size_t>>> neighbours;
  std::multimap<size_t,size_t> orderedTriangles;
  std::set<size_t> orderedVertexes;

  // get all triangles ordered by edges
  for(auto i = poly.triangles_.begin(); i != poly.triangles_.end(); i++)
  {
    orderedTriangles.insert(std::make_pair(getEdgeKey((*i).a,(*i).b), index));
    orderedTriangles.insert(std::make_pair(getEdgeKey((*i).b,(*i).c), index));
    orderedTriangles.insert(std::make_pair(getEdgeKey((*i).a,(*i).c), index));
    index++;
  }

  for(auto i = poly.edges_.begin(); i != poly.edges_.end(); i++)
  {
    // find the neighboring triangles
    auto it = orderedTriangles.equal_range(getEdgeKey((*i).a, (*i).b));
    for(auto j = it.first; j != it.second; j++)
    {
      // get all neighboring vertexes
      orderedVertexes.insert(poly.triangles_[(*j).second].a);
      orderedVertexes.insert(poly.triangles_[(*j).second].b);
      orderedVertexes.insert(poly.triangles_[(*j).second].c);
    }

    // remove known vertexes
    orderedVertexes.erase((*i).a);
    orderedVertexes.erase((*i).b);

    // get remaining vertexes
    auto iterator = orderedVertexes.begin();
    a = (*iterator);
    orderedVertexes.erase(iterator);
    iterator = orderedVertexes.begin();

    // get all four neighboring vertexes
    neighbours.push_back(std::make_pair(std::make_pair((*i).a, (*i).b),
                         std::make_pair(a,(*iterator))));
    orderedVertexes.clear();
  }

  Sphere s;
  index = 0;
  for(auto i = neighbours.begin(); i != neighbours.end(); i++)
  {
    // find 4 vertex circumsphere
    s = findCircumSphere4((*i).first.first,(*i).first.second,(*i).second.first,(*i).second.second);
    // add radius to heap
    heap_.insert(std::make_pair(s.radius, index));
    index++;
  }

  // index=0;
  // for(auto i = _cones.begin(); i != _cones.end(); i++)
  // {
  //   std::cout << "Cone " << index << std::endl;
  //   std::cout << "Neighbouring vertexes: " << neighbours[index].first.first << ' ' << neighbours[index].first.second;
  //   std::cout << ' ' << neighbours[index].second.first << ' ' << neighbours[index].second.second << std::endl;
  //   std::cout << std::endl;
  //   index++;
  // }
}

/*
*   Gets the neighboring edges to each small sphere.
*/
void SchCreator3D::getCones()
{
  std::cout << "Finding neighbours...";

  double cosangle;
  size_t index = 0;
  Eigen::Vector3d edge;
  std::multimap<size_t,size_t> orderedEdges;
  _vertexNeighbours.resize(_numberOfVertexes);

  for(auto i = poly.edges_.begin(); i != poly.edges_.end(); i++)
  {
    // get the toruses
    cosangle = (*i).edge.norm() / (2*_alpha);
    edge = Eigen::Vector3d((*i).edge.m_x,(*i).edge.m_y,(*i).edge.m_z);
    _cones.push_back(Cone(edge,cosangle));

    // order toruses by vertex
    orderedEdges.insert(std::make_pair((*i).a, index));
    orderedEdges.insert(std::make_pair((*i).b, index));
    index++;
  }

  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    auto it = orderedEdges.equal_range(i);
    for(auto j = it.first; j != it.second; j++)
    {
      _vertexNeighbours[i].push_back((*j).second);
    }
  }

  std::cout << " Done." << std::endl;
}

SchCreator3D::Torus SchCreator3D::getTorus(const Sphere &s, size_t a, size_t b)
{
  Eigen::Vector3d center = (_vertexes[a] + _vertexes[b]) / 2; 
  return Torus(center,
              (_vertexes[a] - _vertexes[b]).normalized(),
              (center - s.center).norm());
}

/*
*   Gets the torii from the big spheres
*/
void SchCreator3D::getTorii()
{
  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    _torii.push_back(getTorus((*i).s,(*i).p1,(*i).p2)); 
    _torii.push_back(getTorus((*i).s,(*i).p1,(*i).p3)); 
    _torii.push_back(getTorus((*i).s,(*i).p2,(*i).p3)); 
  }
}

/*
*   Creates the output file given a filename.
*/
void SchCreator3D::writeToFile(const std::string & filename)
{
  // write data into a file
  std::ofstream os;
  os.open(filename.c_str());
  os.precision(10);
  size_t index = 0;

  // radii used for computation on first line
  os << _r << ' ' << _R << std::endl;

  // number of small spheres
  os << ' ' << _smallSpheres.size() << std:: endl;
  // writes small spheres
  for(auto i = _smallSpheres.begin(); i != _smallSpheres.end(); i++)
  {
    os << (*i).radius << ' ' << (*i).center[0] << ' ';
    os << (*i).center[1] << ' ' << (*i).center[2] << std::endl;

    // number of neighbours
    os << _vertexNeighbours[index].size() << std::endl;
    for(auto j = _vertexNeighbours[index].begin(); j != _vertexNeighbours[index].end(); j++)
    {
      // neighbours
      os << (*j) << ' ' << _cones[*j].cosangle << ' ' << _cones[*j].axis[0];
      os << ' ' << _cones[*j].axis[1] << ' ' << _cones[*j].axis[2] << std::endl; 
    }
    index++;
  }

  // number of big spheres
  index = 0;
  os << _bigSpheres.size() << std::endl;
  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    // radius and center of big sphere
    os << (*i).s.radius << ' ' << (*i).s.center[0] << ' ';
    os << (*i).s.center[1] << ' ' << (*i).s.center[2] << std::endl;

    // coordinates of the three points touching the sphere
    os << _vertexes[(*i).p1][0] << ' ' << _vertexes[(*i).p1][1] << ' ' << _vertexes[(*i).p1][2] << ' ';
    os << _vertexes[(*i).p2][0] << ' ' << _vertexes[(*i).p2][1] << ' ' << _vertexes[(*i).p2][2] << ' ';
    os << _vertexes[(*i).p3][0] << ' ' << _vertexes[(*i).p3][1] << ' ' << _vertexes[(*i).p3][2] << std::endl;

    // normals to the sphere planes
    for(auto j = _bigSphereNormals[index].begin(); j != _bigSphereNormals[index].end(); j++)
    {
      os << (*j)[0] << ' ' << (*j)[1] << ' ' << (*j)[2] << std::endl;
    } 

    index++;
  }

  // number of torii
  index = 0;
  os << _torii.size() << std::endl;

  for(auto i = _torii.begin(); i != _torii.end(); i++)
  {
    // write torus external radius, center and normal
    os << (*i).extRadius << ' ' << _R << ' ';
    os << (*i).center[0] << ' ' << (*i).center[1] << ' ' << (*i).center[2] << ' ';
    os << (*i).normal[0] << ' ' << (*i).normal[1] << ' ' << (*i).normal[2] << std::endl;
  }
  
  // close the file
  os.close();
}

void SchCreator3D::computeSCH(const std::string & filename)
{
  // Read points from file with poly_algorithms
  poly.openFromFile(filename);

  // get no. of vertex
  _numberOfVertexes = poly.vertexes_.size();
  _vertexes.reserve(_numberOfVertexes);
  // std::cout << "No. Vertex: " << _numberOfVertexes << std::endl;
  // std::cout << "No. Triangles: " << poly.triangles_.size() << std::endl;
  // std::cout << "No. Edges: " << poly.edges_.size() << std::endl;

  // Ensure alpha is smaller than the max. body distance,
  if(findMaxDistance())
  {
    std::cout << "ERROR, impossible to compute SCH, choose larger R." << std::endl;
    std::cout << "Finished." << std::endl;
    return;
  }

  // get spheres
  getSmallSpheres();
  getBigSpheres();
  getBigSpherePlanes();
  getBigSphereEdges();
  getCones();
  getEdgeNeighbours();
  getTorii();

  std::map<double,size_t>::iterator it = heap_.begin();
  double maxHeap = (*it).first;
  std::cout << "Max Heap: " << maxHeap << std::endl;

  while(maxHeap > _R)
  {
    heap_.erase(it);
    it = heap_.begin();
    maxHeap = (*it).first;
    std::cout << "Max Heap: " << maxHeap << std::endl;
  }
}

} // namespace sch

// Run:
// "/home/amrf/balloon-inflating/sch-creator/build/src/sch-creator-3d"

int main(int argc, char ** argv)
{
  namespace po = boost::program_options;

  double r, R;

  po::options_description desc("Allowed options");
  desc.add_options()
  ("help,h", "produce help message")
  ("r,r", po::value<double>(&r)->default_value(.02),"small sphere radius")
  ("R,R", po::value<double>(&R)->default_value(300.),"big sphere radius")
  ("input-file", po::value<std::string>(),"input file")
  ("output-file", po::value<std::string>(),"output file");

  po::positional_options_description pos;
  pos.add("input-file", 1);
  pos.add("output-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  po::notify(vm);

  if(vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 1;
  }

  std::cout << "\n SCH parameters: r = " << r << ", R = " << R << std::endl << std::endl;

  if(vm.count("input-file") && vm.count("output-file"))
  {
    std::fstream testfile;
    std::string input = vm["input-file"].as<std::string>();
    std::string output = vm["output-file"].as<std::string>();

    std::cout << "Opening " << input << std::endl;
    sch::SchCreator3D sch(r,R);
    sch.computeSCH(input);
    sch.writeToFile(output);
      
// /home/amrf/balloon-inflating/sch-visualization/tests/shared-tests/data/sample_polyhedron.otp 
  }

  return 0;
}