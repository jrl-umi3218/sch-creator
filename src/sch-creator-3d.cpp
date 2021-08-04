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
  os << "----- Sphere -----\n";
  os << "Center: \n" << s.center;
  os << "\nRadius: " << s.radius << '\n';

  return os;
}

/*
*   Prints every vertex along with its neighboring cones
*/
void SchCreator3D::printVertexNeighbours()
{
  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    int index = 0;
    std::cout << "Vertex " << i << std::endl;
    std::cout << "N" << '\t' << "COSANGLE" << "\t\t" << "AXIS COORDINATES" << std::endl;
    for(size_t j = 0; j < _vertexNeighbours[i].size(); j++)
    {
      std::cout << index << '\t' << _vertexNeighbours[i][j].cosangle << '\t';
      std::cout << '\t' << _vertexNeighbours[i][j].axis[0] << '\t' << _vertexNeighbours[i][j].axis[1];
      std::cout << '\t' << _vertexNeighbours[i][j].axis[2] << std::endl;
      index++;
    }
    std::cout << std::endl;
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

  return Plane(base, n);
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
SchCreator3D::Sphere SchCreator3D::findCircumSphere4(const Eigen::Vector3d & a,
                                                     const Eigen::Vector3d & b,
                                                     const Eigen::Vector3d & c,
                                                     const Eigen::Vector3d & d)
{
  // ussing crammer's rule find base matrix
  Eigen::Matrix4d T;
  T << a[0], a[1], a[2], 1, b[0], b[1], b[2], 1, c[0], c[1], c[2], 1, d[0], d[1], d[2], 1;

  // get squared norm values
  Eigen::Vector4d t;
  t << -pow(a.norm(), 2), -pow(b.norm(), 2), -pow(c.norm(), 2), -pow(d.norm(), 2);

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
    // get big sphere
    bs = BigSphere(s,(*i).c,(*i).b,(*i).a);
    // add to big sphere vector
    _bigSpheres.push_back(bs);
    // add to heap
    heap_.insert(std::make_pair(s.radius,currTriangle));
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
*   Gets the neighboring edges to each small sphere.
*/
void SchCreator3D::getCones()
{
  std::cout << "Finding neighbours...";

  Cone edge;
  double cosangle;
  std::multimap<size_t,Cone> orderedEdges;
  _vertexNeighbours.resize(_numberOfVertexes);

  // order edges by vertexes
  for(auto i = poly.edges_.begin(); i != poly.edges_.end(); i++)
  {
    cosangle = (*i).edge.norm() / (2*_alpha);
    edge = Cone(Eigen::Vector3d((*i).edge.m_x,(*i).edge.m_y,(*i).edge.m_z), cosangle);
    orderedEdges.insert(std::make_pair((*i).a, edge));
    orderedEdges.insert(std::make_pair((*i).b, edge));
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

/*
*   Creates the output file given a filename.
*/
void SchCreator3D::writeToFile(const std::string & filename)
{
  // write data into a file
  std::ofstream os;
  os.open(filename.c_str());
  os.precision(16);
  size_t index = 0;

  // radii used for computation on first line
  os << _r << ' ' << _R << std::endl;

  // number of small spheres
  os << ' ' << _smallSpheres.size() << std:: endl;
  // writes small spheres
  for(auto i = _smallSpheres.begin(); i != _smallSpheres.end(); i++)
  {
    std::cout << (*i).radius << ' ' << (*i).center[0] << ' ';
    std::cout << (*i).center[1] << ' ' << (*i).center[2] << std::endl;

    // number of neighbours
    std::cout << _vertexNeighbours[index].size() << std::endl;
    for(auto j = _vertexNeighbours[index].begin(); j != _vertexNeighbours[index].end(); j++)
    {
      // neighbours
      std::cout << (*j).cosangle << ' ' << (*j).axis[0];
      std::cout << ' ' << (*j).axis[1] << ' ' << (*j).axis[2]; 
    }
  }
}

void SchCreator3D::computeSCH(const std::string & filename)
{
  // Read points from file with poly_algorithms
  poly.openFromFile(filename);

  // get no. of vertex
  _numberOfVertexes = poly.vertexes_.size();
  _vertexes.reserve(_numberOfVertexes);
  std::cout << "No. Vertex: " << _numberOfVertexes << std::endl;
  std::cout << "No. Triangles: " << poly.triangles_.size() << std::endl;
  std::cout << "No. Edges: " << poly.edges_.size() << std::endl;

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
  getCones();

  std::multimap<double,size_t>::iterator it = heap_.begin();
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

    

      
// /home/amrf/balloon-inflating/sch-visualization/tests/shared-tests/data/sample_polyhedron.otp 
  }

  return 0;
}