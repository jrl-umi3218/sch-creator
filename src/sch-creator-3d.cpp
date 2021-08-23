#include <sch-creator/sch-creator-3d.h>

namespace sch
{
SchCreator3D::SchCreator3D(double r, double R)
{
  _r = r;
  _R = R;
  _alpha = R - r;
  _epsilon = 1e-8;
} // SchCreator3D

/*
 *   Given a sphere, checks if all points are inside the sphere.
 *   Use this only to check if the algorithm is properly working.
 */
bool SchCreator3D::checkPointsInSphere(const Sphere & s)
{
  double distance, maxDistance = pow(_R, 2) * (1 + _epsilon);
  for(auto i = _vertexes.begin(); i != _vertexes.end(); i++)
  {
    std::cout << " Distance: " << distance << std::endl;

    distance = (s.center - (*i)).squaredNorm();
    if(distance > maxDistance)
    {
      std::cout << std::setprecision(15);
      std::cout << " DOne. Distance: " << distance;
      std::cout << " R: " << maxDistance << std::endl;
      std::cout << std::setprecision(5);
      return false;
    }
  }
  return true;
} // checkPointsInSphere

/*
*   Given a Sphere Center, a vertex B and a radius, gets the derivative 
*   with respect to the radius.
*/
bool SchCreator3D::getDerivative(const SchCreator3D::SphereCenter &currentSphereCenter,
                                 const Eigen::Vector3d &B,
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
} // getDerivative

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

  if(maxD >= 2 * _alpha)
    return true;
  else
    return false;
} // SchCreator3D::findMaxDistance

/*
 *   Builds small spheres from the vertexes obtained with
 *   Polyhedron_algorithm. Radius is equal to _r
 */
void SchCreator3D::getSmallSpheres()
{
  std::cout << "Finding small spheres... ";
  Vector3 temp;
  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    temp = poly.vertexes_[i]->getCoordinates();
    // add vertex to vector
    _vertexes.push_back(Eigen::Vector3d(temp.m_x, temp.m_y, temp.m_z));
    // get small spheres
    _smallSpheres.push_back(Sphere(_vertexes[i], _r));
  }
  std::cout << "Done." << std::endl;
} // SchCreator3D::getSmallSpheres

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
  size_t currTriangle = 0, index = _smallSpheres.size(),
         abID, bcID;
  std::set<size_t> registeredTriangle;
  // bool inSphere = true;

  for(auto i = poly.triangles_.begin(); i != poly.triangles_.end(); i++)
  {
    // get circum sphere
    s = findCircumSphere3((*i).a, (*i).b, (*i).c);
    // add circum radius to heap
    _heap.insert(std::make_pair(s.radius, currTriangle));
    // order sphere vertexes
    abID = getKey((*i).a,(*i).b);
    bcID = getKey((*i).b,(*i).c);
    if((registeredTriangle.count(abID) || registeredTriangle.count(bcID)))
    {
      // get big sphere
      s = findSphereThroughPoints((*i).a, (*i).c, (*i).b);
      bs = BigSphere(s, (*i).a, (*i).c, (*i).b);
      // std::cout << "In sphere: " << checkPointsInSphere(s) << std::endl;
    }
    else
    {
      // get big sphere
      s = findSphereThroughPoints((*i).a, (*i).b, (*i).c);
      bs = BigSphere(s, (*i).a, (*i).b, (*i).c);
      // std::cout << "In sphere: " << checkPointsInSphere(s) << std::endl;
    }
    // register triangle
    registeredTriangle.insert(getKey(bs.p1,bs.p2));
    // registeredTriangle.insert(getEdgeKey(bs.p1,bs.p3));
    registeredTriangle.insert(getKey(bs.p2,bs.p3));
    // add to big sphere vector
    _bigSpheres.push_back(bs);

    // increase the current triangle index
    currTriangle++;
    index++;
  }
  std::cout << " Done." << std::endl;
} // SchCreator3D::getBigSpheres

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

  return Sphere(center, radius);
} // SchCreator3D::findCircumSphere3

/*
*   Given the indexes of four points, finds the circumsphere that
*   touches said indexes.
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
  t << -_vertexes[a].squaredNorm(), -_vertexes[b].squaredNorm(),
       -_vertexes[c].squaredNorm(), -_vertexes[d].squaredNorm();

  // from base matrix get Mn
  Eigen::Matrix4d M1, M2, M3, M4;
  M1 << t, T.rightCols(3);
  M2 << T.col(0), t, T.rightCols(2);
  M3 << T.leftCols(2), t, T.col(3);
  M4 << T.leftCols(3), t;

  // get the determinants
  double Tdet = T.determinant(), D = M1.determinant() / Tdet, 
         E = M2.determinant() / Tdet, F = M3.determinant() / Tdet,
         G = M4.determinant() / Tdet;

  double radius;
  // if(Tdet == 0 || isnan(Tdet)) radius = _epsilon;
  radius = 0.5 * sqrt(pow(D, 2) + pow(E, 2) + pow(F, 2) - 4 * G);
  return Sphere(Eigen::Vector3d(-D / 2, -E / 2, -F / 2), radius);
} // findCircumSphere4


/*
 *   Given 3 points, the function returns the Sphere of radius R
 *   on whose surface lay the points.
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a, size_t b, size_t c)
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
  double sphereRadius = _R-_r;

  // check that sphere's radius is larger than the circle's,
  // else, increase the sphere radius
  while(sphereRadius < circleRadius) sphereRadius += 0.125;

  // get the distance from the center of the sphere
  double distanceFromCenter = sqrt(pow(sphereRadius, 2) - pow(circleRadius, 2));

  // find the center of the sphere
  Eigen::Vector3d sphereCenter = circleCenter3D + distanceFromCenter * n;

  // std::cout << sphereRadius << ' ' << sphereCenter[0] << ' ' << sphereCenter[1] << ' ' << sphereCenter[2] << std::endl;

  return Sphere(sphereCenter, sphereRadius);
} // SchCreator3D::findSphereThroughPoints

/*
 *   Given the index of 3 points, the function finds the plane that touches them
 *   and returns a matrix containing the plane's base vectors.
 */
SchCreator3D::Plane SchCreator3D::findPlaneBase(size_t a, size_t b, size_t c)
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
} // SchCreator3D::findPlaneBase

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
  A << pow(a.norm(), 2), a[0], a[1], 1,
       pow(b.norm(), 2), b[0], b[1], 1, 
       pow(c.norm(), 2), c[0], c[1], 1;

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
} // SchCreator3D::findCircleThroughPoints

/*
 *   Get all torii, its cones and all big sphere plane normals.
 */
void SchCreator3D::getTorii()
{
  size_t id, torusIndex = _smallSpheres.size() + _bigSpheres.size(), index = 0;
  // _bigSphereNormals.resize(_bigSpheres.size());
  std::multimap<size_t, SCHplane> orderedPlanes;
  Eigen::Vector3d planeNormal;
  Face face;
  Torus torus;
  SCHplane plane;

  // go through all edges
  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    // get edge key
    id = getEdgeKey((*i).p1, (*i).p2);
    
    // check if key has already been stored
    if(!toriiKey.count(id))
    {
      // store torus index with its key
      toriiKey.insert(std::make_pair(id, torusIndex));
      // add torus to the vector
      torus = getTorus((*i).s, (*i).p1, (*i).p2);
      _torii.push_back(std::make_pair(torusIndex, torus));
      // compute the cones
      _toriiCones.insert(std::make_pair(torusIndex, getCones((*i).p1, (*i).p2, torus.normal)));
      // increase the torus index
      torusIndex++;
    }

    // get the plane normal
    planeNormal = getPlaneNormal((*i).p1, (*i).p2, (*i).s.center);
    face.plane1 = std::make_pair(toriiKey[id], planeNormal);

    // get edge key
    id = getEdgeKey((*i).p2, (*i).p3);

    // check if key has already been stored
    if(!toriiKey.count(id))
    {

      // store torus index with its key
      toriiKey.insert(std::make_pair(id, torusIndex));
      // add torus to the vector
      torus = getTorus((*i).s, (*i).p2, (*i).p3);
      _torii.push_back(std::make_pair(torusIndex, torus));
      // compute the cones
      _toriiCones.insert(std::make_pair(torusIndex, getCones((*i).p2, (*i).p3, torus.normal)));
      // increase the torus index
      torusIndex++;
    }

    // get the plane normal
    planeNormal = getPlaneNormal((*i).p2, (*i).p3, (*i).s.center);
    face.plane2 = std::make_pair(toriiKey[id], planeNormal);

    // get edge key
    id = getEdgeKey((*i).p3, (*i).p1);
    
    // std::cout << "ID: "<< id << " Estado: " << !toriiKey.count(id) << std::endl;

    // check if key has already been stored
    if(!toriiKey.count(id))
    {
      // store torus index with its key
      toriiKey.insert(std::make_pair(id, torusIndex));
      // add torus to the vector
      torus = getTorus((*i).s, (*i).p3, (*i).p1);
      _torii.push_back(std::make_pair(torusIndex, torus));
      // compute the cones
      _toriiCones.insert(std::make_pair(torusIndex, getCones((*i).p3, (*i).p1, torus.normal)));
      // increase the torus index
      torusIndex++;
    }

    // get the plane normal
    planeNormal = getPlaneNormal((*i).p3, (*i).p1, (*i).s.center);
    face.plane3 = std::make_pair(toriiKey[id], planeNormal);

    // add face to sphere normals
    _bigSphereNormals.push_back(face);
    index++;
  }

  index = _smallSpheres.size();
  for(auto i = _bigSphereNormals.begin(); i != _bigSphereNormals.end(); i++)
  {
    // add all big sphere normals to the map:
    // orderedPlanes<torus index,<triangle index,plane normal>>
    orderedPlanes.insert(std::make_pair((*i).plane1.first, 
                                        std::make_pair(index, 
                                                      (*i).plane1.second)));
    orderedPlanes.insert(std::make_pair((*i).plane2.first, 
                                        std::make_pair(index, 
                                                      (*i).plane2.second)));
    orderedPlanes.insert(std::make_pair((*i).plane3.first, 
                                        std::make_pair(index, 
                                                      (*i).plane3.second)));
    // all keys have exactly two SCHplane elements
    index++;
  }

  for(auto i = toriiKey.begin(); i != toriiKey.end(); i++)
  {
    // get the iterator pointing to the first SCHplane
    auto it = orderedPlanes.equal_range((*i).second);
    // store plane
    plane = std::make_pair((*it.first).second.first,(*it.first).second.second);
    // advance the iterator
    it.first++;
    // add elements to the map
    _toriiPlanes.insert(
        std::make_pair((*it.first).first,
                       // first plane
                       std::make_pair(plane,
                                      // second plane
                                      std::make_pair((*it.first).second.first, 
                                                     (*it.first).second.second))));
  }

  // removeUselessTorii();
}

/*
*   Builds a pair of SCHcones given two vertex indexes and a torus normal
*/
std::pair<SchCreator3D::SCHcone, SchCreator3D::SCHcone> 
          SchCreator3D::getCones(size_t a,size_t b,const Eigen::Vector3d & n)
{
  SCHcone abCone, baCone;
  // compute the cosine
  double cosine = getCosine(a, b);
  // get cone [a,b]
  abCone.first = a;
  abCone.second = Cone(n, cosine);
  // get cone [b,a]
  baCone.first = b;
  baCone.second = Cone(-n, cosine);
  return std::make_pair(abCone, baCone);
}

/*
 *   Computes a Cone's cosine given the indexes of the vertices.
 */
double SchCreator3D::getCosine(size_t a, size_t b)
{
  return (_vertexes[a] - _vertexes[b]).norm() / (2 * _alpha);
}

/*
 *   Gets the normal of a plane given two point indexes and a point.
 */
Eigen::Vector3d SchCreator3D::getPlaneNormal(size_t a, size_t b, Eigen::Vector3d c)
{
  // find the normal to the plane
  Eigen::Vector3d ca = _vertexes[a] - c, cb = _vertexes[b] - c, normal;
  return cb.cross(ca).normalized();
  // else return -ab.cross(bc).normalized();
}

/*
 *   Builds the Torus given a sphere and two points.
 */
SchCreator3D::Torus SchCreator3D::getTorus(const Sphere & s, size_t a, size_t b)
{
  Eigen::Vector3d center = (_vertexes[a] + _vertexes[b]) / 2;
  return Torus(center, (_vertexes[a] - _vertexes[b]).normalized(), 
              (center - s.center).norm()-1e-2);
}

/*
 *   Finds a unique key/ID given two numbers.
 */
size_t SchCreator3D::getEdgeKey(size_t a, size_t b)
{
  return ((a < b) ? (a * _numberOfVertexes + b) : (b * _numberOfVertexes + a));
}

/*
 *   Finds a unique key/ID given two numbers.
 */
size_t SchCreator3D::getKey(size_t a, size_t b)
{
  return (a * _numberOfVertexes + b);
}

/*
*   Finds the neighboring cones to every vertex.
*/
void SchCreator3D::getVertexNeighbours()
{
  size_t id;
  Cone cone;
  std::multimap<size_t, size_t> orderedEdges;
  _vertexNeighbours.resize(_numberOfVertexes);

  for(auto i = poly.edges_.begin(); i != poly.edges_.end(); i++)
  {
    // get edge id
    id = getEdgeKey((*i).a, (*i).b);

    // order edges by vertex
    orderedEdges.insert(std::make_pair((*i).a, toriiKey[id]));
    orderedEdges.insert(std::make_pair((*i).b, toriiKey[id]));
  }

  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    // get all elements in ordered edges with the key i,
    // where i is the index to a vertex
    auto it = orderedEdges.equal_range(i);
    for(auto j = it.first; j != it.second; j++)
    {
      // check if the cone's start point [a,b]
      if((*j).first == _toriiCones[(*j).second].first.first)
        // if it's a, take the first cone from _toriiCones
        cone = _toriiCones[(*j).second].first.second;
      else
        // else, take the second Cone
        cone = _toriiCones[(*j).second].second.second;
      // add the cone to the current vertexe's neighbours
      _vertexNeighbours[i].insert(std::make_pair((*j).second, cone));
    }
  }
} // getVertexNeighbours();

void SchCreator3D::removeUselessTorii()
{
  double d;
  std::vector<size_t> torusToRemove;
  for(auto i = _toriiPlanes.begin(); i != _toriiPlanes.end(); i++)
  {
    d =(_bigSpheres[(*i).second.first.first-_smallSpheres.size()].s.center - 
        _bigSpheres[(*i).second.second.first-_smallSpheres.size()].s.center).squaredNorm();
    // std::cout << (*i).first << ' ' << d << std::endl;
    if(d < 1e-5){
      _removeTorii.insert(std::make_pair((*i).first,false));
      torusToRemove.push_back((*i).first);
    } 
    else _removeTorii.insert(std::make_pair((*i).first,true));
  }

  SCHplane temp;
  for(auto i = torusToRemove.begin(); i != torusToRemove.end(); i++)
  {
    _toriiCones.erase(*i);

    std::cout << _toriiPlanes[*i].first.first << ' ' << _toriiPlanes[*i].second.first << ' ';
    temp = _toriiPlanes[*i].first;
    _toriiPlanes[*i].first = _toriiPlanes[*i].second;
    _toriiPlanes[*i].second = temp;
    std::cout << _toriiPlanes[*i].first.first << ' ' << _toriiPlanes[*i].second.first << std::endl;

  }
}


void SchCreator3D::getHeap()
{
  size_t a, b, c, d;
  std::set<size_t> indexes;
  Sphere s;
  for(auto i = _torii.begin(); i != _torii.end(); i++)
  {
    // get the four vertexes
    a = _toriiCones[(*i).first].first.first;
    b = _toriiCones[(*i).first].second.first;
    c = findVertex(_bigSpheres[_toriiPlanes[(*i).first].first.first-_smallSpheres.size()],a,b);
    d = findVertex(_bigSpheres[_toriiPlanes[(*i).first].second.first-_smallSpheres.size()],a,b);
    // compute the circum sphere
    s = findCircumSphere4(a,b,c,d);
    // insert circumradius to the heap
    _heap.insert(std::make_pair(s.radius,(*i).first));
  }

} // getVertexNeighbours

size_t SchCreator3D::findVertex(const BigSphere &bs, size_t a, size_t b)
{
  bool inSet;
  std::set<size_t> indexes;

  // add known vertexes to the set
  indexes.insert(a);
  indexes.insert(b);

  // add first sphere vertex
  inSet = indexes.insert(bs.p1).second;
  // if p1 is not in the set, return p1
  if(inSet) return bs.p1;
  // else, try adding p2 to the set
  else inSet = indexes.insert(bs.p2).second;

  // if p2 isn't in the set, return p2
  if(inSet) return bs.p2;
  // else, return p3
  else return bs.p3;
} // insertToSet

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
  os << ' ' << _smallSpheres.size() << std::endl;
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
      os << (*j).first << ' ' << (*j).second.cosangle << ' ';
      os << (*j).second.axis[0] << ' ';
      os << (*j).second.axis[1] << ' ';
      os << (*j).second.axis[2] << std::endl;
    }
    index++;
  }

  // number of big spheres
  index = 0;
  os << _bigSpheres.size() << std::endl;
  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    // radius and center of big sphere
    os << _R << ' ' << (*i).s.center[0] << ' ';
    os << (*i).s.center[1] << ' ' << (*i).s.center[2] << std::endl;

    // coordinates of the three points touching the sphere
    os << _vertexes[(*i).p1][0] << ' ' << _vertexes[(*i).p1][1] << ' ' << _vertexes[(*i).p1][2] << ' ';
    os << _vertexes[(*i).p2][0] << ' ' << _vertexes[(*i).p2][1] << ' ' << _vertexes[(*i).p2][2] << ' ';
    os << _vertexes[(*i).p3][0] << ' ' << _vertexes[(*i).p3][1] << ' ' << _vertexes[(*i).p3][2] << std::endl;

    // normals to the sphere planes
    os << _bigSphereNormals[index].plane1.first << ' ';
    os << _bigSphereNormals[index].plane1.second[0] << ' ';
    os << _bigSphereNormals[index].plane1.second[1] << ' ';
    os << _bigSphereNormals[index].plane1.second[2] << std::endl;

    os << _bigSphereNormals[index].plane2.first << ' ';
    os << _bigSphereNormals[index].plane2.second[0] << ' ';
    os << _bigSphereNormals[index].plane2.second[1] << ' ';
    os << _bigSphereNormals[index].plane2.second[2] << std::endl;

    os << _bigSphereNormals[index].plane3.first << ' ';
    os << _bigSphereNormals[index].plane3.second[0] << ' ';
    os << _bigSphereNormals[index].plane3.second[1] << ' ';
    os << _bigSphereNormals[index].plane3.second[2] << std::endl;

    index++;
  }

  // number of torii
  os << _torii.size() << std::endl;
  for(auto i = _torii.begin(); i != _torii.end(); i++)
  {
    // remove torus boolean
    os << 1 << std::endl;

    // write torus' external radius, center and normal coordinates
    os << (*i).second.extRadius << ' ' << _R << ' ';

    os << (*i).second.center[0] << ' ';
    os << (*i).second.center[1] << ' ';
    os << (*i).second.center[2] << ' ';

    os << (*i).second.normal[0] << ' ';
    os << (*i).second.normal[1] << ' ';
    os << (*i).second.normal[2] << std::endl;

    // write torus' cones
    // cone [a,b]: a index, cosangle, axis coordinates
    os << _toriiCones[(*i).first].first.first << ' ';
    os << _toriiCones[(*i).first].first.second.cosangle << ' ';
    os << _toriiCones[(*i).first].first.second.axis[0] << ' ';
    os << _toriiCones[(*i).first].first.second.axis[1] << ' ';
    os << _toriiCones[(*i).first].first.second.axis[2] << std::endl;

    // cone [b,a]: b index, cosangle, axis coordinates
    os << _toriiCones[(*i).first].second.first << ' ';
    os << _toriiCones[(*i).first].second.second.cosangle << ' ';
    os << _toriiCones[(*i).first].second.second.axis[0] << ' ';
    os << _toriiCones[(*i).first].second.second.axis[1] << ' ';
    os << _toriiCones[(*i).first].second.second.axis[2] << std::endl;
    
    // write torus' planes
    // first triangle: triangle index, normal coordinates
    os << _toriiPlanes[(*i).first].first.first << ' ';
    os << _toriiPlanes[(*i).first].first.second[0] << ' ';
    os << _toriiPlanes[(*i).first].first.second[1] << ' ';
    os << _toriiPlanes[(*i).first].first.second[2] << std::endl;

    // second triangle: triangle index, normal coordinates
    os << _toriiPlanes[(*i).first].second.first << ' ';
    os << _toriiPlanes[(*i).first].second.second[0] << ' ';
    os << _toriiPlanes[(*i).first].second.second[1] << ' ';
    os << _toriiPlanes[(*i).first].second.second[2] << std::endl;
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
  getTorii();
  getVertexNeighbours();
  getHeap();

  std::map<double, size_t>::iterator it = _heap.begin();
  double maxHeap = (*it).first;
  std::cout << "Max Heap: " << maxHeap << std::endl;

  while(maxHeap > _R)
  {
    _heap.erase(it);
    it = _heap.begin();
    maxHeap = (*it).first;
    std::cout << "Max Heap: " << maxHeap << std::endl;
  }
}
} // namespace sch

// ./build/src/sch-creator-3d /home/amrf/balloon-inflating/sch-visualization/tests/shared-tests/data/sample_polyhedron.otp 
int main(int argc, char ** argv)
{
  namespace po = boost::program_options;

  double r, R;

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")(
      "r,r", po::value<double>(&r)->default_value(.02),
      "small sphere radius")("R,R", po::value<double>(&R)->default_value(300.),
                             "big sphere radius")("input-file", po::value<std::string>(),
                                                  "input file")("output-file", po::value<std::string>(), "output file");

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
    sch::SchCreator3D sch(r, R);
    sch.computeSCH(input);
    sch.writeToFile(output);

    // /home/amrf/balloon-inflating/sch-visualization/tests/shared-tests/data/sample_polyhedron.otp
  }

  return 0;
}