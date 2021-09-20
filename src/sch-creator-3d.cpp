#include <sch-creator/sch-creator-3d.h>

namespace sch
{
SchCreator3D::SchCreator3D(double r, double R)
{
  _r = r;
  _R = R;
  _desiredAlpha = R - r;
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
    // std::cout << " Distance: " << distance << std::endl;

    distance = (s.center - (*i)).squaredNorm();
    if(distance - maxDistance > _epsilon)
    {
      std::cout << std::setprecision(15);
      std::cout << " Done. Distance: " << distance;
      std::cout << " R: " << maxDistance << std::endl;
      std::cout << std::setprecision(5);
      return false;
    }
  }
  return true;
} // checkPointsInSphere

/*
*   Given the indexes of the vertexes of the testing sphere and the dissapearing testing point,
*   gets the derivative with respect to the radius.
*/
bool SchCreator3D::getDerivative(size_t v1, size_t v2, size_t v3, size_t v4)
{
  // get the plane through the points
  SchCreator3D::Plane p = findPlaneBase(v1, v3, v2);

  // get the circle on the plane that goes thorugh the points
  Eigen::Vector2d circle2D = findCircleThroughPoints(v1, v3, v2,p);

  // get the thistance from the center of the circel to the sphere center
  double dToSphereCenter = sqrt(_alpha*_alpha - circle2D.squaredNorm());

  // convert the circle's center coordinates to 3D
  Eigen::Vector3d circle3D;
  circle3D(0) = _vertexes[v1](0) + circle2D.dot(p.base.col(0));
  circle3D(1) = _vertexes[v1](1) + circle2D.dot(p.base.col(1));
  circle3D(2) = _vertexes[v1](2) + circle2D.dot(p.base.col(2));

  // get the derivative of the sphere center coordinates
  Eigen::Vector3d dsc = (_alpha/dToSphereCenter)*p.normal;

  return 2*((circle3D+dToSphereCenter*p.normal).dot(dsc)-_vertexes[v4].dot(dsc)-_alpha) > 0;
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

  if(maxD >= 2 * _desiredAlpha)
    return true;
  else
    return false;
} // SchCreator3D::findMaxDistance

void SchCreator3D::initialize()
{
  Vector3 vertexTemp;
  size_t id, index = 0, torusIndex = 0, e1, e2, e3;
  std::vector<size_t> tempX;
  std::multimap<size_t, size_t> orderedEdges;
  std::map<size_t, size_t> edgesByKey;
  std::multimap<size_t, size_t> orderedFaces;

  for(auto i = poly.triangles_.begin(); i != poly.triangles_.end(); i++)
  {
    id = getEdgeKey((*i).a,(*i).b);
    e1 = id;
    orderedFaces.insert(std::make_pair(id,index));
    if(edgesByKey.insert(std::make_pair(id,torusIndex)).second)
      torusIndex++;

    id = getEdgeKey((*i).b,(*i).c);
    e2 = id;
    orderedFaces.insert(std::make_pair(id,index));
    if(edgesByKey.insert(std::make_pair(id,torusIndex)).second)
      torusIndex++;

    id = getEdgeKey((*i).c,(*i).a);
    e3 = id;
    orderedFaces.insert(std::make_pair(id,index));
    if(edgesByKey.insert(std::make_pair(id,torusIndex)).second)
      torusIndex++;


    _SCHtriangles.push_back(SCHtriangle((*i).a,(*i).b,(*i).c,
                            edgesByKey[e1],edgesByKey[e2],edgesByKey[e3]));

    index++;
  }

  // index = 0;
  // std::cout << "\nTriangles" << std::endl;
  // for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  // {
  //   std::cout << index << ' ' << (*i).p1 << ' ' << (*i).p2
  //             << ' ' << (*i).p3 << std::endl;
  //   index++;
  // }

  index = 0;
  for(auto i = poly.edges_.begin(); i != poly.edges_.end(); i++)
  {
    // get edge id
    id = getEdgeKey((*i).a, (*i).b);

    // order edges by vertex
    orderedEdges.insert(std::make_pair((*i).a,edgesByKey[id]));
    orderedEdges.insert(std::make_pair((*i).b,edgesByKey[id]));


    auto it = orderedFaces.equal_range(id);
    for(auto j = it.first; j != it.second; j++)
      tempX.push_back((*j).second);

    _SCHedges.push_back(SCHedge((*i).a,(*i).b,tempX[0],tempX[1]));

    index++;
    tempX.clear();
  }

  // index = 0;
  // std::cout << "\nEdge neighbours" << std::endl;
  // for(auto i = _SCHedges.begin(); i != _SCHedges.end(); i++)
  // {
  //   std::cout << index << ' ' << (*i).vertex1 << ' ' << (*i).vertex2
  //             << ' ' << (*i).face1 << ' ' << (*i).face2 << std::endl;
  //   index++;
  // }

  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    vertexTemp = poly.vertexes_[i]->getCoordinates();
    // add vertex to vector
    _vertexes.push_back(Eigen::Vector3d(vertexTemp.m_x,
                                        vertexTemp.m_y,
                                        vertexTemp.m_z));

    // get all neighbours for the current vertex
    auto neighbours = orderedEdges.equal_range(i);
    // add them to a vector
    for(auto j = neighbours.first; j!= neighbours.second; j++)
      tempX.push_back((*j).second);
    // add SCHvertex to vector
    _SCHvertexes.push_back(SCHvertex(_vertexes[i],tempX));
    tempX.clear();
  }

  // std::cout << "\nVertex neighbours" << std::endl;
  // index = 0;
  // for(auto i = _SCHvertexes.begin(); i != _SCHvertexes.end(); i++)
  // {
  //   std::cout << index << ' ';
  //   for(auto j = (*i).neighbours.begin(); j != (*i).neighbours.end(); j++)
  //     std::cout << *j << ' ';
  //   std::cout << std::endl;
  //   index++;
  // }

  // find the center of the point cloud
  initialCenter = Eigen::Vector3d(0,0,0);
  for(auto i = _vertexes.begin(); i != _vertexes.end(); i++)
  {
    initialCenter += *i;
  }
  initialCenter/=_numberOfVertexes;
} // initialize

/*
 *   Builds small spheres from the vertexes obtained with
 *   Polyhedron_algorithm. Radius is equal to _r
 */
void SchCreator3D::getSmallSpheres()
{
  std::cout << "Finding small spheres... ";
  for(auto i = _SCHvertexes.begin(); i != _SCHvertexes.end(); i++)
  {
    if((*i).inHull)
      // get small spheres
      _smallSpheres.push_back(Sphere((*i).vertex,_r));
  }
  std::cout << _smallSpheres.size() << ' ';
  std::cout << "Done." << std::endl;
} // getSmallSpheres

void SchCreator3D::updateVertexesIndex()
{
  size_t index;
  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    index = i;
    for(auto j = difference.begin(); j != difference.end(); j++)
      if(index > *j) index--;
      else break;
    
    v.push_back(index);
  }

  // std::cout << "O" << '\t' << "v[i]" << std::endl;
  // for(size_t i = 0; i < _numberOfVertexes; i++)
  //   std::cout << i << '\t' << v[i] << std::endl;
}

/*
 *   Builds the big spheres from the triangles obtained with
 *   Polyhedron_algorithms and adds them to a vector.
 *   Adds the circum radius and the triangle index to the heap.
 */
void SchCreator3D::getBigSpheres()
{
  std::cout << "Finding big spheres..." << std::endl;
  Sphere s;
  BigSphere bs;
  double diff;
  size_t index = _smallSpheres.size();

  for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  {
    if(!(*i).inHull) continue;

    s = findSphereThroughPoints((*i).p1, (*i).p2, (*i).p3);
    
    diff = (s.center-initialCenter).squaredNorm()-(_desiredAlpha*_desiredAlpha*(1-_epsilon));
    // std::cout << index << ' ' << diff;
    // std::cout << ' ' << pow(diff,2) << ' ';
    // std::cout << checkPointsInSphere(s) << std::endl;
    if(diff > 0 || pow(diff,2) < 1e-5)
    {
      s = findSphereThroughPoints((*i).p1, (*i).p3, (*i).p2);
      // std::cout << "Switch sphere: ";
      // std::cout << index << ' ' << diff << std::endl;
      // std::cout << checkPointsInSphere(s) << std::endl;

      bs = BigSphere(s,v[(*i).p1],v[(*i).p3],v[(*i).p2]);
    } else
    {
      bs = BigSphere(s,v[(*i).p1],v[(*i).p2],v[(*i).p3]);
    }

    // add to big sphere vector
    _bigSpheres.push_back(bs);

    // increase the current triangle index
    index++;
  }
  std::cout << _bigSpheres.size();
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
  double Tdet = T.determinant();

  if(Tdet == 0) Tdet = 1e-10;
  
  double D = M1.determinant() / Tdet, E = M2.determinant() / Tdet, 
         F = M3.determinant() / Tdet, G = M4.determinant() / Tdet;

  double radius = 0.5 * sqrt(pow(D, 2) + pow(E, 2) + pow(F, 2) - 4 * G);
  return Sphere(Eigen::Vector3d(-D / 2, -E / 2, -F / 2), radius);
} // findCircumSphere4


/*
 *   Given 3 points, the function returns the Sphere of radius R
 *   on whose surface lay the points.
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a, size_t b, size_t c)
{
  // get plane base and normal
  SchCreator3D::Plane p = findPlaneBase(a, c, b);
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
  // double circleRadius = circleCenter2D.norm();
  double sphereRadius = _desiredAlpha;

  // check that sphere's radius is larger than the circle's,
  // else, increase the sphere radius
  // while(sphereRadius < circleRadius) sphereRadius += 0.125;

  // get the distance from the center of the sphere
  double distanceFromCenter = sqrt(pow(sphereRadius, 2) - circleCenter2D.squaredNorm());

  // find the center of the sphere
  Eigen::Vector3d sphereCenter = circleCenter3D + (distanceFromCenter + _epsilon) * n;

  // std::cout << sphereRadius << ' ' << sphereCenter[0] << ' ' << sphereCenter[1] << ' ' << sphereCenter[2] << std::endl;

  return Sphere(sphereCenter, _R);
} // SchCreator3D::findSphereThroughPoints

SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a,SchCreator3D::Plane p,
                                                           Eigen::Vector2d circleCenter2D)
{
  return findSphereThroughPoints(a,p,circleCenter2D,_alpha);
} // SchCreator3D::findSphereThroughPoints

SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a,SchCreator3D::Plane p,
                                                           Eigen::Vector2d circleCenter2D,
                                                           double R)
{
  // convert coordinates to 3D
  Eigen::Vector3d circleCenter3D;
  circleCenter3D(0) = _vertexes[a](0) + circleCenter2D.dot(p.base.col(0));
  circleCenter3D(1) = _vertexes[a](1) + circleCenter2D.dot(p.base.col(1));
  circleCenter3D(2) = _vertexes[a](2) + circleCenter2D.dot(p.base.col(2));

  // get sphere and circle radius
  // double circleRadius = circleCenter2D.norm();
  double sphereRadius = R;

  // check that sphere's radius is larger than the circle's,
  // else, increase the sphere radius
  // while(sphereRadius < circleRadius) sphereRadius += 0.125;

  // get the distance from the center of the sphere
  double distanceFromCenter = sqrt(pow(sphereRadius, 2) - circleCenter2D.squaredNorm());

  // find the center of the sphere
  Eigen::Vector3d sphereCenter = circleCenter3D + (distanceFromCenter + _epsilon) * p.normal;

  // std::cout << sphereRadius << ' ' << sphereCenter[0] << ' ' << sphereCenter[1] << ' ' << sphereCenter[2] << std::endl;

  return Sphere(sphereCenter, R);
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

Eigen::Vector2d SchCreator3D::findCircleThroughPoints(size_t a,
                                                      size_t b,
                                                      size_t c,
                                                      SchCreator3D::Plane p)
{
  // find the 2D projection of the points using the base
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], ac = _vertexes[c] - _vertexes[a];
  Eigen::Vector2d a2d(0, 0), b2d = p.base * ab, c2d = p.base * ac;
  findCircleThroughPoints(a2d, b2d, c2d);
  // Using crammer's rule, find the base matrix
  Eigen::MatrixXd A(3, 4);
  A << pow(a2d.norm(), 2), a2d[0], a2d[1], 1,
       pow(b2d.norm(), 2), b2d[0], b2d[1], 1,
       pow(c2d.norm(), 2), c2d[0], c2d[1], 1;

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
  std::cout << "Finding torii..." << std::endl;
  size_t id, torusIndex = _smallSpheres.size() + _bigSpheres.size(),
         index = _smallSpheres.size();
  std::multimap<size_t, SCHplane> orderedPlanes;
  Eigen::Vector3d planeNormal;
  Face face;
  Torus torus;
  SCHplane plane;

  // go through all edges
  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    // std::cout << (*i).p1 << ' ' << (*i).p2 << ' ' << (*i).p3 << std::endl;
    // get edge key
    id = getEdgeKey((*i).p1, (*i).p2);

    // check if key has already been stored
    if(!toriiKey.count(id))
    {
      // store torus index with its key
      toriiKey.insert(std::make_pair(id, torusIndex));
      // add torus to the vector
      torus = getTorus((*i).s, (*i).p1, (*i).p2);
      _torii.push_back(torus);
      // compute the cones
      _toriiCones.insert(std::make_pair(torusIndex,getCones((*i).p1,(*i).p2,torus.normal)));
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
      // std::cout << torusIndex << '\t' << (*i).p2 << '\t' << (*i).p3 << std::endl;
      _torii.push_back(torus);
      // compute the cones
      _toriiCones.insert(std::make_pair(torusIndex,getCones((*i).p2,(*i).p3, torus.normal)));
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
      // std::cout << torusIndex << '\t' << (*i).p3 << '\t' << (*i).p1 << std::endl;
      _torii.push_back(torus);
      // compute the cones
      _toriiCones.insert(std::make_pair(torusIndex,getCones((*i).p3,(*i).p1,torus.normal)));
      // increase the torus index
      torusIndex++;
    }

    // get the plane normal
    planeNormal = getPlaneNormal((*i).p3, (*i).p1, (*i).s.center);
    face.plane3 = std::make_pair(toriiKey[id], planeNormal);

    // add face to sphere normals
    _bigSphereNormals.insert(std::make_pair(index,face));
    index++;
  }

  index = _smallSpheres.size();
  for(auto i = _bigSphereNormals.begin(); i != _bigSphereNormals.end(); i++)
  {
    // add all big sphere normals to the map:
    // orderedPlanes<torus index,<triangle index,plane normal>>
    orderedPlanes.insert(std::make_pair((*i).second.plane1.first,
                                        std::make_pair(index,
                                                      (*i).second.plane1.second)));
    orderedPlanes.insert(std::make_pair((*i).second.plane2.first,
                                        std::make_pair(index,
                                                      (*i).second.plane2.second)));
    orderedPlanes.insert(std::make_pair((*i).second.plane3.first,
                                        std::make_pair(index,
                                                      (*i).second.plane3.second)));
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

  // std::cout << "Torii cones " << std::endl;
  // for(auto i = _toriiCones.begin(); i != _toriiCones.end(); i++)
  //   std::cout << (*i).first << ' ' << (*i).second.first.first << ' ' << (*i).second.second.first << std::endl;
  
  // std::cout << "Sphere normals " << std::endl;
  // for(auto i = _bigSphereNormals.begin(); i != _bigSphereNormals.end(); i++)
  // {
  //   std::cout << (*i).first << ' ' << (*i).second.plane1.first << ' '
  //             << (*i).second.plane2.first << ' '
  //             << (*i).second.plane3.first << std::endl;
  // }
  std::cout << " Done." << std::endl;
} // getTorii

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
  return (_smallSpheres[a].center - _smallSpheres[b].center).norm() / (2 * _desiredAlpha);
}

/*
 *   Gets the normal of a plane given two point indexes and a point.
 */
Eigen::Vector3d SchCreator3D::getPlaneNormal(size_t a, size_t b, Eigen::Vector3d c)
{
  // find the normal to the plane
  Eigen::Vector3d ca = _smallSpheres[a].center - c, 
                  cb = _smallSpheres[b].center - c;
  return -cb.cross(ca).normalized();
  // else return -ab.cross(bc).normalized();
}

/*
 *   Builds the Torus given a sphere and two points.
 */
SchCreator3D::Torus SchCreator3D::getTorus(const Sphere & s, size_t a, size_t b)
{
  Eigen::Vector3d center = (_smallSpheres[a].center + _smallSpheres[b].center) / 2;
  return Torus(center, (_smallSpheres[a].center - _smallSpheres[b].center).normalized(),
              (center - s.center).norm());
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
  size_t id, index = 0;
  Cone cone;
  std::multimap<size_t, size_t> orderedEdges;
  _vertexNeighbours.resize(_smallSpheres.size());

  for(auto i = _SCHedges.begin(); i != _SCHedges.end(); i++)
  {
    if(!(*i).inHull) continue;
    // get edge id
    id = getEdgeKey(v[(*i).vertex1], v[(*i).vertex2]);

    // order edges by vertex
    orderedEdges.insert(std::make_pair(v[(*i).vertex1], toriiKey[id]));
    orderedEdges.insert(std::make_pair(v[(*i).vertex2], toriiKey[id]));

    // std::cout << toriiKey[id] << ' ' << v[(*i).vertex1] << ' ' << v[(*i).vertex2] << std::endl;
  }

  for(auto i = _smallSpheres.begin(); i != _smallSpheres.end(); i++)
  {
    // get all elements in ordered edges with the key i,
    // where i is the index to a vertex
    auto it = orderedEdges.equal_range(index);
    // std::cout << index << std::endl;
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
      _vertexNeighbours[index].push_back(std::make_pair((*j).second, cone));
      // std::cout << ' ' << (*j).second;
    }

    // std::cout << std::endl;

    index++;
  }
} // getVertexNeighbours();

void SchCreator3D::getHeap()
{
  std::set<size_t> indexes;
  Sphere s;
  std::vector<SCHheap> tempHeap;

  size_t index = 0;
  for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  {
    // get circum sphere
    s = findCircumSphere3((*i).p1, (*i).p2, (*i).p3);
    tempHeap.push_back(SCHheap(s.radius,type(triangle),index));
    temp.push_back(SCHheap(s.radius,type(triangle),index));
    std::cout << (*i).p1 << ' ' << (*i).p2 << ' ' << (*i).p3 << ' ' << s.radius << std::endl;

    index++;
  }

  index = 0;
  for(auto i = _SCHedges.begin(); i != _SCHedges.end(); i++)
  {
    // insert circumradius to the heap
    tempHeap.push_back(SCHheap(getEdgeHeap(*i),type(edge),index));
    temp.push_back(SCHheap(getEdgeHeap(*i),type(edge),index));
    index++;
  }

  _heap = std::priority_queue<SCHheap>(tempHeap.begin(),tempHeap.end());
} // getHeap

double SchCreator3D::getEdgeHeap(const SchCreator3D::SCHedge &e)
{
  size_t vertex3, vertex4;
  Sphere s;
  // get the four vertexes
  vertex3 = findVertex(_SCHtriangles[e.face1],e.vertex1,e.vertex2);
  vertex4 = findVertex(_SCHtriangles[e.face2],e.vertex1,e.vertex2);
  // compute the circum sphere
  s = findCircumSphere4(e.vertex1,e.vertex2,vertex3,vertex4);
  std::cout << e.vertex1 << ' ' << e.vertex2 << ' ' << vertex3 << ' ' << vertex4
            << ' ' << s.radius << std::endl;
  return s.radius;
} // getEdgeHeap

size_t SchCreator3D::findVertex(const SCHtriangle &t,size_t a,size_t b)
{
  bool inSet;
  std::set<size_t> indexes;

  // add known vertexes to the set
  indexes.insert(a);
  indexes.insert(b);

  // add first sphere vertex
  inSet = indexes.insert(t.p1).second;
  // if p1 is not in the set, return p1
  if(inSet) return t.p1;
  // else, try adding p2 to the set
  else inSet = indexes.insert(t.p2).second;
  // if p2 isn't in the set, return p2
  if(inSet) return t.p2;
  // else, return p3
  else return t.p3;
} // findVertex

void SchCreator3D::changeTopology()
{
  size_t f1, f2;
  size_t a, b, c, d, i;
  SCHneighbours fg, hj;

  a = _SCHedges[_heap.top().index].vertex1;
  b = _SCHedges[_heap.top().index].vertex2;
  f1 = _SCHedges[_heap.top().index].face1;
  f2 = _SCHedges[_heap.top().index].face2;
  c = findVertex(f1,a,b);
  d = findVertex(f2,a,b);

  if(c == d)
  {
    std::cout << "Same triangle" << std::endl;
  }

  fg = findEdge(a, c, d);
  hj = findEdge(b, c, d);

  std::cout << "RELEVANT DATA" << std::endl;
  std::cout << a << ' ' << b << ' ' << c << ' ' << d << std::endl;
  std::cout << fg.first << ' ' << fg.second << ' ' << hj.first << ' ' << hj.second << std::endl;

  std::cout << "Check if same neighbour: "; //<< checkSameNeighbour(c,d,type(edge)) << std::endl;
  if(!checkSameNeighbour(c,d,type(edge)))
  {
    std::cout << "Invert Edge\n" << std::endl;
    invertEdge(c,d,fg,hj);

  } else
  {
    i = findEdge(c,d);
    std::cout << "Check which vertex dissapears" << std::endl;

    std::cout << "Check i: " << i << std::endl;;
    if(checkSameNeighbour(fg.second,i,type(triangle)))
    {
      std::cout << "Active neighbours: " << findActiveNeighbours(a) << std::endl;
      std::cout << "Check if a Dissapears under: " << std::endl;
      bool dissapear = getDerivative(c,d,b,a);
      std::cout << dissapear << std::endl;

      if(!dissapear)
      {
        std::cout << "Invert Edge\n" << std::endl;
        invertEdge(c,d,fg,hj);
      } else
      {
        std::cout << "Dissapear vertex\n" << std::endl;
        dissapearVertex(a,c,d,b,fg,hj,i);
      }
    } else
    {
      std::cout << "Active neighbours: " << findActiveNeighbours(a) << std::endl;
      std::cout << "Check if b Dissapears under " << std::endl;
      bool dissapear = getDerivative(c,d,a,b);
      std::cout << dissapear << std::endl;
      if(!dissapear)
      {
        std::cout << "Invert Edge\n" << std::endl;
        invertEdge(c,d,fg,hj);

      } else
      {
        std::cout << "Dissapear vertex\n" << std::endl;
        dissapearVertex(b,c,d,a,hj,fg,i);
      }
    }
  }
} // changeTopology

size_t SchCreator3D::findActiveNeighbours(size_t v)
{
  size_t activeNeighbours = 0;
  for(auto i = _SCHvertexes[v].neighbours.begin(); 
      i !=  _SCHvertexes[v].neighbours.end(); i++)
  {
    if(_SCHedges[*i].inHull) activeNeighbours++;
  }
  return activeNeighbours;
}

void SchCreator3D::substituteByVertex()
{ 
} // substituteByVertex

void SchCreator3D::dissapearVertex(size_t v1, size_t v2, size_t v3 , size_t v4,
                                        SCHneighbours e12,SCHneighbours e34, size_t e)
{
  // insert dissapearing vertex to queue
  difference.insert(v1);
  // remove elements from hull
  _SCHvertexes[v1].removeFromHull();
  _SCHedges[_heap.top().index].removeFromHull();
  _SCHedges[e12.first].removeFromHull();
  _SCHedges[e12.second].removeFromHull();
  _SCHtriangles[_SCHedges[_heap.top().index].face1].removeFromHull();
  _SCHtriangles[_SCHedges[_heap.top().index].face2].removeFromHull();

  // remove common face
  if(!_SCHtriangles[_SCHedges[e12.first].face1].inHull)
    _SCHtriangles[_SCHedges[e12.first].face2].removeFromHull();
  else _SCHtriangles[_SCHedges[e12.first].face1].removeFromHull();
  if(!_SCHtriangles[_SCHedges[e12.second].face1].inHull)
    _SCHtriangles[_SCHedges[e12.second].face2].removeFromHull();
  else _SCHtriangles[_SCHedges[e12.second].face1].removeFromHull();

  // make new triangle
  _SCHtriangles.push_back(SCHtriangle(v2,v3,v4,e34.first,e,e34.second));

  // update triangles
  updateNeighbours(e34,_SCHtriangles.size()-1);
  updateNeighbours(e,_SCHtriangles.size()-1);

  // insert triangle to heap
  _heap.push(SCHheap(findCircumSphere3(v2,v3,v4).radius,type(triangle),_SCHtriangles.size()-1));
  temp.push_back(SCHheap(findCircumSphere3(v2,v3,v4).radius,type(triangle),_SCHtriangles.size()-1));

  // insert edges to heap
  _heap.push(SCHheap(getEdgeHeap(_SCHedges[e34.first]),type(edge),e34.first));
  _heap.push(SCHheap(getEdgeHeap(_SCHedges[e]),type(edge),e));
  _heap.push(SCHheap(getEdgeHeap(_SCHedges[e34.second]),type(edge),e34.second));

  temp.push_back(SCHheap(getEdgeHeap(_SCHedges[e34.first]),type(edge),e34.first));
  temp.push_back(SCHheap(getEdgeHeap(_SCHedges[e]),type(edge),e));
  temp.push_back(SCHheap(getEdgeHeap(_SCHedges[e34.second]),type(edge),e34.second));
} // dissapearVertex

size_t SchCreator3D::findEdge(size_t v1, size_t v2)
{
  std::vector<size_t>::iterator i;
  for(i = _SCHvertexes[v1].neighbours.begin();
      i != _SCHvertexes[v1].neighbours.end(); i++)
  {
    // std::cout << _SCHedges[*i].vertex1 << ' ' << _SCHedges[*i].vertex2 << std::endl;
    if(_SCHedges[*i].vertex1 == v2 || _SCHedges[*i].vertex2 == v2)
      // std::cout << '\n' << *i << '\n' << std::endl;
      break;
  }

  return *i;
} // findEdge

void SchCreator3D::invertEdge(size_t v3, size_t v4, SCHneighbours e12, SCHneighbours e34)
{

  // remove edges and traingles from heap
  _SCHedges[_heap.top().index].removeFromHull();
  _SCHtriangles[_SCHedges[_heap.top().index].face1].removeFromHull();
  _SCHtriangles[_SCHedges[_heap.top().index].face2].removeFromHull();

  // add the new edge
  size_t index = _SCHtriangles.size();
  _SCHvertexes[v3].neighbours.push_back(_SCHedges.size());
  _SCHvertexes[v4].neighbours.push_back(_SCHedges.size());
  SCHedge e = SCHedge(v3,v4,index,index+1);
  _SCHedges.push_back(SCHedge(v3,v4,index,index+1));
  
  
  // make new edges
  updateNeighbours(e12,index);
  updateNeighbours(e34,index+1);

  // get the new triangles
  index = _SCHedges.size()-1;
  _SCHtriangles.push_back(SCHtriangle(_SCHedges[_heap.top().index].vertex1,
                                      v3, v4,index+1,index,index+2));
  _SCHtriangles.push_back(SCHtriangle(_SCHedges[_heap.top().index].vertex2,
                                      v3, v4,index+3,index,index+4));
  // insert new triangles to the heap
  std::cout << findCircumSphere3(_SCHedges[_heap.top().index].vertex1,v3, v4).radius
            << ' ' << type(triangle) << ' ' << _SCHtriangles.size()-2 << std::endl;
  std::cout << findCircumSphere3(_SCHedges[_heap.top().index].vertex2,v3, v4).radius
            << ' ' << type(triangle) << ' ' << _SCHtriangles.size()-1<< std::endl;
  _heap.push(SCHheap(findCircumSphere3(_SCHedges[_heap.top().index].vertex1,
                                      v3, v4).radius,type(triangle),_SCHtriangles.size()-2));
  _heap.push(SCHheap(findCircumSphere3(_SCHedges[_heap.top().index].vertex2,
                                      v3, v4).radius,type(triangle),_SCHtriangles.size()-1));
  temp.push_back(SCHheap(findCircumSphere3(_SCHedges[_heap.top().index].vertex1,
                                      v3, v4).radius,type(triangle),_SCHtriangles.size()-2));
  temp.push_back(SCHheap(findCircumSphere3(_SCHedges[_heap.top().index].vertex2,
                                      v3, v4).radius,type(triangle),_SCHtriangles.size()-1));

  // std::cout << e.vertex1 << ' ' << e.vertex2 << ' '
  //           << findVertex(_SCHtriangles[e.face1],e.vertex1,e.vertex2) << ' '
  //           << findVertex(_SCHtriangles[e.face2],e.vertex1,e.vertex2) 
  //           << ' ' << getEdgeHeap(e) << std::endl; 
  // add edges to the hull and heap again
  e = _SCHedges[e12.first];
  _SCHedges.push_back(e);
  // std::cout << e.vertex1 << ' ' << e.vertex2 << ' '
  //           << findVertex(_SCHtriangles[e.face1],e.vertex1,e.vertex2) << ' '
  //           << findVertex(_SCHtriangles[e.face2],e.vertex1,e.vertex2) 
  //           << ' ' << getEdgeHeap(e) << std::endl; 
  updateNeighbours(e12.first,_SCHtriangles.size()-2,_SCHedges.size()-1);
  _heap.push(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  temp.push_back(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  _SCHvertexes[e.vertex1].neighbours.push_back(_SCHedges.size()-1);
  _SCHvertexes[e.vertex2].neighbours.push_back(_SCHedges.size()-1);

  e = _SCHedges[e12.second];
  _SCHedges.push_back(e);
  // std::cout << e.vertex1 << ' ' << e.vertex2 << ' '
  //           << findVertex(_SCHtriangles[e.face1],e.vertex1,e.vertex2) << ' '
  //           << findVertex(_SCHtriangles[e.face2],e.vertex1,e.vertex2) 
  //           << ' ' << getEdgeHeap(e) << std::endl;
  updateNeighbours(e12.second,_SCHtriangles.size()-2,_SCHedges.size()-1);
  _heap.push(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  temp.push_back(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  _SCHvertexes[e.vertex1].neighbours.push_back(_SCHedges.size()-1);
  _SCHvertexes[e.vertex2].neighbours.push_back(_SCHedges.size()-1);

  e = _SCHedges[e34.first];
  _SCHedges.push_back(e);
  // std::cout << e.vertex1 << ' ' << e.vertex2 << ' '
  //           << findVertex(_SCHtriangles[e.face1],e.vertex1,e.vertex2) << ' '
  //           << findVertex(_SCHtriangles[e.face2],e.vertex1,e.vertex2) 
  //           << ' ' << getEdgeHeap(e) << std::endl;
  updateNeighbours(e34.first,_SCHtriangles.size()-1,_SCHedges.size()-1);
  _heap.push(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  temp.push_back(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  _SCHvertexes[e.vertex1].neighbours.push_back(_SCHedges.size()-1);
  _SCHvertexes[e.vertex2].neighbours.push_back(_SCHedges.size()-1);

  e = _SCHedges[e34.second];
  _SCHedges.push_back(e);
  // std::cout << e.vertex1 << ' ' << e.vertex2 << ' '
  //           << findVertex(_SCHtriangles[e.face1],e.vertex1,e.vertex2) << ' '
  //           << findVertex(_SCHtriangles[e.face2],e.vertex1,e.vertex2) 
  //           << ' ' << getEdgeHeap(e) << std::endl;
  updateNeighbours(e34.second,_SCHtriangles.size()-1,_SCHedges.size()-1);
  _heap.push(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  temp.push_back(SCHheap(getEdgeHeap(e),type(edge),_SCHedges.size()-1));
  _SCHvertexes[e.vertex1].neighbours.push_back(_SCHedges.size()-1);
  _SCHvertexes[e.vertex2].neighbours.push_back(_SCHedges.size()-1);

  // remove edges from hull
  _SCHedges[e12.first].removeFromHull();
  _SCHedges[e12.second].removeFromHull();
  _SCHedges[e34.first].removeFromHull();
  _SCHedges[e34.second].removeFromHull();


} // invertEdge

void SchCreator3D::updateNeighbours(SchCreator3D::SCHneighbours e, size_t index)
{
  updateNeighbours(e.first, index);
  updateNeighbours(e.second, index);
} // updateNeighbours

void SchCreator3D::updateNeighbours(size_t e, size_t index)
{
  if(!_SCHtriangles[_SCHedges[e].face1].inHull)
    _SCHedges[e].face1 = index;
  else
    _SCHedges[e].face2 = index;
} // updateNeighbours


void SchCreator3D::updateNeighbours(size_t oldE, size_t f, size_t e)
{
  size_t face;
  // find the face to update
  if(_SCHedges[e].face1 == f) face = _SCHedges[e].face2;
  else face = _SCHedges[e].face1;  
  // find the old edge and update it
  if(_SCHtriangles[face].e1 == oldE) _SCHtriangles[face].e1 = e;
  else if(_SCHtriangles[face].e2 == oldE) _SCHtriangles[face].e2 = e;
  else if(_SCHtriangles[face].e3 == oldE) _SCHtriangles[face].e3 = e;
  else std::cout << ":(" << std::endl;
} // updateNeighbours

bool SchCreator3D::checkSameNeighbour(size_t v1, size_t v2, type t)
{
  if(t==0)
  {
    for(auto i = _SCHvertexes[v1].neighbours.begin();
        i != _SCHvertexes[v1].neighbours.end(); i++)
    {
      if(_SCHedges[*i].vertex1 == v2 || _SCHedges[*i].vertex2 == v2)
      {
        std::cout << "i: " << *i << std::endl;
        return _SCHedges[*i].inHull;
      }
    }
  } else
  {
    std::cout << "Faces: ";
    if(_SCHedges[v1].face1 == _SCHedges[v2].face1)
    {
      std::cout << _SCHedges[v1].face1 << ' ' << _SCHedges[v2].face1 << std::endl;
      return true;
    }
    else if(_SCHedges[v1].face1 == _SCHedges[v2].face2)
    {
      std::cout << _SCHedges[v1].face1 << ' ' << _SCHedges[v2].face2 << std::endl;
      return true;
    }
    else if(_SCHedges[v1].face2 == _SCHedges[v2].face1)
    {
      std::cout << _SCHedges[v1].face2 << ' ' << _SCHedges[v2].face1 << std::endl;
      return true;
    } else if(_SCHedges[v1].face2 == _SCHedges[v2].face2)
    {
      std::cout << _SCHedges[v1].face2 << ' ' << _SCHedges[v2].face2 << std::endl;
      return true;
    }
  }

  return false;
} // checkSameNeighbour


SchCreator3D::SCHneighbours SchCreator3D::findEdge(size_t v1, size_t v2, size_t v3)
{
  SCHneighbours edges;
  std::vector<size_t> tempEdges = {_SCHtriangles[_SCHedges[_heap.top().index].face1].e1,
                              _SCHtriangles[_SCHedges[_heap.top().index].face1].e2,
                              _SCHtriangles[_SCHedges[_heap.top().index].face1].e3,
                              _SCHtriangles[_SCHedges[_heap.top().index].face2].e1,
                              _SCHtriangles[_SCHedges[_heap.top().index].face2].e2,
                              _SCHtriangles[_SCHedges[_heap.top().index].face2].e3};
  // std::cout << "FIND EDGE: " << ' ' << tempEdges[0] << ' ' << tempEdges[1] << ' ' << tempEdges[2]
  //           << ' ' << tempEdges[3] << ' ' << tempEdges[4] << ' ' << tempEdges[5] <<'\n';
  for(auto i = tempEdges.begin(); i != tempEdges.end(); i++)
  {
    // std::cout << _SCHedges[*i].vertex1 << ' ' << _SCHedges[*i].vertex2 << std::endl;
    if(_SCHedges[*i].vertex1 == v1 || _SCHedges[*i].vertex2 == v1)
      if(_SCHedges[*i].vertex1 == v2 || _SCHedges[*i].vertex2 == v2)
        // std::cout << '\n' << *i << '\n' << std::endl;
        edges.first = *i;
      else if(_SCHedges[*i].vertex1 == v3 || _SCHedges[*i].vertex2 == v3)
        // std::cout << '\n' << *i << '\n' << std::endl;
        edges.second = *i;
  }

  return edges;
} // findEdge

size_t SchCreator3D::findVertex(size_t f,size_t a,size_t b)
{
  bool inSet;
  std::set<size_t> indexes;

  // add known vertexes to the set
  indexes.insert(a);
  indexes.insert(b);

  // add first sphere vertex
  inSet = indexes.insert(_SCHtriangles[f].p1).second;
  // if p1 is not in the set, return p1
  if(inSet) return _SCHtriangles[f].p1;
  // else, try adding p2 to the set
  else inSet = indexes.insert(_SCHtriangles[f].p2).second;

  // if p2 isn't in the set, return p2
  if(inSet) return _SCHtriangles[f].p2;
  // else, return p3
  else return _SCHtriangles[f].p3;
} // insertToSet

/*
 *   Creates the output file given a filename.
 */
void SchCreator3D::writeToFile(const std::string & filename)
{
  size_t index = 0;
  // std::vector<std::pair<size_t,Eigen::Vector3d>> SCHvertexes;
  std::vector<std::pair<size_t,BigSphere>> SCHspheres;
  std::vector<std::pair<size_t,Torus>> SCHtorii;

  // for(auto i = _SCHvertexes.begin(); i != _SCHvertexes.end(); i++)
  // {
  //   if((*i).inHull)
  //     SCHvertexes.push_back(std::make_pair(index,(*i).vertex));
  //   index++;
  // }

  index = _smallSpheres.size();
  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    SCHspheres.push_back(std::make_pair(index,(*i)));
    index++;
  }

  index = _smallSpheres.size() + _bigSpheres.size();
  for(auto i = _torii.begin(); i != _torii.end(); i++)
  {
    SCHtorii.push_back(std::make_pair(index,(*i)));
    index++;
  }

  // write data into a file
  std::ofstream os;
  os.open(filename.c_str());
  os.precision(16);

  // radii used for computation on first line
  os << _r << ' ' << _R << std::endl;

  // number of small spheres
  os << ' ' << _smallSpheres.size() << std::endl;
  // writes small spheres
  index = 0;
  for(auto i = _smallSpheres.begin(); i != _smallSpheres.end(); i++)
  {
    os << (*i).radius << ' ';
    os << (*i).center[0] << ' ';
    os << (*i).center[1] << ' ';
    os << (*i).center[2] << std::endl;

    // number of neighbours
    os << _vertexNeighbours[index].size() << std::endl;
    for(auto j = _vertexNeighbours[index].begin();
        j != _vertexNeighbours[index].end(); j++)
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
  os << SCHspheres.size() << std::endl;
  for(auto i = SCHspheres.begin(); i != SCHspheres.end(); i++)
  {
    // radius and center of big sphere
    os << _R << ' ' << (*i).second.s.center[0] << ' ';
    os << (*i).second.s.center[1] << ' '
       << (*i).second.s.center[2] << std::endl;

    // coordinates of the three points touching the sphere
    os << _smallSpheres[(*i).second.p1].center[0] << ' '
       << _smallSpheres[(*i).second.p1].center[1] << ' '
       << _smallSpheres[(*i).second.p1].center[2] << ' ';
    os << _smallSpheres[(*i).second.p2].center[0] << ' '
       << _smallSpheres[(*i).second.p2].center[1] << ' '
       << _smallSpheres[(*i).second.p2].center[2] << ' ';
    os << _smallSpheres[(*i).second.p3].center[0] << ' '
       << _smallSpheres[(*i).second.p3].center[1] << ' '
       << _smallSpheres[(*i).second.p3].center[2] << std::endl;

    // normals to the sphere planes
    os << _bigSphereNormals[(*i).first].plane1.first << ' ';
    os << _bigSphereNormals[(*i).first].plane1.second[0] << ' '
       << _bigSphereNormals[(*i).first].plane1.second[1] << ' '
       << _bigSphereNormals[(*i).first].plane1.second[2] << std::endl;

    os << _bigSphereNormals[(*i).first].plane2.first << ' ';
    os << _bigSphereNormals[(*i).first].plane2.second[0] << ' '
       << _bigSphereNormals[(*i).first].plane2.second[1] << ' '
       << _bigSphereNormals[(*i).first].plane2.second[2] << std::endl;

    os << _bigSphereNormals[(*i).first].plane3.first << ' ';
    os << _bigSphereNormals[(*i).first].plane3.second[0] << ' '
       << _bigSphereNormals[(*i).first].plane3.second[1] << ' '
       << _bigSphereNormals[(*i).first].plane3.second[2] << std::endl;
  }

  // number of torii
  os << SCHtorii.size() << std::endl;
  for(auto i = SCHtorii.begin(); i != SCHtorii.end(); i++)
  {
    // remove torus boolean
    os << 1 << std::endl;

    // write torus' external radius, center and normal coordinates
    os << (*i).second.extRadius << ' ' << _R << ' ';

    os << (*i).second.center[0] << ' '
       << (*i).second.center[1] << ' '
       << (*i).second.center[2] << ' ';

    os << (*i).second.normal[0] << ' '
       << (*i).second.normal[1] << ' '
       << (*i).second.normal[2] << std::endl;

    // write torus' cones
    // cone [a,b]: a index, cosangle, axis coordinates
    os << _toriiCones[(*i).first].first.first << ' ';
    os << _toriiCones[(*i).first].first.second.cosangle << ' ';
    os << _toriiCones[(*i).first].first.second.axis[0] << ' '
       << _toriiCones[(*i).first].first.second.axis[1] << ' '
       << _toriiCones[(*i).first].first.second.axis[2] << std::endl;

    // cone [b,a]: b index, cosangle, axis coordinates
    os << _toriiCones[(*i).first].second.first << ' ';
    os << _toriiCones[(*i).first].second.second.cosangle << ' ';
    os << _toriiCones[(*i).first].second.second.axis[0] << ' '
       << _toriiCones[(*i).first].second.second.axis[1] << ' '
       << _toriiCones[(*i).first].second.second.axis[2] << std::endl;

    // write torus' planes
    // first triangle: triangle index, normal coordinates
    os << _toriiPlanes[(*i).first].first.first << ' ';
    os << _toriiPlanes[(*i).first].first.second[0] << ' '
       << _toriiPlanes[(*i).first].first.second[1] << ' '
       << _toriiPlanes[(*i).first].first.second[2] << std::endl;

    // second triangle: triangle index, normal coordinates
    os << _toriiPlanes[(*i).first].second.first << ' ';
    os << _toriiPlanes[(*i).first].second.second[0] << ' '
       << _toriiPlanes[(*i).first].second.second[1] << ' '
       << _toriiPlanes[(*i).first].second.second[2] << std::endl;
  }

  // close the file
  os.close();

  std::cout << "\nSCH Created, output file " << filename.c_str() << std::endl;
}

void SchCreator3D::printEdges()
{
  size_t index = 0;
  std::cout << "\nedges " << _SCHedges.size() << std::endl;
  std::cout << "n\tVertex\tFaces\tInHull\n";

  for(auto i = _SCHedges.begin(); i != _SCHedges.end(); i++)
  {
    std::cout << index << '\t' << (*i).vertex1 << ' '
              << (*i).vertex2 << '\t'
              << (*i).face1 << ' '
              << (*i).face2 << '\t'
              << (*i).inHull << '\n';
    index++;
  }
}

void SchCreator3D::printVertexes()
{
  std::cout << "\nvertexes " << std::endl;
  size_t index = 0;
  std::cout << "\nIndex\tInHull " << std::endl;
  for(auto i = _SCHvertexes.begin(); i != _SCHvertexes.end(); i++)
  {
    std::cout << index << '\t' << (*i).inHull << std::endl;
    std::cout << "neighbours: ";
    for(auto j = (*i).neighbours.begin(); j != (*i).neighbours.end(); j++)
      std::cout << *j << ' ';

    std::cout << std::endl;
    index++;
  }
}

void SchCreator3D::printTriangles()
{
  size_t index = 0;
  std::cout << "\nfaces " << _SCHtriangles.size() << std::endl;
  std::cout << "n\tVertexes\tEdges   \tInHull" << std::endl;

  for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  {
    std::cout << index << '\t' << (*i).p1 << ' '  << (*i).p2 << ' '  << (*i).p3 << "\t\t";
    std::cout << (*i).e1 << ' '  << (*i).e2 << ' '  << (*i).e3 << "   \t";
    std::cout << (*i).inHull << std::endl;
    index++;
  }

}

void SchCreator3D::printHeap()
{
  size_t x =  _heap.size();
  for(size_t i = 0; i < x; i++)
  {
    std::cout << i << ' ' << _heap.top().radius
              << ' ' << _heap.top().t
              << ' ' << _heap.top().index << std::endl;
    _heap.pop();
  }
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

  initialize();
  getHeap();


  printVertexes();
  printEdges();
  printTriangles();

  std::cout << _heap.size() << std::endl;

  double maxHeap = _heap.top().radius;
  _alpha = maxHeap;
  
  std::cout << "\n-------------------------------" << std::endl;

  std::cout << "Max Heap: " << maxHeap << ' '
              << _heap.top().t << ' ' << _heap.top().index << std::endl;

    std::cout << "Triangles: " << _SCHtriangles.size();
    std::cout << "\tEdges: " << _SCHedges.size() << std::endl;

  while(maxHeap > _R)
  {
    if(_heap.top().t == 0)
    {
      changeTopology();
      _heap.pop();
    }
    

    bool check = (_heap.top().t == 0) ? _SCHedges[_heap.top().index].inHull:_SCHtriangles[_heap.top().index].inHull;
    while(!check)
    {
      std::cout << "Max Heap: " << maxHeap << ' '
                << _heap.top().t << ' ' << _heap.top().index << std::endl;
      std::cout << "Check: " << check << std::endl;
      _heap.pop();
      check = (_heap.top().t == 0) ? _SCHedges[_heap.top().index].inHull:_SCHtriangles[_heap.top().index].inHull;
    }

    maxHeap = _heap.top().radius;
    _alpha = maxHeap;

    std::cout << "\n-------------------------------" << std::endl;

    std::cout << "Max Heap: " << maxHeap << ' '
              << _heap.top().t << ' ' << _heap.top().index << std::endl;
  }


  printVertexes();
  printEdges();
  printTriangles();

  updateVertexesIndex();
  // get spheres
  getSmallSpheres();
  getBigSpheres();
  getTorii();
  getVertexNeighbours();
  
  std::priority_queue<SCHheap> vHeap = std::priority_queue<SCHheap>(temp.begin(),temp.end());
  std::cout << "HEAP" << std::endl;
  for(size_t i = 0; i < vHeap.size(); i++)
  {
    std::cout << vHeap.top().radius << ' ' << vHeap.top().index << ' ' << vHeap.top().t;
    bool check = (vHeap.top().t == 0) ? _SCHedges[vHeap.top().index].inHull:_SCHtriangles[vHeap.top().index].inHull;
    std::cout << ' ' << check << std::endl;
    vHeap.pop();
  }
} // computeSCH
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