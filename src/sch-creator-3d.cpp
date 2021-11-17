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

/**
 * @brief Checks if all points are inside the given sphere by comparing
 * the distance of all vertexes to its center.
 * 
 * @param s an SchCreator3D Sphere
 * @return true if all the vertexes are inside the sphere
 * @return false if any of the vertexes is not in the sphere
 */
bool SchCreator3D::checkPointsInSphere(const Sphere & s)
{
  // set the max distance
  double distance, maxDistance = _R * _R * (1 + _epsilon);
  for(auto i = _vertexes.begin(); i != _vertexes.end(); i++)
  {
    // get the squared distance from the center of the sphere 
    // to the vertex
    distance = (s.center - (*i)).squaredNorm();

    // check if the difference between the distance and the max
    // distance is within the boundary
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

/**
 * @brief Checks to see if a vertex is inside the circumsphere of three points
 * 
 * @param a The vertex to check.
 * 
 * The circumsphere vertexes
 * @param b 
 * @param c 
 * @param d 
 * @return true if the distance between the sphere center and a is within the boundary
 * @return else, false 
 */
bool SchCreator3D::checkPointsInSphere(size_t a, size_t b, size_t c, size_t d)
{
  Sphere s = findCircumSphere3(b, c, d);
  double distance = (s.center - _vertexes[a]).squaredNorm();
  return distance < (s.radius * s.radius);
}

/**
 * @brief Finds the derivative of the distance between the center of a sphere going through 
 * three points and a vertex, with respect to the sphere's radius
 * 
 * The sphere's vertexes
 * @param v1 
 * @param v2 
 * @param v3 
 * The vertex from which we measure the distance
 * @param v4  
 * @return true if the derivative is positive. false otherwise 
 */
bool SchCreator3D::getDerivative(size_t v1, size_t v2, size_t v3, size_t v4)
{
  // get the plane through the points
  SchCreator3D::Plane p = findPlaneBase(v1, v2, v3);

  // get the circle on the plane that goes thorugh the points
  Eigen::Vector2d circle2D = findCircleThroughPoints(v1, v2, v3, p);

  // get the thistance from the center of the circel to the sphere center
  double dToSphereCenter = sqrt(_alpha * _alpha - circle2D.squaredNorm());

  // convert the circle's center coordinates to 3D
  Eigen::Vector3d circle3D;
  circle3D(0) = _vertexes[v1](0) + circle2D.dot(p.base.col(0));
  circle3D(1) = _vertexes[v1](1) + circle2D.dot(p.base.col(1));
  circle3D(2) = _vertexes[v1](2) + circle2D.dot(p.base.col(2));

  // get the derivative of the sphere center coordinates
  Eigen::Vector3d dsc = (_alpha / dToSphereCenter) * p.normal;

  return 2 * ((circle3D + dToSphereCenter * p.normal).dot(dsc) - _vertexes[v4].dot(dsc) - _alpha) > 0;
} // getDerivative

/**
 * @brief Finds the maximum body distance from the polyhedron's edges.
 * 
 * @return true if the desired circumference is equal or larger than the max distance.
 * false otherwise.
 */
bool SchCreator3D::findMaxDistance()
{
  double d;

  for(size_t i = 0; i < poly.edges_.size(); i++)
  {
    // get the length of current edge
    d = poly.edges_[i].edge.norm();
    // update largest distance
    if(maxBodyDistance < d) 
      maxBodyDistance = d;
    // update the smallest distance
    if(noise > d) 
      noise = d;
  }

  // obtain the noise from the smallest distance
  noise /= 1000;

  std::cout << "Maximum body distance: " << maxBodyDistance << std::endl;
  std::cout << "Minimum body distance: " << noise << std::endl;

  if(maxBodyDistance >= 2 * _desiredAlpha)
    return true;
  else
    return false;
} // SchCreator3D::findMaxDistance

/**
 * @brief Randomly returns a positive or negative noise
 * 
 */
double SchCreator3D::addNoise()
{
  int n = rand() % 10;
  if(n % 2 == 0)
    return -noise;
  else
    return noise;
} // addNoise

/**
 * @brief Uses poly to get the vertexes, edges and triangles of the
 * input polyhedron and creates SCH objects with them.
 * 
 */
void SchCreator3D::initialize()
{
  Vector3 vertexTemp;
  Eigen::Vector3d vertex;
  size_t id, index = 0, torusIndex = 0, e1, e2, e3;
  std::vector<size_t> tempX;
  std::multimap<size_t, size_t> orderedEdges;
  std::map<size_t, size_t> edgesByKey;
  std::multimap<size_t, size_t> orderedFaces;

  // initialize the center of the point cloud at the origin
  initialCenter = Eigen::Vector3d(0, 0, 0);

  // create the vector of vertexes
  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    // read the vertex
    vertexTemp = poly.vertexes_[i]->getCoordinates();
    // convert to Eigen Vector
    vertex = Eigen::Vector3d(vertexTemp.m_x, vertexTemp.m_y, vertexTemp.m_z);
    // add vertex to vector
    _vertexes.push_back(vertex);

    initialCenter += vertex;
  }

  // get the mean of the point cloud
  initialCenter /= int(_numberOfVertexes);

  size_t a, b, c;
  for(auto i = poly.triangles_.begin(); i != poly.triangles_.end(); i++)
  {
    // get the triangle's vertexes
    a = (*i).a;
    b = (*i).b;
    c = (*i).c;

    // order the vertexes in counter clockwise direction
    orderTriangle(a, c, b);

    // get info. on the three edges
    id = getEdgeKey(a, b);
    e1 = id;
    orderedFaces.insert(std::make_pair(id, index));
    if(edgesByKey.insert(std::make_pair(id, torusIndex)).second)
    {
      // order edges by vertex
      orderedEdges.insert(std::make_pair(a, torusIndex));
      orderedEdges.insert(std::make_pair(b, torusIndex));
      torusIndex++;
    }

    id = getEdgeKey(b, c);
    e2 = id;
    orderedFaces.insert(std::make_pair(id, index));
    if(edgesByKey.insert(std::make_pair(id, torusIndex)).second)
    {
      // order edges by vertex
      orderedEdges.insert(std::make_pair(b, torusIndex));
      orderedEdges.insert(std::make_pair(c, torusIndex));
      torusIndex++;
    }

    id = getEdgeKey(c, a);
    e3 = id;
    orderedFaces.insert(std::make_pair(id, index));
    if(edgesByKey.insert(std::make_pair(id, torusIndex)).second)
    {
      // order edges by vertex
      orderedEdges.insert(std::make_pair(c, torusIndex));
      orderedEdges.insert(std::make_pair(a, torusIndex));
      torusIndex++;
    }

    // add ccw triangle to the vector
    _SCHtriangles.push_back(SCHtriangle(a, b, c, edgesByKey[e1], edgesByKey[e2], edgesByKey[e3]));

    // increase triangle index
    index++;
  }

  // get the edges
  edgesByKey.clear();
  for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  {
    // get edge id
    id = getEdgeKey((*i).p1, (*i).p2);
    if(edgesByKey.insert(std::make_pair(id, torusIndex)).second)
    {
      // get faces
      auto it = orderedFaces.equal_range(id);
      for(auto j = it.first; j != it.second; j++) tempX.push_back((*j).second);
      // add edge to vector
      _SCHedges.push_back(SCHedge((*i).p1, (*i).p2, tempX[0], tempX[1]));
      tempX.clear();
    }

    // get edge id
    id = getEdgeKey((*i).p2, (*i).p3);
    if(edgesByKey.insert(std::make_pair(id, torusIndex)).second)
    {
      // get faces
      auto it = orderedFaces.equal_range(id);
      for(auto j = it.first; j != it.second; j++) tempX.push_back((*j).second);
      // add edge to vector
      _SCHedges.push_back(SCHedge((*i).p2, (*i).p3, tempX[0], tempX[1]));
      tempX.clear();
    }

    // get edge id
    id = getEdgeKey((*i).p3, (*i).p1);
    if(edgesByKey.insert(std::make_pair(id, torusIndex)).second)
    {
      // get faces
      auto it = orderedFaces.equal_range(id);
      for(auto j = it.first; j != it.second; j++) tempX.push_back((*j).second);
      // add edges to vector
      _SCHedges.push_back(SCHedge((*i).p3, (*i).p1, tempX[0], tempX[1]));
      tempX.clear();
    }
  }
  
  // get vertex neighbours
  for(size_t i = 0; i < _numberOfVertexes; i++)
  {
    // get all neighbours for the current vertex
    auto neighbours = orderedEdges.equal_range(i);
    // add them to a vector
    for(auto j = neighbours.first; j != neighbours.second; j++) tempX.push_back((*j).second);
    // add SCHvertex to vector
    _SCHvertexes.push_back(SCHvertex(_vertexes[i], tempX));
    tempX.clear();
  }

} // initialize

/**
 * @brief Creates a vector of spheres containing the vertexes in the hull.
 * Radius is equal to _r for all spheres.
 * 
 */
void SchCreator3D::getSmallSpheres()
{
  std::cout << "Finding small spheres... ";
  size_t index = 0;
  for(auto i = _SCHvertexes.begin(); i != _SCHvertexes.end(); i++)
  {
    if((*i).inHull)
    {
      // get small spheres
      _smallSpheres.push_back(Sphere((*i).vertex, _r));
      i->ssIndex = index;
      index++;
    }
  }
  std::cout << "Done." << std::endl;
} // getSmallSpheres

/**
 * @brief Given three vertex indexes, it finds the sphere going through them and compares
 * checs if the vertexes are ccw on the plane going through them, verifies the angle between
 * the vector going from the center to the sphere to the plane normal, and the vector from 
 * the point cloud's center to a and the plane normal are both less than 90°
 * 
 * All three parameters are vertex indexes
 * @param a 
 * @param b 
 * @param c 
 * @return * If the points are ccw, true. Else, false
 */
bool SchCreator3D::checkOrientation(size_t a, size_t b, size_t c)
{
  // find the sphere going through the points
  Sphere s = findSphereThroughPoints(a, b, c);
  // get the vector from the center of the sphere to a (Ca) 
  // and from the center of the point cloud to a (Oa)
  Eigen::Vector3d Ca = _vertexes[a] - s.center, Oa = _vertexes[a] - initialCenter;

  // get the plane going through the points
  Plane p = findPlaneBase(a, b, c);
  // find the 2D projection of the points using the base
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], ac = _vertexes[c] - _vertexes[a];
  Eigen::Vector2d a2d(0, 0), b2d = p.base * ab, c2d = p.base * ac, ab2d = b2d - a2d, ac2d = c2d - a2d;
  Eigen::Matrix2d m;
  m << ab2d.transpose(), ac2d.transpose();
  // get the determinant of m
  // this value determines whether or not the points are ccw on the plane
  double d = m.determinant();

  // if the dot products are less than 0, the angles are smaller than 90°
  return !(d > 0 && ((p.normal.dot(Ca) > 0) && (p.normal.dot(Oa) > 0)));
} // checkOrientation

/**
 * @brief Builds the big spheres in the hull from _SCHtriangles.
 * All spheres have a radius _R.
 * All the triangles must be ccw for the sphere to be computed correctly.
 * 
 */
void SchCreator3D::getBigSpheres()
{
  std::cout << "Finding big spheres...";
  Sphere s;
  BigSphere bs;
  size_t index = _smallSpheres.size();

  for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  {
    if(i->inHull)
    {
      // find the sphere going through the points
      s = findSphereThroughPoints((*i).p1, (*i).p2, (*i).p3);
      // create big sphere
      bs = BigSphere(s, _SCHvertexes[(*i).p1].ssIndex, _SCHvertexes[(*i).p2].ssIndex, _SCHvertexes[(*i).p3].ssIndex);

      // add to big sphere vector
      _bigSpheres.push_back(bs);
      i->bsIndex = index;

      // increase the current triangle index
      index++;
    }
  }

  std::cout << " Done." << std::endl;
} // SchCreator3D::getBigSpheres

/*
 *   Given 3 points, the function returns the circumsphere
 */

/**
 * @brief Finds the circumsphere of three points
 * 
 * All parameters are vertex indexes
 * @param a 
 * @param b 
 * @param c 
 * @return SchCreator3D::Sphere 3-point circumsphere
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

  return Sphere(_vertexes[a] + k1 * ab + k2 * ac, (k1 * ab + k2 * ac).norm());
} // SchCreator3D::findCircumSphere3

/**
 * @brief Finds the circumsphere of four points
 * 
 * All parameters are vertex indexes
 * @param a 
 * @param b 
 * @param c 
 * @param d 
 * @return SchCreator3D::Sphere 4-point circumsphere
 */
SchCreator3D::Sphere SchCreator3D::findCircumSphere4(size_t a, size_t b, size_t c, size_t d)
{
  // ussing crammer's rule find base matrix
  Eigen::Matrix4d T;
  T << _vertexes[a][0], _vertexes[a][1], _vertexes[a][2], 1, _vertexes[b][0], _vertexes[b][1], _vertexes[b][2], 1,
      _vertexes[c][0], _vertexes[c][1], _vertexes[c][2], 1, _vertexes[d][0], _vertexes[d][1], _vertexes[d][2], 1;

  // get squared norm values
  Eigen::Vector4d t;
  t << -_vertexes[a].squaredNorm(), -_vertexes[b].squaredNorm(), -_vertexes[c].squaredNorm(),
      -_vertexes[d].squaredNorm();

  // from base matrix get Mn
  Eigen::Matrix4d M1, M2, M3, M4;
  M1 << t, T.rightCols(3);
  M2 << T.col(0), t, T.rightCols(2);
  M3 << T.leftCols(2), t, T.col(3);
  M4 << T.leftCols(3), t;

  // get the determinants
  double Tdet = T.determinant();

  if(Tdet == 0) Tdet = 1e-10;

  double D = M1.determinant() / Tdet, E = M2.determinant() / Tdet, F = M3.determinant() / Tdet,
         G = M4.determinant() / Tdet;

  double radius = 0.5 * sqrt((D*D) + (E*E) + (F*F) - 4 * G);
  return Sphere(Eigen::Vector3d(-D / 2, -E / 2, -F / 2), radius);
} // findCircumSphere4

/**
 * @brief Finds the sphere going through 3 points with radius _desiredAlpha
 * 
 * All parameters are vertex indexes
 * @param a 
 * @param b 
 * @param c 
 * @return SchCreator3D::Sphere sphere touching 3 points
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a, size_t b, size_t c)
{
  // get plane base and normal
  SchCreator3D::Plane p = findPlaneBase(a, c, b);
  Eigen::MatrixXd base = p.base;

  // find the 2D projection of the points using the base
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], ac = _vertexes[c] - _vertexes[a];
  Eigen::Vector2d a2d(0, 0), b2d = base * ab, c2d = base * ac;

  // get the center of the circle (in 2D)
  Eigen::Vector2d circleCenter2D = findCircleThroughPoints(a2d, b2d, c2d);

  Sphere s = findSphereThroughPoints(a,p,circleCenter2D,_desiredAlpha);

  return Sphere(s.center, _R);
} // SchCreator3D::findSphereThroughPoints

/**
 * @brief Finds the sphere going through 3 points with radius _alpha
 * 
 * @param t  an SCHtriangle
 * @return SchCreator3D::Sphere 
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(const SchCreator3D::SCHtriangle & t)
{
  // get the plane through the points
  Plane p = findPlaneBase(t.p1, t.p3, t.p2);
  // get the circle on the plane that goes thorugh the points
  Eigen::Vector2d circle2D = findCircleThroughPoints(t.p1, t.p2, t.p3, p);

  return findSphereThroughPoints(t.p1, p, circle2D, _alpha);
}

/**
 * @brief Finds the sphere going through 3 points with radius _alpha
 * 
 * @param a vertex index
 * @param p plane going through the points
 * @param circleCenter2D the center of the circle on the plane touching the points
 * @return SchCreator3D::Sphere 
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a,
                                                           SchCreator3D::Plane p,
                                                           Eigen::Vector2d circleCenter2D)
{
  return findSphereThroughPoints(a, p, circleCenter2D, _alpha);
} // SchCreator3D::findSphereThroughPoints

/**
 * @brief Finds the sphere going through 3 points with radius _alpha
 * 
 * @param a vertex index
 * @param p plane going through the points
 * @param circleCenter2D the center of the circle on the plane touching the points
 * @param R the radius of the sphere
 * @return SchCreator3D::Sphere 
 */
SchCreator3D::Sphere SchCreator3D::findSphereThroughPoints(size_t a,
                                                           SchCreator3D::Plane p,
                                                           Eigen::Vector2d circleCenter2D,
                                                           double R)
{
  // convert coordinates to 3D
  Eigen::Vector3d circleCenter3D;
  circleCenter3D(0) = _vertexes[a](0) + circleCenter2D.dot(p.base.col(0));
  circleCenter3D(1) = _vertexes[a](1) + circleCenter2D.dot(p.base.col(1));
  circleCenter3D(2) = _vertexes[a](2) + circleCenter2D.dot(p.base.col(2));

  // get sphere and circle radius
  double sphereRadius = R;
  // get the distance from the center of the sphere
  double distanceFromCenter;
  if(circleCenter2D.norm() > sphereRadius)
  {
    // we do this to avoid any NaN values in the sphere center
    distanceFromCenter = 0;
  }
  else
  {
    distanceFromCenter = sqrt(sphereRadius * sphereRadius - circleCenter2D.squaredNorm());
  }

  // find the center of the sphere
  Eigen::Vector3d sphereCenter = circleCenter3D + (distanceFromCenter + _epsilon) * p.normal;

  return Sphere(sphereCenter, R);
} // SchCreator3D::findSphereThroughPoints

/**
 * @brief Finds the plane going through tree points
 * 
 * All parameters are indexes of the vertexes
 * @param a 
 * @param b 
 * @param c 
 * @return SchCreator3D::Plane 
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

/**
 * @brief Finds the circumcircle of three points in a 2D space
 * 
 * 2D vectors representing the reflections of the 3D vertexes on the plane going 
 * through them
 * @param a 
 * @param b 
 * @param c 
 * @return Eigen::Vector2d 
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
} // SchCreator3D::findCircleThroughPoints


/**
 * @brief Finds the circumcircle of three points on the plane
 * 
 * @param a vertex index
 * @param b vertex index
 * @param c vertex index
 * @param p plane going through the points
 * @return Eigen::Vector2d 
 */
Eigen::Vector2d SchCreator3D::findCircleThroughPoints(size_t a, size_t b, size_t c, SchCreator3D::Plane p)
{
  // find the 2D projection of the points using the base
  Eigen::Vector3d ab = _vertexes[b] - _vertexes[a], ac = _vertexes[c] - _vertexes[a];
  Eigen::Vector2d a2d(0, 0), b2d = p.base * ab, c2d = p.base * ac;

  return findCircleThroughPoints(a2d, b2d, c2d);
} // SchCreator3D::findCircleThroughPoints

/**
 * @brief Builds the torii from _SCHedges. This includes building the cones and the plane normals
 * corresponding to each torus.
 * 
 */
void SchCreator3D::getTorii()
{
  std::cout << "Finding torii...";
  size_t torusIndex = _smallSpheres.size() + _bigSpheres.size(), index = _smallSpheres.size();
  std::multimap<size_t, SCHplane> orderedPlanes;
  Eigen::Vector3d planeNormal;
  Face face;
  Torus torus;
  SCHplane plane;
  size_t diff = _smallSpheres.size();

  // go through all edges
  for(auto i = _SCHedges.begin(); i != _SCHedges.end(); i++)
  {
    // check if the edge is part of the hull
    if(i->inHull)
    {
      // get the torus
      torus = getTorus(_bigSpheres[_SCHtriangles[i->face1].bsIndex - diff].s, i->vertex1, i->vertex2);
      // add torus to _torii
      _torii.push_back(torus);
      // get torus cones                torusIndex, pair(SCHcone,SCHcone)
      _toriiCones.insert(std::make_pair(
          torusIndex, getCones(_SCHvertexes[i->vertex1].ssIndex, _SCHvertexes[i->vertex2].ssIndex, torus.normal)));
      // set edge.torusIndex to torusIndex
      i->torusIndex = torusIndex;

      torusIndex++;
    }
  }

  index = 0;
  for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  {
    if(i->inHull)
    {
      // get the torus index
      torusIndex = _SCHedges[i->e1].torusIndex;

      // get the first plane normal
      planeNormal = getPlaneNormal(_SCHvertexes[(*i).p1].ssIndex, _SCHvertexes[(*i).p2].ssIndex,
                                   _bigSpheres[i->bsIndex - diff].s.center);
      face.plane1 = std::make_pair(torusIndex, planeNormal);

      torusIndex = _SCHedges[i->e2].torusIndex;

      // get the second plane normal
      planeNormal = getPlaneNormal(_SCHvertexes[(*i).p2].ssIndex, _SCHvertexes[(*i).p3].ssIndex,
                                   _bigSpheres[i->bsIndex - diff].s.center);
      face.plane2 = std::make_pair(torusIndex, planeNormal);

      torusIndex = _SCHedges[i->e3].torusIndex;

      // get the third plane normal
      planeNormal = getPlaneNormal(_SCHvertexes[(*i).p3].ssIndex, _SCHvertexes[(*i).p1].ssIndex,
                                   _bigSpheres[i->bsIndex - diff].s.center);
      face.plane3 = std::make_pair(torusIndex, planeNormal);

      // add face to sphere normals
      _bigSphereNormals.insert(std::make_pair(i->bsIndex, face));
    }
    index++;
  }

  index = _smallSpheres.size();
  for(auto i = _bigSphereNormals.begin(); i != _bigSphereNormals.end(); i++)
  {
    // add all big sphere normals to the map:
    // orderedPlanes<torus index,<triangle index,plane normal>>
    orderedPlanes.insert(std::make_pair((*i).second.plane1.first, std::make_pair(index, (*i).second.plane1.second)));
    orderedPlanes.insert(std::make_pair((*i).second.plane2.first, std::make_pair(index, (*i).second.plane2.second)));
    orderedPlanes.insert(std::make_pair((*i).second.plane3.first, std::make_pair(index, (*i).second.plane3.second)));
    // all keys have exactly two SCHplane elements
    index++;
  }

  index = _smallSpheres.size() + _bigSpheres.size();
  for(auto i = _toriiCones.begin(); i != _toriiCones.end(); i++)
  {
    // get the iterator pointing to planes of torus i->first
    auto it = orderedPlanes.equal_range(i->first);

    // get the plane normal of the first plane
    plane = (*it.first).second;
    // advance the iterator to get second plane
    it.first++;
    // insert to _toriiPlanes
    // torus index
    _toriiPlanes.insert(std::make_pair(i->first,
                                       // first SCHplane, second SCHplane
                                       std::make_pair(plane, (*it.first).second)));
  }

  std::cout << " Done." << std::endl;
} // getTorii

/**
 * @brief Builds the cones corresponding to a torus
 * 
 * @param a small sphere index
 * @param b small sphere index
 * @param n torus normal
 * @return std::pair<SchCreator3D::SCHcone, SchCreator3D::SCHcone> 
 */
std::pair<SchCreator3D::SCHcone, SchCreator3D::SCHcone> SchCreator3D::getCones(size_t a,
                                                                               size_t b,
                                                                               const Eigen::Vector3d & n)
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

/**
 * @brief Computes the cosine of the cone
 * 
 * @param a small sphere index
 * @param b small sphere index
 * @return double 
 */
double SchCreator3D::getCosine(size_t a, size_t b)
{
  return (_smallSpheres[a].center - _smallSpheres[b].center).norm() / (2 * _desiredAlpha);
}

/**
 * @brief Computes the plane normals
 * 
 * @param a small sphere index
 * @param b small sphere index
 * @param c big sphere center
 * @return Eigen::Vector3d 
 */
Eigen::Vector3d SchCreator3D::getPlaneNormal(size_t a, size_t b, Eigen::Vector3d c)
{
  // find the normal to the plane
  Eigen::Vector3d ca = _smallSpheres[a].center - c, 
                  cb = _smallSpheres[b].center - c;
  return ca.cross(cb).normalized();
}

/**
 * @brief Computes the plane normals
 * 
 * @param a vertex index
 * @param b vertex index
 * @param f SCHtriangle (face)
 * @return Eigen::Vector3d 
 */
Eigen::Vector3d SchCreator3D::getPlaneNormal(size_t a, size_t b, size_t f)
{
  // find the sphere center
  Eigen::Vector3d c = findSphereThroughPoints(_SCHtriangles[f]).center;
  Eigen::Vector3d ca = _vertexes[a] - c, 
                  cb = _vertexes[b] - c;
  return ca.cross(cb).normalized();
}

/**
 * @brief Computes a torus
 * 
 * @param s SchCreator3D sphere
 * @param a vertex index
 * @param b vertex index
 * @return SchCreator3D::Torus 
 */
SchCreator3D::Torus SchCreator3D::getTorus(const Sphere & s, size_t a, size_t b)
{
  Eigen::Vector3d center = (_vertexes[a] + _vertexes[b]) / 2;
  return Torus(center, (_vertexes[a] - _vertexes[b]).normalized(), (center - s.center).norm());
}

/**
 * @brief Finds a unique key/id given two numvers
 * 
 * @param a vertex index
 * @param b vertex index
 * @return size_t 
 */
size_t SchCreator3D::getEdgeKey(size_t a, size_t b)
{
  return ((a < b) ? (a * _numberOfVertexes + b) : (b * _numberOfVertexes + a));
}

/**
 * @brief Finds a unique key/id given two numvers
 * 
 * @param edge edge index
 * @return size_t 
 */
size_t SchCreator3D::getEdgeKey(size_t edge)
{
  return getEdgeKey(_SCHedges[edge].vertex1,_SCHedges[edge].vertex2);
}

/**
 * @brief Finds a key given two numbers
 * 
 * @param a 
 * @param b 
 * @return size_t 
 */
size_t SchCreator3D::getKey(size_t a, size_t b)
{
  return (a * _numberOfVertexes + b);
}

/**
 * @brief Finds and orders all active neighbours for the vertexes
 * 
 */
void SchCreator3D::getVertexNeighbours()
{
  Cone cone;
  size_t index = 0;
  std::multimap<size_t, size_t> orderedEdges;
  _vertexNeighbours.resize(_smallSpheres.size());
  std::vector<std::vector<size_t>> edgeIDs;
  edgeIDs.resize(_smallSpheres.size());

  // go through all vertexes in hull
  for(auto i = _SCHvertexes.begin(); i != _SCHvertexes.end(); i++)
  {
    if((*i).inHull) 
    {
      // go through all of its neighbours 
      for(auto j = (*i).neighbours.begin(); j != (*i).neighbours.end(); j++)
      {
        if(_SCHedges[*j].inHull)
        {
          // if the neighbour is in the hull, add it to edgeIDs
          edgeIDs[index].push_back(_SCHedges[*j].torusIndex);
        }
      }
      index++;
    }
  }

  index = 0;
  size_t difference = _smallSpheres.size() + _bigSpheres.size();
  // go through all small spheres
  for(auto i = _smallSpheres.begin(); i != _smallSpheres.end(); i++)
  {
    // get the neighbours from edgeIDs
    for(auto j = edgeIDs[index].begin(); j != edgeIDs[index].end(); j++)
    {
      // check if the index corresponds to the first or second SCHCone
      if(index == _toriiCones[*j].first.first)
      {
        cone = _toriiCones[*j].first.second;
      }
      else
      {
        cone = _toriiCones[*j].second.second;
      }
      // add the corresponding cone to the vertexe's neighbours
      _vertexNeighbours[index].push_back(std::make_pair(*j, cone));
    }

    index++;
  }
} // getVertexNeighbours();

/**
 * @brief Computes the initial heap for all edges and triangles.
 * 
 */
void SchCreator3D::getHeap()
{
  Sphere s;
  std::set<size_t> indexes;
  std::vector<SCHheap> tempHeap;

  size_t index = 0;
  for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
  {
    // get circumsphere
    s = findCircumSphere3((*i).p1, (*i).p2, (*i).p3);
    tempHeap.push_back(SCHheap(s.radius, type(triangle), index));

    index++;
  }

  index = 0;
  for(auto i = _SCHedges.begin(); i != _SCHedges.end(); i++)
  {
    // insert circumradius to the heap
    tempHeap.push_back(SCHheap(getEdgeHeap(*i), type(edge), index));
    index++;
  }

  // build the heap
  _heap = std::priority_queue<SCHheap>(tempHeap.begin(), tempHeap.end());
} // getHeap


/**
 * @brief Finds an edge's heap
 * 
 * @param e SCH edge
 * @return double 
 */
double SchCreator3D::getEdgeHeap(const SchCreator3D::SCHedge & e)
{
  Sphere s;
  size_t vertex3, vertex4;

  // get the four vertexes
  vertex3 = findVertex(_SCHtriangles[e.face1], e.vertex1, e.vertex2);
  vertex4 = findVertex(_SCHtriangles[e.face2], e.vertex1, e.vertex2);
  // compute the circum sphere
  s = findCircumSphere4(e.vertex1, e.vertex2, vertex3, vertex4);
#ifdef DISPLAY_INFO
  std::cout << e.vertex1 << ' ' << e.vertex2 << ' ' << vertex3 << ' ' << vertex4 << ' ' << s.radius << std::endl;
#endif
  return s.radius;
} // getEdgeHeap

/**
 * @brief Finds the missing vertex
 * 
 * @param t an SCHtriangle
 * @param a triangle's known vertex index
 * @param b triangle's known vertex index
 * @return size_t 
 */
size_t SchCreator3D::findVertex(const SCHtriangle & t, size_t a, size_t b)
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
  else
    inSet = indexes.insert(t.p2).second;
  // if p2 isn't in the set, return p2
  if(inSet) return t.p2;
  // else, return p3
  else
    return t.p3;
} // findVertex

/**
 * @brief Simple swap function for size_t types
 * 
 * @param a 
 * @param b 
 */
void SchCreator3D::swap(size_t & a, size_t & b)
{
  size_t temp = a;
  a = b;
  b = temp;
} // swap

/**
 * @brief Simple swap function for SCHneighbours types
 * 
 * @param a 
 * @param b 
 */
void SchCreator3D::swap(SchCreator3D::SCHneighbours & a, SchCreator3D::SCHneighbours & b)
{
  SCHneighbours temp = a;
  a = b;
  b = temp;
} // swap

/**
 * @brief Given the indexes of three points, it checks if the orientation of the triangle is 
 * counter clockwise and swaps the last two elements if it is not
 * 
 * All three parameters are vertex indexes
 * @param a 
 * @param b 
 * @param c 
 */
void SchCreator3D::orderTriangle(size_t a, size_t & b, size_t & c)
{
  if(!checkOrientation(a, b, c))
  {
    swap(b, c);
  } 
} // orderTriangle

/**
 * @brief Check for limit case. Verifies the four last points are all inside the 
 * heap's circumsphere. Sets _limitCase to true and changes _desiredAlpha and _R's 
 * value if they aren't.
 * 
 * 
 * @param a the test vertex
 * @param b sphere vertex
 * @param c sphere vertex
 * @param d sphere vertex
 * @param r max heap radius
 * @return true if we've reached the limit case. false otherwise
 */
bool SchCreator3D::checkLimitCase(size_t a, size_t b, size_t c, size_t d, double r)
{
  if(!checkPointsInSphere(a, b, c, d))
  {
    std::cout << "\nNo more vertexes can be removed. Stoping the algorithm at radius " 
              << r << "...\n" << std::endl;
    _limitCase = true;
    _desiredAlpha = r;
    return true;
  }

  return false;
}

/**
 * @brief Executes a change in topology, it either inverts the edge, disappears the
 * vertex under the sphere or does nothing.
 * 
 * @param heap max heap
 */
void SchCreator3D::changeTopology(SCHheap heap)
{
  size_t f1, f2;
  size_t a, b, c, d, i, c_copy, d_copy;
  SCHneighbours fg, hj, fg_copy, hj_copy;

  // get the edge vertexes
  a = _SCHedges[heap.index].vertex1;
  b = _SCHedges[heap.index].vertex2;
  // get the neighbour faces
  f1 = _SCHedges[heap.index].face1;
  f2 = _SCHedges[heap.index].face2;
  // get the neighbour vertexes from the faces
  c = findVertex(f1, a, b);
  d = findVertex(f2, a, b);
  // find the edges conecting a and b to c and d
  fg = findEdge(f1, f2, a, c, d);
  hj = findEdge(f1, f2, b, c, d);
  // security copies of the relevant data
  c_copy = c;
  d_copy = d;
  fg_copy = fg;
  hj_copy = hj;

  // check for thickness of edge in the heap
#ifdef DISPLAY_INFO
  std::cout << "Check thickness of edge " << heap.index << ". Is torus thin? ";
#endif
  if(torusThicknessCheck(_SCHedges[heap.index].vertex1, _SCHedges[heap.index].vertex2, _SCHedges[heap.index].face1,
                         _SCHedges[heap.index].face2))
  {
#ifdef DISPLAY_INFO
    std::cout << "1.\nRELEVANT DATA" << std::endl;
    std::cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d << std::endl;
    std::cout << "f: " << fg.first << " g: " << fg.second 
              << " h: " << hj.first << " j: " << hj.second << std::endl;

    std::cout << "Check if cd have the same neighbour: ";
#endif
    // check if c and d are connected by an edge
    if(!checkSameNeighbour(c, d, type(edge)))
    {
      // if they are not connected, invert the edge
      invertEdge(heap, c, d, fg, hj);

#ifdef DISPLAY_INFO
      std::cout << "0. Invert Edge\n" << std::endl;
#endif
    }
    else
    {
      // find the edge connecting c and d
      i = findEdge(c, d);
#ifdef DISPLAY_INFO
      std::cout << "1. Edge i: " << i << std::endl;
#endif

      // check for limit case
      if(_activeVertexes == 4)
      {
        if(checkLimitCase(a, b, c, d, heap.radius))
        {
          if(checkLimitCase(b, a, c, d, heap.radius))
          {
            // no more vertexes can be removed
            return;
          }
          else
          {
#ifdef DISPLAY_INFO
            std::cout << "disappear vertex " << b << ". LIMIT CASE. \n" << std::endl;
#endif
            disappearVertex(heap, b, c, d, a, hj, fg, i);
            _limitCase = true;
            _desiredAlpha = heap.radius;
          }
        }
        else
        {
#ifdef DISPLAY_INFO
          std::cout << "disappear vertex " << a << ". LIMIT CASE. \n" << std::endl;
#endif
          disappearVertex(heap, a, d, c, b, fg, hj, i);
          _limitCase = true;
          _desiredAlpha = heap.radius;
        }
      }
      else
      {
#ifdef DISPLAY_INFO
        std::cout << "Check which vertex disappears" << std::endl;
        std::cout << "Do f and i have the same face as neighbour?" << std::endl;
#endif
        // to check which vertex disappears first 
        // check if edges fg.second and i share a face as neighbour
        if(checkSameNeighbour(fg.second, i, type(triangle)))
        {
          // a could disappear
#ifdef DISPLAY_INFO
          std::cout << " Active neighbours of vertex" << a 
                    << ": " << findActiveNeighbours(a) << std::endl;
          std::cout << "Vertex a disappears." << std::endl;
#endif
        }
        else
        {
          // b disappears
#ifdef DISPLAY_INFO
          std::cout << " Active neighbours of vertex" << b 
                    << ": " << findActiveNeighbours(b) << std::endl;
          std::cout << "Vertex b disappears. Updating values..." << std::endl;
#endif
          // update/swap elements
          swap(a, b);
          swap(c, d);
          swap(fg, hj);
#ifdef DISPLAY_INFO
          // print the updated data
          std::cout << "RELEVANT DATA" << std::endl;
          std::cout << "a: " << a << " b: " << b 
                    << " c: " << c << " d: " << d << std::endl;
          std::cout << "f: " << fg.first << " g: " << fg.second 
                    << " h: " << hj.first << " j: " << hj.second << std::endl;
#endif
        }
        // check thicknes of the edges connected to a
#ifdef DISPLAY_INFO
        std::cout << "Check if vertex should disappear\n" 
                  << "Check thickness of torii " << fg.first 
                  << " and " << fg.second << std::endl;
#endif
        if(checkTorii(fg.first, fg.second))
        {
          // disappear the vertex
#ifdef DISPLAY_INFO
          std::cout << "disappear vertex " << a << '\n' << std::endl;
#endif
          disappearVertex(heap, a, d, c, b, fg, hj, i);
        }
        else
        {
          // invert the edge
#ifdef DISPLAY_INFO
          std::cout << "Invert Edge. Torus Check. \n" << std::endl;
#endif
          // here we must use the copy of the original c, d, fg and hj
          // just in case that the original values were swapped,
          // this ensures the new triangles are counter-clockwise.
          invertEdge(heap, c_copy, d_copy, fg_copy, hj_copy);
        }
      }
    }
  }
  else
  {
    // do nothing
    // continue to the next element in the heap
#ifdef DISPLAY_INFO
    std::cout << " 0.\n";
#endif
  }

} // changeTopology

size_t SchCreator3D::findCommonFace(size_t e1, size_t e2)
{
  std::set<size_t> faces;
  faces.insert(_SCHedges[e1].face1);
  faces.insert(_SCHedges[e1].face2);

  if(faces.insert(_SCHedges[e2].face1).second)
    return _SCHedges[e2].face1;
  else if(faces.insert(_SCHedges[e2].face2).second)
    return _SCHedges[e2].face2;
  else
    return _SCHtriangles.size();
}

size_t SchCreator3D::findActiveNeighbours(size_t v)
{
  size_t activeNeighbours = 0;
  for(auto i = _SCHvertexes[v].neighbours.begin(); i != _SCHvertexes[v].neighbours.end(); i++)
  {
    if(_SCHedges[*i].inHull) activeNeighbours++;
  }
  return activeNeighbours;
}

/**
 * @brief Change of topology. The vertex disappears under the edge.
 * 
 * @param heap max heap
 */
void SchCreator3D::disappearUnderEdge(SchCreator3D::SCHheap heap)
{
#ifdef DISPLAY_INFO
  std::cout << "\ndisappear under edge. ";
#endif

  // a is the vertex to disappear. b and c are its immediate neighbours
  size_t a, b, c, d, e;
  size_t e1 = _SCHtriangles[heap.index].e1, e2 = _SCHtriangles[heap.index].e2, e3 = _SCHtriangles[heap.index].e3, e4;
  // triangle in heap
  size_t f1 = heap.index, f3, f4;
  // get the missing neighbour triangles
  size_t f2 = _SCHedges[e1].face1 == f1 ? _SCHedges[e1].face2 : _SCHedges[e1].face1;
  size_t temp = _SCHedges[e2].face1 == f1 ? _SCHedges[e2].face2 : _SCHedges[e2].face1;
  // if f2 and f3 are the same, the vertex in common is the vertex to disappear
  if(f2 == temp)
  {
    a = findVertex(e1, e2);
    c = _SCHedges[e1].vertex1 == a ? _SCHedges[e1].vertex2 : _SCHedges[e1].vertex1;
    b = _SCHedges[e2].vertex1 == a ? _SCHedges[e2].vertex2 : _SCHedges[e2].vertex1;
    f3 = _SCHedges[e3].face1 == f1 ? _SCHedges[e3].face2 : _SCHedges[e3].face1;
    e4 = findEdge(f2, e1, e2);
  }
  else
  {
    f3 = temp;
    // get the face for the last edge
    temp = _SCHedges[e3].face1 == f1 ? _SCHedges[e3].face2 : _SCHedges[e3].face1;
    // if f2 and temp are the same, the vertex in common is the vertex to disappear
    if(f2 == temp)
    {
      a = findVertex(e1, e3);
      b = _SCHedges[e1].vertex1 == a ? _SCHedges[e1].vertex2 : _SCHedges[e1].vertex1;
      c = _SCHedges[e3].vertex1 == a ? _SCHedges[e3].vertex2 : _SCHedges[e3].vertex1;
      f2 = temp;
      e4 = findEdge(f2, e1, e3);
    }
    else
    {
      // else, the vertex in common for f3 and temp is the vertex to disappear
      a = findVertex(e2, e3);
      b = _SCHedges[e3].vertex1 == a ? _SCHedges[e3].vertex2 : _SCHedges[e3].vertex1;
      c = _SCHedges[e2].vertex1 == a ? _SCHedges[e2].vertex2 : _SCHedges[e2].vertex1;
      f2 = f3;
      f3 = temp;
      e4 = findEdge(f2, e2, e3);
    }
  }

  // find the missing vertex from f3
  d = findVertex(f3, b, c);

  // find the face conecting bce/f4 and e
  f4 = _SCHedges[e4].face1 == f2 ? _SCHedges[e4].face2 : _SCHedges[e4].face1;
  e = findVertex(f4, b, c);

#ifdef DISPLAY_INFO
  std::cout << "RELEVANT DATA\n" << a << ' ' << b << ' ' << c << ' ' << d << ' ' << e << std::endl;
  std::cout << f1 << ' ' << f2 << ' ' << f3 << ' ' << f4 << std::endl;
  std::cout << e1 << ' ' << e2 << ' ' << e3 << ' ' << e4 << std::endl;
#endif

  // remove vertex from hull
  _SCHvertexes[a].removeFromHull();

  // remove triangles from hull
  _SCHtriangles[f1].removeFromHull();
  _SCHtriangles[f2].removeFromHull();
  _SCHtriangles[f3].removeFromHull();
  _SCHtriangles[f4].removeFromHull();

  // get edges for face f3 and f4
  SCHneighbours ef3 = findEdges(f3, d, b, e2);
  SCHneighbours ef4 = findEdges(f4, e, b, e4);

  size_t triangleIndex = _SCHtriangles.size(), edgeIndex = _SCHedges.size();

#ifdef DISPLAY_INFO
  std::cout << ef3.first << ' ' << ef3.second << ' ' << ef4.first << ' ' << ef4.second << std::endl;
  std::cout << "disappear vertex " << a << " under edge " << edgeIndex + 4 << std::endl;

  std::cout << "\nMake new triangles:\n";
#endif

  // make new triangles
  _SCHtriangles.push_back(SCHtriangle(b, d, c, edgeIndex, edgeIndex + 1, edgeIndex + 4));
  _SCHtriangles.push_back(SCHtriangle(b, c, e, edgeIndex + 4, edgeIndex + 3, edgeIndex + 2));

#ifdef DISPLAY_INFO
  std::cout << triangleIndex << ' ' << b << ' ' << d << ' ' << c << ' ' << edgeIndex << ' ' << edgeIndex + 1
            << edgeIndex + 4 << std::endl;
  std::cout << triangleIndex + 1 << ' ' << b << ' ' << c << ' ' << e << ' ' << edgeIndex + 4 << ' ' << edgeIndex + 3
            << edgeIndex + 2 << std::endl;

  std::cout << "\nMake new edges:\n";
#endif

  // update old faces
  updateNeighbours(ef3, triangleIndex);
  updateNeighbours(ef4, triangleIndex + 1);

  // make new edges
  makeNewEdge(_SCHedges[ef3.first], ef3.first, edgeIndex, triangleIndex, heap.radius);
  makeNewEdge(_SCHedges[ef3.second], ef3.second, edgeIndex + 1, triangleIndex, heap.radius);
  makeNewEdge(_SCHedges[ef4.first], ef4.first, edgeIndex + 2, triangleIndex + 1, heap.radius);
  makeNewEdge(_SCHedges[ef4.second], ef4.second, edgeIndex + 3, triangleIndex + 1, heap.radius);

  // remove edges
  _SCHedges[e1].removeFromHull();
  _SCHedges[e2].removeFromHull();
  _SCHedges[e3].removeFromHull();
  _SCHedges[e4].removeFromHull();

  _SCHedges[ef3.first].removeFromHull();
  _SCHedges[ef3.second].removeFromHull();
  _SCHedges[ef4.first].removeFromHull();
  _SCHedges[ef4.second].removeFromHull();

  // make the edge under which the vertex disappears
  SCHedge newEdge = SCHedge(b, c, triangleIndex + 1, triangleIndex);
  _SCHedges.push_back(newEdge);
#ifdef DISPLAY_INFO
  std::cout << edgeIndex + 4 << ' ';
#endif
  double radius = getEdgeHeap(newEdge);
  if(radius < heap.radius)
    _heap.push(SCHheap(radius, type(edge), edgeIndex + 4));
  else
  {
    std::cout << "This radius is larger than current max heap" << std::endl;
  }

  addToVertexNeighbours(edgeIndex + 4);

  _activeVertexes--;
} // disappearUnderEdge

/**
 * @brief Find a face's edges that are not e. 
 * 
 * @param f face index
 * @param a vertex index
 * @param b vertex index
 * @param e unwanted edge index
 * @return SchCreator3D::SCHneighbours 
 */
SchCreator3D::SCHneighbours SchCreator3D::findEdges(size_t f, size_t a, size_t b, size_t e)
{
  SCHneighbours edges;
  size_t key = getEdgeKey(a, b);
  size_t e1 = _SCHtriangles[f].e1 == e ? _SCHtriangles[f].e3 : _SCHtriangles[f].e1,
         e2 = _SCHtriangles[f].e2 == e ? _SCHtriangles[f].e3 : _SCHtriangles[f].e2;
  if(getEdgeKey(e1) == key)
  {
    edges.first = e1;
    edges.second = e2;
  }
  else
  {
    edges.first = e2;
    edges.second = e1;
  }

  return edges;
} // findEdges

/**
 * @brief Find a triangle's missing edge.
 * 
 * @param f face index
 * @param e1 edge index
 * @param e2 edge index
 * @return size_t edge index
 */
size_t SchCreator3D::findEdge(size_t f, size_t e1, size_t e2)
{
  bool inSet;
  std::set<size_t> indexes;

  // insert the known edges
  indexes.insert(e1);
  indexes.insert(e2);

  // insert the first edge
  inSet = indexes.insert(_SCHtriangles[f].e1).second;
  // if the edge is inserted to the set, return its index
  if(inSet) return _SCHtriangles[f].e1;
  // else insert the next edge
  else
    inSet = indexes.insert(_SCHtriangles[f].e2).second;

  // if the edge is inserted to the set, return its index
  if(inSet) return _SCHtriangles[f].e2;
  // else return the last edge
  else
    return _SCHtriangles[f].e3;
} // findEdge

/**
 * @brief Removes all of a vertex's neighbours from the hull
 * 
 * @param v vertex index
 */
void SchCreator3D::removeNeighboursFromHull(size_t v)
{
  for(auto i = _SCHvertexes[v].neighbours.begin(); i != _SCHvertexes[v].neighbours.end(); i++)
    _SCHedges[*i].removeFromHull();

} // removeNeighboursFromHull

/**
 * @brief Change of topology. Disppears a vertex under a sphere/face.
 * 
 * @param heap max heap
 * @param v1 disappearing vertex
 * @param v2 sphere vertex
 * @param v3 sphere vertex
 * @param v4 sphere vertex
 * @param e12 pair of edges conected to v1
 * @param e34 pair of edges conected to v4
 * @param e edge connecting v1 and v2
 */
void SchCreator3D::disappearVertex(SCHheap heap,
                                   size_t v1,
                                   size_t v2,
                                   size_t v3,
                                   size_t v4,
                                   SCHneighbours e12,
                                   SCHneighbours e34,
                                   size_t e)
{
  // remove elements from hull
  _SCHvertexes[v1].removeFromHull();
  removeNeighboursFromHull(v1);
  _SCHtriangles[_SCHedges[heap.index].face1].removeFromHull();
  _SCHtriangles[_SCHedges[heap.index].face2].removeFromHull();

  // remove common face
  if(!_SCHtriangles[_SCHedges[e12.first].face1].inHull)
    _SCHtriangles[_SCHedges[e12.first].face2].removeFromHull();
  else
    _SCHtriangles[_SCHedges[e12.first].face1].removeFromHull();

  if(!_SCHtriangles[_SCHedges[e12.second].face1].inHull)
    _SCHtriangles[_SCHedges[e12.second].face2].removeFromHull();
  else
    _SCHtriangles[_SCHedges[e12.second].face1].removeFromHull();

  // make new triangle
#ifdef DISPLAY_INFO
  std::cout << "New triangle: " << std::endl;
#endif
  size_t newT = _SCHtriangles.size(), index = _SCHedges.size();

  _SCHtriangles.push_back(SCHtriangle(v4, v3, v2, index + 2, index + 1, index));

  // update triangles
  updateNeighbours(e34, newT);
  updateNeighbours(e, newT);

  // insert triangle to heap
#ifdef DISPLAY_INFO
  std::cout << v4 << ' ' << v3 << ' ' << v2 << ' ' << findCircumSphere3(v2, v3, v4).radius << std::endl;
#endif
  double ccradius = findCircumSphere3(v2, v3, v4).radius;
  if(ccradius > heap.radius)
  { 
    std::cout << "New heap is larger than current max heap" << std::endl;
  }
  else
  {
    _heap.push(SCHheap(ccradius, type(triangle), newT));
  }

  // insert edge to heap, update edge
#ifdef DISPLAY_INFO
  std::cout << "New edges: " << std::endl;
#endif
  makeNewEdge(_SCHedges[e34.first], e34.first, index, newT, heap.radius);
  makeNewEdge(_SCHedges[e], e, index + 1, newT, heap.radius);
  makeNewEdge(_SCHedges[e34.second], e34.second, index + 2, newT, heap.radius);

  // remove old edges from hull
  _SCHedges[e34.first].removeFromHull();
  _SCHedges[e].removeFromHull();
  _SCHedges[e34.second].removeFromHull();

  _activeVertexes--;
} // disapearVertex

/**
 * @brief Adds the edge to its vertexes neighbours
 * 
 * @param e an edge index
 */
void SchCreator3D::addToVertexNeighbours(size_t e)
{
  _SCHvertexes[_SCHedges[e].vertex1].neighbours.push_back(e);
  _SCHvertexes[_SCHedges[e].vertex2].neighbours.push_back(e);
} // addToVertexNeighbours

/**
 * @brief Finds the edge connecting two vertexes by looking through their neighbours
 * 
 * @param v1 vertex index
 * @param v2 vertex index
 * @return size_t, edge index
 */
size_t SchCreator3D::findEdge(size_t v1, size_t v2)
{
  std::vector<size_t>::iterator i;
  for(i = _SCHvertexes[v1].neighbours.begin(); i != _SCHvertexes[v1].neighbours.end(); i++)
  {
    if(_SCHedges[*i].vertex1 == v2 || _SCHedges[*i].vertex2 == v2)
      if(_SCHedges[*i].inHull) break;
  }

  return *i;
} // findEdge

/**
 * @brief Change in topology. Inverts the edge in the heap and creates a new edge from v3 to v4.
 * 
 * @param heap max heap
 * @param v3 vertex index
 * @param v4 vertex index
 * @param e12 pair of edges from the heap's vertexes to v3
 * @param e34 pair of edges from the heap's vertexes to v4
 */
void SchCreator3D::invertEdge(SCHheap heap, size_t v3, size_t v4, SCHneighbours e12, SCHneighbours e34)
{
  // remove edges and traingles from heap
  _SCHedges[heap.index].removeFromHull();
  _SCHtriangles[_SCHedges[heap.index].face1].removeFromHull();
  _SCHtriangles[_SCHedges[heap.index].face2].removeFromHull();

  // add the new edge
  size_t index = _SCHtriangles.size();
  _SCHedges.push_back(SCHedge(v3, v4, index, index + 1));
  addToVertexNeighbours(_SCHedges.size() - 1);

  // make new edges
  updateNeighbours(e12, index + 1);
  updateNeighbours(e34, index);

  // get the new triangles
  index = _SCHedges.size() - 1;
#ifdef DISPLAY_INFO
  std::cout << "New triangles: " << std::endl;
  std::cout << _SCHtriangles.size() << ' ' << _SCHedges[heap.index].vertex2 << ' ' << v3 << ' ' << v4 << ' '
            << findCircumSphere3(_SCHedges[heap.index].vertex2, v3, v4).radius << std::endl;
#endif
  _SCHtriangles.push_back(SCHtriangle(_SCHedges[heap.index].vertex2, v3, v4, index + 3, index, index + 4));

#ifdef DISPLAY_INFO
  std::cout << _SCHtriangles.size() << ' ' << _SCHedges[heap.index].vertex1 << ' ' << v4 << ' ' << v3 << ' '
            << findCircumSphere3(_SCHedges[heap.index].vertex1, v3, v4).radius << std::endl;
#endif
  _SCHtriangles.push_back(SCHtriangle(_SCHedges[heap.index].vertex1, v4, v3, index + 2, index, index + 1));

  // insert new triangles to the heap
  _heap.push(SCHheap(findCircumSphere3(_SCHedges[heap.index].vertex1, v3, v4).radius, type(triangle),
                     _SCHtriangles.size() - 1));
  _heap.push(SCHheap(findCircumSphere3(_SCHedges[heap.index].vertex2, v3, v4).radius, type(triangle),
                     _SCHtriangles.size() - 2));

#ifdef DISPLAY_INFO
  std::cout << "New edges: " << std::endl;
  std::cout << _SCHedges.size() - 1 << ' ' << v3 << ' ' << v4 << std::endl;
#endif
  // add edges to the hull and heap again
  makeNewEdge(_SCHedges[e12.first], e12.first, _SCHedges.size(), _SCHtriangles.size() - 1, heap.radius);
  makeNewEdge(_SCHedges[e12.second], e12.second, _SCHedges.size(), _SCHtriangles.size() - 1, heap.radius);
  makeNewEdge(_SCHedges[e34.first], e34.first, _SCHedges.size(), _SCHtriangles.size() - 2, heap.radius);
  makeNewEdge(_SCHedges[e34.second], e34.second, _SCHedges.size(), _SCHtriangles.size() - 2, heap.radius);

  // remove edges from hull
  _SCHedges[e12.first].removeFromHull();
  _SCHedges[e12.second].removeFromHull();
  _SCHedges[e34.first].removeFromHull();
  _SCHedges[e34.second].removeFromHull();

} // invertEdge

/**
 * @brief Creates a new edge from an old one, adds it to _SCHedges, updates its neighbours 
 * and adds it to the heap.
 * 
 * @param e the new edge
 * @param prevE index of the previous edge
 * @param newE index of the new edge
 * @param t new triangle index
 * @param maxHeap max heap radius
 */
void SchCreator3D::makeNewEdge(const SchCreator3D::SCHedge & e, size_t prevE, 
                               size_t newE, size_t t, double maxHeap)
{
  double radius;
  // add the new edge
  _SCHedges.push_back(e);
  // updates the edge's faces
  updateNeighbours(prevE, t, newE);
#ifdef DISPLAY_INFO
  std::cout << newE << ' ';
#endif
  // gets the edges circumradius
  radius = getEdgeHeap(_SCHedges[newE]);
  // checks to see if the radius is smaller than the maxHeap
  if(radius < maxHeap)
  {
    _heap.push(SCHheap(radius, type(edge), newE));
  }
  else
  {
#ifdef DISPLAY_INFO
    std::cout << "New heap is larger than current max heap" << std::endl;
#endif
  }
  // adds the new edge to its vertexes neighbours
  addToVertexNeighbours(newE);
} // makeNewHeap

/**
 * @brief Updates the face neighbours of a pair of edges
 * 
 * @param e edge index pair
 * @param index face index
 */
void SchCreator3D::updateNeighbours(SchCreator3D::SCHneighbours e, size_t index)
{
  updateNeighbours(e.first, index);
  updateNeighbours(e.second, index);
} // updateNeighbours

/**
 * @brief Updates the face neighbour of an edge
 * 
 * @param e  edge index
 * @param index face/triangle index
 */
void SchCreator3D::updateNeighbours(size_t e, size_t index)
{
  if(!_SCHtriangles[_SCHedges[e].face1].inHull)
    _SCHedges[e].face1 = index;
  else
    _SCHedges[e].face2 = index;
} // updateNeighbours

/**
 * @brief Updates an edge's face neigbours.
 * 
 * @param oldE old edge index
 * @param f face/triangle index
 * @param e new edge index
 */
void SchCreator3D::updateNeighbours(size_t oldE, size_t f, size_t e)
{
  size_t face;
  // find the face to update
  if(_SCHedges[e].face1 == f)
    face = _SCHedges[e].face2;
  else
    face = _SCHedges[e].face1;
#ifdef DISPLAY_INFO
  std::cout << "Face to update: " << face << " edge: " << e << " old edge: " << oldE << std::endl;
#endif

  // find the old edge and update it
  if(_SCHtriangles[face].e1 == oldE)
    _SCHtriangles[face].e1 = e;
  else if(_SCHtriangles[face].e2 == oldE)
    _SCHtriangles[face].e2 = e;
  else if(_SCHtriangles[face].e3 == oldE)
    _SCHtriangles[face].e3 = e;
} // updateNeighbours

/**
 * @brief Checks if two edges have the same face as neighbour
 * 
 * @param e1 edge index
 * @param e2 edge index
 * @param f1 face index
 * @return true if the edges share a face that is not f1. false otherwise
 */
bool SchCreator3D::checkSameNeighbour(size_t e1, size_t e2, size_t f1)
{
  bool check = (_SCHedges[e1].face1 == _SCHedges[e2].face1) && _SCHedges[e1].face1 != f1;
  check += (_SCHedges[e1].face1 == _SCHedges[e2].face2) && _SCHedges[e1].face1 != f1;
  check += (_SCHedges[e1].face2 == _SCHedges[e2].face1) && _SCHedges[e1].face2 != f1;
  check += (_SCHedges[e1].face2 == _SCHedges[e2].face2) && _SCHedges[e1].face2 != f1;

  return check;
}

/**
 * @brief Checks if two vertexes or edges have the same neigbour.
 * 
 * @param v1 vetex/edge index
 * @param v2 vetex/edge index
 * @param t type. edge to find vertex neigbour, triangle to find edge neighbour
 * @return true if the elements share a neighbour that is active in the hull.
 * @return false otherwise
 */
bool SchCreator3D::checkSameNeighbour(size_t v1, size_t v2, type t)
{
  if(t == 0)
  {
    // go through vertex v1 neighbours
    for(auto i = _SCHvertexes[v1].neighbours.begin(); i != _SCHvertexes[v1].neighbours.end(); i++)
    {
      // check if any of the neighbours also has v2
      if(_SCHedges[*i].vertex1 == v2 || _SCHedges[*i].vertex2 == v2)
      {
        // if the edge is not in the hull
        if(!_SCHedges[*i].inHull && i != (_SCHvertexes[v1].neighbours.end() - 1))
        {
          continue;
        }
        else
          return _SCHedges[*i].inHull;
      }
    }
  }
  else
  {
    if(_SCHedges[v1].face1 == _SCHedges[v2].face1)
    {
#ifdef DISPLAY_INFO
      std::cout << "1. Face: ";
      std::cout << _SCHedges[v1].face1 << ' ' << _SCHedges[v2].face1 << std::endl;
#endif
      return true;
    }
    else if(_SCHedges[v1].face1 == _SCHedges[v2].face2)
    {
#ifdef DISPLAY_INFO
      std::cout << "1. Face: ";
      std::cout << _SCHedges[v1].face1 << ' ' << _SCHedges[v2].face2 << std::endl;
#endif
      return true;
    }
    else if(_SCHedges[v1].face2 == _SCHedges[v2].face1)
    {
#ifdef DISPLAY_INFO
      std::cout << "1. Face: ";
      std::cout << _SCHedges[v1].face2 << ' ' << _SCHedges[v2].face1 << std::endl;
#endif
      return true;
    }
    else if(_SCHedges[v1].face2 == _SCHedges[v2].face2)
    {
#ifdef DISPLAY_INFO
      std::cout << "1. Face: ";
      std::cout << _SCHedges[v1].face2 << ' ' << _SCHedges[v2].face2 << std::endl;
#endif
      return true;
    }
  }

  return false;
} // checkSameNeighbour

/**
 * @brief Finds the pair of edges connecting v2 and v3 to v1.
 * 
 * @param f1 face index
 * @param f2 face index
 * @param v1 central vertex index
 * @param v2 vertex index
 * @param v3 vertex index
 * @return SchCreator3D::SCHneighbours 
 */
SchCreator3D::SCHneighbours SchCreator3D::findEdge(size_t f1, size_t f2, size_t v1, size_t v2, size_t v3)
{
  SCHneighbours edges;
  std::map<size_t, size_t> tempEdges = {{getEdgeKey(_SCHtriangles[f1].e1), _SCHtriangles[f1].e1},
                                        {getEdgeKey(_SCHtriangles[f1].e2), _SCHtriangles[f1].e2},
                                        {getEdgeKey(_SCHtriangles[f1].e3), _SCHtriangles[f1].e3},
                                        {getEdgeKey(_SCHtriangles[f2].e1), _SCHtriangles[f2].e1},
                                        {getEdgeKey(_SCHtriangles[f2].e2), _SCHtriangles[f2].e2},
                                        {getEdgeKey(_SCHtriangles[f2].e3), _SCHtriangles[f2].e3}};

  edges.first = tempEdges[getEdgeKey(v1, v2)];
  edges.second = tempEdges[getEdgeKey(v1, v3)];

  return edges;
} // findEdge

/**
 * @brief Finds the common vertex for two edges
 * 
 * @param e1 edge index
 * @param e2 edge index
 * @return size_t vertex index 
 */
size_t SchCreator3D::findVertex(size_t e1, size_t e2)
{
  bool inSet;
  std::set<size_t> indexes;

  // add vertexes of the first edge
  indexes.insert(_SCHedges[e1].vertex1);
  indexes.insert(_SCHedges[e1].vertex2);

  // add first vertex of second edge
  inSet = indexes.insert(_SCHedges[e2].vertex1).second;
  // if this vertex is not in the set, return the last vertex
  if(inSet) 
  {
    return _SCHedges[e2].vertex2;
  }
  // otherwise, return this vertex
  else
  {
    return _SCHedges[e2].vertex1;
  } 
} // findVertex

/**
 * @brief Finds the missing vertex in a triangle
 * 
 * @param f face index
 * @param a vertex index
 * @param b vertex index
 * @return size_t vertex index
 */
size_t SchCreator3D::findVertex(size_t f, size_t a, size_t b)
{
  bool inSet;
  std::set<size_t> indexes;

  // add known vertexes to the set
  indexes.insert(a);
  indexes.insert(b);

  // add first sphere vertex
  inSet = indexes.insert(_SCHtriangles[f].p1).second;
  // if p1 is not in the set, return p1
  if(inSet) 
  {
    return _SCHtriangles[f].p1;
  }
  // else, try adding p2 to the set
  else
  {
    inSet = indexes.insert(_SCHtriangles[f].p2).second;
  }

  // if p2 isn't in the set, return p2
  if(inSet) 
  {
    return _SCHtriangles[f].p2;
  }
  // else, return p3
  else
  {
    return _SCHtriangles[f].p3;
  }
} // insertToSet

/**
 * @brief Creates the output file
 * 
 * @param filename output filename
 */
void SchCreator3D::writeToFile(const std::string & filename)
{
  size_t index = 0;

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
  os << _bigSpheres.size() << std::endl;
  for(auto i = _bigSpheres.begin(); i != _bigSpheres.end(); i++)
  {
    // radius and center of big sphere
    os << _R << ' ' << (*i).s.center[0] << ' ';
    os << (*i).s.center[1] << ' ' << (*i).s.center[2] << std::endl;

    // coordinates of the three points touching the sphere
    os << _smallSpheres[(*i).p1].center[0] << ' ' 
       << _smallSpheres[(*i).p1].center[1] << ' '
       << _smallSpheres[(*i).p1].center[2] << ' ';
    os << _smallSpheres[(*i).p2].center[0] << ' ' 
       << _smallSpheres[(*i).p2].center[1] << ' '
       << _smallSpheres[(*i).p2].center[2] << ' ';
    os << _smallSpheres[(*i).p3].center[0] << ' ' 
       << _smallSpheres[(*i).p3].center[1] << ' '
       << _smallSpheres[(*i).p3].center[2] << std::endl;

    // normals to the sphere planes
    os << _bigSphereNormals[index].plane1.first << ' ';
    os << _bigSphereNormals[index].plane1.second[0] << ' ' 
       << _bigSphereNormals[index].plane1.second[1] << ' '
       << _bigSphereNormals[index].plane1.second[2] << std::endl;

    os << _bigSphereNormals[index].plane2.first << ' ';
    os << _bigSphereNormals[index].plane2.second[0] << ' ' 
       << _bigSphereNormals[index].plane2.second[1] << ' '
       << _bigSphereNormals[index].plane2.second[2] << std::endl;

    os << _bigSphereNormals[index].plane3.first << ' ';
    os << _bigSphereNormals[index].plane3.second[0] << ' ' 
       << _bigSphereNormals[index].plane3.second[1] << ' '
       << _bigSphereNormals[index].plane3.second[2] << std::endl;
    
    index++;
  }

  // number of torii
  os << _torii.size() << std::endl;
  for(auto i = _torii.begin(); i != _torii.end(); i++)
  {
    // remove torus boolean
    os << 1 << std::endl;

    // write torus' external radius, center and normal coordinates
    os << (*i).extRadius << ' ' << _R << ' ';

    os << (*i).center[0] << ' ' << (*i).center[1] << ' ' << (*i).center[2] << ' ';

    os << (*i).normal[0] << ' ' << (*i).normal[1] << ' ' << (*i).normal[2] << std::endl;

    // write torus' cones
    // cone [a,b]: a index, cosangle, axis coordinates
    os << _toriiCones[index].first.first << ' ';
    os << _toriiCones[index].first.second.cosangle << ' ';
    os << _toriiCones[index].first.second.axis[0] << ' ' 
       << _toriiCones[index].first.second.axis[1] << ' '
       << _toriiCones[index].first.second.axis[2] << std::endl;

    // cone [b,a]: b index, cosangle, axis coordinates
    os << _toriiCones[index].second.first << ' ';
    os << _toriiCones[index].second.second.cosangle << ' ';
    os << _toriiCones[index].second.second.axis[0] << ' ' 
       << _toriiCones[index].second.second.axis[1] << ' '
       << _toriiCones[index].second.second.axis[2] << std::endl;

    // write torus' planes
    // first triangle: triangle index, normal coordinates
    os << _toriiPlanes[index].first.first << ' ';
    os << _toriiPlanes[index].first.second[0] << ' '
       << _toriiPlanes[index].first.second[1] << ' '
       << _toriiPlanes[index].first.second[2] << std::endl;

    // second triangle: triangle index, normal coordinates
    os << _toriiPlanes[index].second.first << ' ';
    os << _toriiPlanes[index].second.second[0] << ' '
       << _toriiPlanes[index].second.second[1] << ' '
       << _toriiPlanes[index].second.second[2] << std::endl;
    
    index++;
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
    std::cout << index << '\t' << (*i).vertex1 << ' ' << (*i).vertex2 << '\t' << (*i).face1 << ' ' << (*i).face2 << '\t'
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
    if(i->inHull)
    {
      std::cout << index << '\t' << (*i).inHull << std::endl;
      std::cout << "neighbours: ";
      for(auto j = (*i).neighbours.begin(); j != (*i).neighbours.end(); j++) std::cout << *j << ' ';

      std::cout << std::endl;
    }

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
    std::cout << index << '\t' << (*i).p1 << ' ' << (*i).p2 << ' ' << (*i).p3 << "\t\t";
    std::cout << (*i).e1 << ' ' << (*i).e2 << ' ' << (*i).e3 << "   \t";
    std::cout << (*i).inHull << std::endl;
    index++;
  }
}

void SchCreator3D::printHeap()
{
  std::priority_queue<SCHheap> heap(_heap);
  size_t x = _heap.size();
  for(size_t i = 0; i < x; i++)
  {
    std::cout << i << ' ' << heap.top().radius << ' ' << heap.top().t << ' ' << heap.top().index << std::endl;
    heap.pop();
  }
}

bool SchCreator3D::checkHeap()
{
  bool check = true;
  switch(_heap.top().t)
  {
    case 0:
      check *= _SCHedges[_heap.top().index].inHull;
      check *= _SCHvertexes[_SCHedges[_heap.top().index].vertex1].inHull;
      check *= _SCHvertexes[_SCHedges[_heap.top().index].vertex2].inHull;
      check *= _SCHtriangles[_SCHedges[_heap.top().index].face1].inHull;
      check *= _SCHtriangles[_SCHedges[_heap.top().index].face2].inHull;
      break;
    default:
      check *= _SCHtriangles[_heap.top().index].inHull;
      break;
  }

  return check;
} // checkHeap

bool SchCreator3D::checkTorii(size_t e1, size_t e2)
{
  bool check;
#ifdef DISPLAY_INFO
  std::cout << "\nChecking thickness of edge " << e1 << ": " << std::endl;
#endif
  check = torusThicknessCheck(_SCHedges[e1].vertex1, _SCHedges[e1].vertex2, _SCHedges[e1].face1, _SCHedges[e1].face2);

#ifdef DISPLAY_INFO
  std::cout << check << "\nChecking thickness of edge " << e2 << ": " << std::endl;
#endif
  check *= torusThicknessCheck(_SCHedges[e2].vertex1, _SCHedges[e2].vertex2, _SCHedges[e2].face1, _SCHedges[e2].face2);
#ifdef DISPLAY_INFO
  std::cout << check << std::endl;
#endif
  return check;
}

/*  
 *  Given the indexes of two vertexex and two faces that contain them, this function computes
 *  the dot product of their plane normals.
 *
 *  If the dot product is equal or close to -1, the function returns true.
 */
bool SchCreator3D::torusThicknessCheck(size_t v1, size_t v2, size_t f1, size_t f2)
{
  Eigen::Vector3d n1 = getPlaneNormal(v1, v2, f1), n2 = getPlaneNormal(v2, v1, f2);

  double dotProduct = n1.dot(n2), d = dotProduct + 1;
  return (d * d) <= _epsilon;
} // torusThicknessCheck

/*  
 *  Given a filename, computes the SCH for the hull in the file with the desired
 *  value of alpha.
 */
void SchCreator3D::computeSCH(const std::string & filename)
{
  // read points from file with poly_algorithms
  poly.openFromFile(filename);
  // get no. of vertex
  _numberOfVertexes = poly.vertexes_.size();
  _vertexes.reserve(_numberOfVertexes);

  // ensure alpha is smaller than the max. body distance,
  if(findMaxDistance())
  {
    std::cout << "ERROR, impossible to compute SCH, choose larger R." << std::endl;
    std::cout << "Finished." << std::endl;
    return;
  }

  // get the vertexes, triangles and edges
  initialize();
  // get the initial number of active vertexes
  _activeVertexes = _vertexes.size();
  // get the initial heap
  getHeap();

#ifdef DISPLAY_INFO
  printVertexes();
  printEdges();
  printTriangles();
#endif

  // get the max heap and assign it to alpha
  SCHheap heap = _heap.top();
  _alpha = heap.radius;

#ifdef DISPLAY_INFO
  std::cout << "\n-------------------------------" << std::endl;

  std::cout << "Triangles: " << _SCHtriangles.size();
  std::cout << "\tEdges: " << _SCHedges.size() << std::endl;
  std::cout << "Max Heap: " << _alpha << ' ' << _heap.top().t << ' ' << _heap.top().index << std::endl;

#endif

  bool check; 
  // compare alpha to the desired alpha value (R-r)
  while(_alpha > _desiredAlpha)
  {

#ifdef DISPLAY_INFO
    // check that all triangles in the hull have all of its edges in the hull as well
    size_t tindex = 0;
    for(auto i = _SCHtriangles.begin(); i != _SCHtriangles.end(); i++)
    {
      if(i->inHull)
      {
        if(!_SCHedges[i->e1].inHull || !_SCHedges[i->e2].inHull || !_SCHedges[i->e3].inHull)
        {
          std::cout << "The triangle index " << tindex << ':' << i->inHull << '\t' 
                    << "Edges: " << i->e1 << ':' << _SCHedges[i->e1].inHull << ' ' 
                    << i->e2 << ':' << _SCHedges[i->e2].inHull << ' ' 
                    << i->e3 << ':' << _SCHedges[i->e3].inHull << std::endl;
          exit(1);
        }
      }
      tindex++;
    }
#endif
    
    // if the heap is empty or we have only four points left and 
    // a triangle in the max heap leave the loop
    if(_heap.size() == 0 || (_activeVertexes <= 4 && heap.t == 1)) 
      break;

    // remove the max heap
    _heap.pop();

    // check if the removed max heap is an edge or a triangle
    if(heap.t == 0)
    {
      // if it is an edge:
      //    check if the vertex disappears under the sphere
      //    or if the edge is inverted
      changeTopology(heap);
    }
    else
    {
      // if it is a triangle:
      //    disappear the vertex under the sphere
      disappearUnderEdge(heap);
    }

    // after removal check if we've reached a limit case
    if(_limitCase)
    {
#ifdef DISPLAY_INFO
      printVertexes();
      printEdges();
      printTriangles();
#endif

      std::cout << "Limit case reached. Setting new minimum value for alpha at " 
                << _desiredAlpha << "..." << std::endl;
      // set R to its new value
      _R = _r + _desiredAlpha;
      break;
    }

#ifdef DISPLAY_INFO
    std::cout << "\nCheck if current heap exists... \n";
#endif
    // check if all the elements in the new heap are in the hull
    check = checkHeap();
    while(!check)
    {
#ifdef DISPLAY_INFO
      std::cout << "Max Heap: " << _heap.top().radius << ' ' << _heap.top().t << ' ' << _heap.top().index
                << " Target alpha: " << _desiredAlpha;
      std::cout << " Check: " << check << std::endl;
#endif
      // if any of its elements is not in the hull, we remove it
      _heap.pop();
      // check again for the new heap
      check = checkHeap();
    }
#ifdef DISPLAY_INFO
    std::cout << " Done." << std::endl;
#endif

    // update the value of alpha
    heap = _heap.top();
    _alpha = heap.radius;

#ifdef DISPLAY_INFO
    printVertexes();
    printEdges();
    printTriangles();

    std::cout << "\n-------------------------------" << std::endl;

    std::cout << "Max Heap: " << _alpha << ' ' << heap.t << ' ' << heap.index << " Target alpha: " << _desiredAlpha
              << std::endl;
#endif
  }

  // get small and big spheres, torus and vertex neighbours
  getSmallSpheres();
  getBigSpheres();
  getTorii();
  getVertexNeighbours();

} // computeSCH
} // namespace sch

// ./build/src/sch-creator-3d
// /home/amrf/balloon-inflating/sch-visualization/tests/shared-tests/data/sample_polyhedron.otp
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