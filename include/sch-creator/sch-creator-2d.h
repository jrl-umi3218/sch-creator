/*! \file sch-creator-2d.h
 *  \brief Declaration file of the Class sch-creator-2d
 *  \author Ana Ruiz
 *  \version 0.0.0
 *	\date 21.04.07
 */

#ifndef SCH_CREATOR_2D_H
#define SCH_CREATOR_2D_H

#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <memory>
#include <string>
#include <vector>

/*! \namespace SCH
 *  \brief Strictly Convex Hull Namespace
 */

namespace SCH
{

/*! \class SchCreator2D
 *	\brief %Class SchCreator2D
 *	\author Escande Adrien
 *   \author Cochet-Grasset Amelie
 *	\version 0.0.0
 *	\date 07.10.23
 *	\bug None
 *	\warning None
 *
 * This class read a 3d points cloud and compute its smooth convex hull made of spheres and tori
 */

class SchCreator2D
{
  /*! \class Point
   *	\brief %Class Point
   *
   * A point structure containing the coordinates of the point and its
   * state in the Hull (boolean)
   */
public:
  struct Point
  {
    Eigen::Vector2d point;
    bool inHull = true;

    Point(Eigen::Vector2d p)
    {
      point = p;
    }

    void removeFromHull()
    {
      inHull = false;
    }
  };

  /*! \class Triangle
   *	\brief %Class Triangle
   *
   * A triangle structure containing the vertices of the triangle, its
   * corresponding Circumcircle radius and its state in the heap.
   */

  struct Triangle
  {
    Eigen::Vector2d a, b, c;
    double d; // circumcircle diameter
    bool inHeap = true;

    Triangle(const Eigen::Vector2d & A, const Eigen::Vector2d & B, const Eigen::Vector2d & C)
    {
      a = A;
      b = B;
      c = C;
      d = findCircumcircleRadius();
    }

    double findCircumcircleRadius()
    {
      double A = (a - b).norm();
      double B = (b - c).norm();
      double C = (a - c).norm();
      double s = (A + B + C) / 2;
      return (A * B * C) / (4 * sqrt(s * (s - A) * (s - B) * (s - C)));
    }

    void removeFromHeap()
    {
      inHeap = false;
    }
  };

  /*! \class Radius
   *	\brief %Class Radius
   *
   * A radius structure containing the indexes of the triangle's vertices,
   * the index of the triangle itself and the circumcircle radius.
   */

  struct Radius
  {
    int frontpointIndex, midpointIndex, endpointIndex, triangleIndex;
    double radius;

    Radius(int fpIndex, int mpIndex, int epIndex, int tIndex, double R)
    {
      frontpointIndex = fpIndex;
      midpointIndex = mpIndex;
      endpointIndex = epIndex;
      triangleIndex = tIndex;
      radius = R;
    }
  };

public:
  SchCreator2D(const std::string &points);

public:
  std::vector<Eigen::Vector2d> FindSch2D(double alpha);
  void makeTriangles(const std::vector<Point> & points,
                     std::list<Triangle> & triangles,
                     std::vector<Radius> & heap,
                     size_t previousMidpoint);
  void updateTriangles(const std::vector<Point> & points, std::list<Triangle> & triangles, std::vector<Radius> & heap);
  bool checkHull(const std::vector<Eigen::Vector2d> & points);
  void readPointsFromFile();
  

private:
  void swap(Radius & a, Radius & b);
  std::list<Triangle> listTriangles(std::vector<Point> & points, std::vector<Radius> & heap);
  void insertElementToHeap(std::vector<Radius> & heap, Radius & a);
  void removePointFromHull(std::vector<Point> & points, std::list<Triangle> & triangles, const Radius & heap);
  void compareToChild(std::vector<Radius> & heap);
  void removeHeap(std::vector<Radius> & heap);
  bool checkIfMaxHeapIsInHull(std::vector<Radius> & heap, std::vector<Point> & points, std::list<Triangle> & triangles);
  size_t findPreviousPoint(const std::vector<Point> & points, size_t pointIndex);
  size_t findNextPoint(const std::vector<Point> & points, size_t pointIndex);
  

  std::vector<Eigen::Vector2d> _points;
  double _alpha;
  std::string _pointsPath;
};

} // namespace SCH
#endif