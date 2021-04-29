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
#include <functional>
#include <queue>

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
      size_t frontpointIndex, midpointIndex, endpointIndex, triangleIndex;
      double radius;

      Radius(size_t fpIndex, size_t mpIndex, size_t epIndex, size_t tIndex, double R)
      {
        frontpointIndex = fpIndex;
        midpointIndex = mpIndex;
        endpointIndex = epIndex;
        triangleIndex = tIndex;
        radius = R;
      }
      bool operator<(const Radius &b) {
        return (radius < b.radius);
      }
    };

  public:
    SchCreator2D(const std::string &points);

  public:
    std::vector<Eigen::Vector2d> FindSch2D(double alpha);
    void makeTriangles(size_t previousMidpoint);
    void updateTriangles(std::list<Triangle> & triangles);
    bool checkHull(const std::vector<Eigen::Vector2d> & points);
    void readPointsFromFile();
    

  private:
    void listTriangles();
    void removePointFromHull(const Radius & heap);
    bool checkIfMaxHeapIsInHull();
    size_t findPreviousPoint(size_t pointIndex);
    size_t findNextPoint(size_t pointIndex);
    

    std::vector<Eigen::Vector2d> _points;
    std::vector<Point> _pointsStructure;
    double _alpha;
    std::string _pointsPath;
    std::priority_queue<Radius> _heap;
    std::list<SchCreator2D::Triangle> _triangles;
  }; //class SchCreator2D

} // namespace SCH

#endif