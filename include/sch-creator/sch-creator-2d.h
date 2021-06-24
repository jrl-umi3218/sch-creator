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
#include <yaml-cpp/yaml.h>
// #include <sch/S_Object/S_Sphere.h>

/*! \namespace sch
 *  \brief Strictly Convex Hull Namespace
 */

namespace sch
{

  /*! \class SchCreator2D
  *	\brief %Class SchCreator2D
  *	\author Ruiz Ana
  *	\version 0.0.0
  *	\date 21.04.07
  *	\bug None
  *	\warning None
  *
  * This class read a 2d points cloud and compute its strictly convex hull made of arcs 
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
      double r; // circumcircle diameter
      bool inHeap = true;

      Triangle(const Eigen::Vector2d & A, const Eigen::Vector2d & B, const Eigen::Vector2d & C)
      {
        a = A;
        b = B;
        c = C;
        r = findCircumcircleRadius();
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
    * A structure containing the indexes of the triangle's vertices,
    * the index of the triangle itself and the circumcircle radius.
    */

    struct Radius
    {
      size_t startpointIndex, midpointIndex, 
             endpointIndex, prevTriangleIndex,
             triangleIndex, nextTriangleIndex;
      double radius;

      Radius(size_t spIndex, size_t mpIndex, 
             size_t epIndex, size_t prevtIndex, 
             size_t tIndex, size_t nexttIndex, double R)
      {
        startpointIndex = spIndex;
        midpointIndex = mpIndex;
        endpointIndex = epIndex;
        prevTriangleIndex = prevtIndex;
        triangleIndex = tIndex;
        nextTriangleIndex = nexttIndex;
        radius = R;
      }

      Radius() {}
      
      bool operator<(const Radius &b) const {
        return (radius < b.radius);
      }
    };

  public:
    SchCreator2D(const std::string &points);

  public:
    void FindSch2D(double alpha);
    bool checkHull();
    void printTriangles();
    void printPoints();

  private:
    void listTriangles();
    void makeTriangles(size_t previousMidpoint);
    void updateTriangles(std::list<Triangle> & triangles);
    void removePointFromHull(const Radius & heap);
    bool checkIfMaxHeapIsInHull();
    size_t findPreviousPoint(size_t pointIndex);
    size_t findNextPoint(size_t pointIndex);
    size_t findPreviousTriangle(size_t triangleIndex);
    size_t findNextTriangle(size_t triangleIndex);
    void findNewAlpha();
    void readPointsFromFile();
    void makeYAML();
    void findClosestPoint();
  
  public:
    std::vector<Eigen::Vector2d> _points;
    std::vector<Eigen::Vector2d> _schPoints;
    std::vector<Point> _pointsStructure;
    size_t _eliminatedPoints;
    double _alpha;
    double _initialAlpha;
    std::string _pointsPath;
    std::priority_queue<Radius> _heap;
    std::vector<Radius> _eliminatedVertex;
    std::list<SchCreator2D::Triangle> _triangles;
  }; //class SchCreator2D

} // namespace sch

#endif