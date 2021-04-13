/*! \file sch-creator-2d.h
 *  \brief Declaration file of the Class sch-creator-2d
 *  \author Ana Ruiz
 *  \version 0.0.0
 *	\date 21.04.07
 */

#ifndef SCH_H
#define SCH_H

#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <memory>
#include <math.h>
#include <Eigen/Core>
#include <fstream>

/*! \namespace SCH
 *  \brief Strictly Convex Hull Namespace
 */

namespace SCH
{
    /*! \class Point
   *	\brief %Class Point
   *
   * A point structure containing the coordinates of the point and its 
   * state in the Hull (boolean)
   */

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

        Triangle(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d C)
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
   * the index of the triangle itself and the circumcircle raidus.
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
        private:
            std::vector<Eigen::Vector2d> _points;
            double _alpha; 
        public:
        // WARNING! Be sure that points entered belong to convex hull
            SchCreator2D(std::vector<Eigen::Vector2d> points, double alpha);

        public:
            std::list<Triangle> SchCreator2D::listTriangles(std::vector<Point> points, std::vector<Radius> &heap);
            void makeTriangles(const std::vector<Point> &points, std::list<Triangle> &triangles, std::vector<Radius> &heap, int previousMidpoint);
            void updateTriangles(const std::vector<Point> &points, std::list<Triangle> &triangles, std::vector<Radius> &heap);
            std::vector<Eigen::Vector2d> FindSch2D();
            bool checkHull(const std::vector<Eigen::Vector2d> &points);

        private:
            void swap(Radius &a, Radius &b);
            void insertElementToHeap(std::vector<Radius> &heap, Radius a);
            void removePointFromHull(std::vector<Point> &points, std::list<Triangle> &triangles, const Radius &heap);
            void compareToChild(std::vector<Radius> &heap);
            void removeHeap(std::vector<Radius> &heap);
            bool checkIfMaxHeapIsInHull(std::vector<Radius> &heap, std::vector<Point> points, std::list<Triangle> &triangles);
            int findPreviousPoint(const std::vector<Point> &points, int pointIndex);
            int findNextPoint(const std::vector<Point> &points, int pointIndex);
            
    };


}
#endif