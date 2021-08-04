/*! \file sch-creator-2d.h
 *  \brief Declaration file of the Class sch-creator-2d
 *  \author Ana Ruiz
 *  \version 0.0.0
 *	\date 21.07.09
 */

#ifndef SCH_CREATOR_3D_H
#define SCH_CREATOR_3D_H

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
// #include <math.h>
// #include <memory>
#include <string>
#include <vector>
#include <map>
#include <sch/S_Polyhedron/Polyhedron_algorithms.h>
#include <boost/program_options.hpp>


/*! \namespace sch
 *  \brief Strictly Convex Hull Namespace
 */

namespace sch
{

  /*! \class SchCreator3D
  *	\brief %Class SchCreator3D
  *	\author Ruiz Ana
  *	\version 0.0.0
  *	\date 21.07.09
  *	\bug None
  *	\warning None
  *
  * This class read a 3d points cloud and compute its strictly convex hull made of arcs 
  */

  class SchCreator3D
  {
  public:
    struct Plane
    {
        Eigen::MatrixXd base;
        Eigen::Vector3d normal;

        Plane(const Eigen::MatrixXd &b, const Eigen::Vector3d &n)
        {
            base = b;
            normal = n;
        }

        Plane() {}
    };

    struct Sphere
    {
        Eigen::Vector3d center;
        double radius;

        Sphere(const Eigen::Vector3d &c, double r)
        {
            center = c;
            radius = r;
        }

        Sphere() {}

        friend std::ostream& operator<<(std::ostream& os, const Sphere &s);
        
    };

    struct BigSphere
    {
      Sphere s;
      size_t p1,p2,p3;

      BigSphere(const Sphere &sphere,
                size_t a, size_t b, size_t c)
      {
        s = sphere;
        p1 = a;
        p2 = b;
        p3 = c;
      }
      
      BigSphere() {}

    };

    struct SphereCenter
    {
        double circleRadius;
        Eigen::Vector3d circleCenter;
        Eigen::Vector3d planeNormal;

        SphereCenter(double r, const Eigen::Vector3d &cc, const Eigen::Vector3d &pn)
        {
            circleRadius = r;
            circleCenter = cc;
            planeNormal = pn;
        }

    };

    struct Cone
    {
      Eigen::Vector3d axis;
      double cosangle;
      
      Cone(const Eigen::Vector3d &e, double ca)
      {
        axis = e;
        cosangle = ca;
      }

      Cone() {}
    };

  public:
    SchCreator3D(double r, double R);
    
  private:
    Plane findPlaneBase(size_t a, 
                        size_t b,
                        size_t c);
    Plane findPlaneBase(size_t a, 
                        size_t b,
                        const Eigen::Vector3d &c);
    Eigen::Vector2d findCircleThroughPoints(
                        const Eigen::Vector2d &a, 
                        const Eigen::Vector2d &b,
                        const Eigen::Vector2d &c);
    bool getDerivative(const SphereCenter &currentSphereCenter, 
                        const Eigen::Vector3d B, double radius);
    Sphere findCircumSphere3(size_t a, 
                             size_t b, 
                             size_t c);
    Sphere findSphereThroughPoints(
                        size_t a, 
                        size_t b,
                        size_t c);
    Sphere findSphereThroughPoints(
                        size_t a, 
                        size_t b,
                        size_t c,
                        double radius);
    Sphere findCircumSphere4(
                        const Eigen::Vector3d &a, 
                        const Eigen::Vector3d &b,
                        const Eigen::Vector3d &c, 
                        const Eigen::Vector3d &d);
    bool findMaxDistance();  
    void getSmallSpheres();
    void getBigSpheres();
    void getCones();
    void getBigSpherePlanes();
    bool checkPointsInSphere(const Eigen::Vector3d &c, double r);
  public:
    void computeSCH(const std::string &filename);
    void writeToFile(const std::string &filename);
    void printVertexNeighbours();

    Polyhedron_algorithms poly = Polyhedron_algorithms();

  private:
    double _r, _R, _alpha, _epsilon;
    std::vector<Eigen::Vector3d> _vertexes;
    size_t _numberOfVertexes;
    std::vector<Sphere> _smallSpheres;
    std::vector<BigSphere> _bigSpheres;
    std::vector<std::vector<Cone>> _vertexNeighbours;
    std::vector<std::vector<Eigen::Vector3d>> _bigSphereNormals;
    std::vector<sch::SchCreator3D::SphereCenter> _sphereCenters;
    std::multimap<double,size_t,std::greater<double>> heap_;
  }; //class SchCreator3D

} // namespace sch

#endif