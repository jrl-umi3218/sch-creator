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
        size_t p1, p2, p3;

        Plane(const Eigen::MatrixXd &b, const Eigen::Vector3d &n,
              size_t A, size_t B, size_t C)
        {
            base = b;
            normal = n;
            p1 = A;
            p2 = B;
            p3 = C;
        }

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

    struct Torus
    {
      Eigen::Vector3d center, normal;
      double extRadius;

      Torus(const Eigen::Vector3d &c, const Eigen::Vector3d &n, double r)
      {
        center = c;
        normal = n;
        extRadius = r;
      }

      Torus() {}
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
                        size_t a, 
                        size_t b,
                        size_t c, 
                        size_t d);
    bool findMaxDistance();  
    void getSmallSpheres();
    void getBigSpheres();
    void getCones();
    Torus getTorus(const Sphere &s, size_t a, size_t b);
    void removeUselessTorii();
    void getEdgeNeighbours();
    void getTorusNeighbours();
    void getBigSpherePlanes();
    void getBigSphereEdges();
    size_t getEdgeKey(size_t a, size_t b);
    bool checkPointsInSphere(const Eigen::Vector3d &c, double r);
  public:
    void computeSCH(const std::string &filename);
    void writeToFile(const std::string &filename);
    void printVertexNeighbours();
    void printSphereEdges();
    void printBigSpherePlanes();

    Polyhedron_algorithms poly = Polyhedron_algorithms();

  private:
    double _r, _R, _alpha, _epsilon;
    std::vector<Eigen::Vector3d> _vertexes;
    size_t _numberOfVertexes;
    std::vector<Sphere> _smallSpheres;
    std::vector<BigSphere> _bigSpheres;
    std::vector<Cone> _cones;
    std::vector<std::pair<size_t,Torus>> _torii;
    std::vector<std::vector<size_t>> _vertexNeighbours;
    std::vector<std::pair<std::pair<size_t,size_t>,std::pair<size_t,size_t>>> _edgeNeighbours;
    std::vector<std::vector<std::pair<size_t,Eigen::Vector3d>>> _bigSphereNormals;
    std::vector<std::vector<size_t>> _bigSphereEdgess;
    std::map<size_t, size_t> ids;
    std::map<size_t,std::pair<std::pair<size_t,size_t>,std::pair<size_t,size_t>>> _toriiNeighbours;
    std::vector<sch::SchCreator3D::SphereCenter> _sphereCenters;
    std::map<double,size_t,std::greater<double>> heap_;
  }; //class SchCreator3D

} // namespace sch

#endif