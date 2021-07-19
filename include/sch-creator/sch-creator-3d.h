/*! \file sch-creator-2d.h
 *  \brief Declaration file of the Class sch-creator-2d
 *  \author Ana Ruiz
 *  \version 0.0.0
 *	\date 21.07.09
 */

#ifndef SCH_CREATOR_3D_H
#define SCH_CREATOR_3D_H

#include <Eigen/Dense>
// #include <fstream>
#include <iostream>
// #include <math.h>
// #include <memory>
// #include <string>
#include <vector>
// #include <sch/S_Polyhedron/Polyhedron_algorithms.h>


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

    struct Plane
    {
        Eigen::MatrixXd base;
        Eigen::Vector3d normal;

        Plane(const Eigen::MatrixXd &b, const Eigen::Vector3d &n)
        {
            base = b;
            normal = n;
        }
    };

  public:
    SchCreator3D();
    
  private:
    Plane findPlaneBase(const Eigen::Vector3d &a, 
                        const Eigen::Vector3d &b,
                        const Eigen::Vector3d &c);
    Eigen::Vector2d findCircleThroughPoints(
                        const Eigen::Vector2d &a, 
                        const Eigen::Vector2d &b,
                        const Eigen::Vector2d &c);
    bool getDerivative(const SphereCenter &currentSphereCenter, 
                        const Eigen::Vector3d B, double radius);
    
  
  public:
    Sphere findCircumSphere3(Eigen::Vector3d a, 
                             Eigen::Vector3d b, 
                             Eigen::Vector3d c);
    Sphere findSphereThroughPoints(
                        const Eigen::Vector3d &a, 
                        const Eigen::Vector3d &b,
                        const Eigen::Vector3d &c,
                        const double radius);
    Sphere findCircumSphere4(
                        const Eigen::Vector3d &a, 
                        const Eigen::Vector3d &b,
                        const Eigen::Vector3d &c, 
                        const Eigen::Vector3d &d);
    void computeVVR(const std::string &filename);

  private:
    Eigen::Vector3d _points;
    std::vector<sch::SchCreator3D::SphereCenter> _sphereCenters;
    // Polyhedron_algorithms poly;
    
  }; //class SchCreator3D

} // namespace sch

#endif