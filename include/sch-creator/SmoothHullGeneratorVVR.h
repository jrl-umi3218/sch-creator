/*! \file SmoothHullGeneratorVVR.h
 *  \brief Declaration file of the Class SmoothHullGeneratorVVR
 *  \author Escande Adrien
 *  \author Cochet-Grasset Amelie
 *  \author Benallegue Mehdi
 *  \author Sylvain Miossec
 *  \version 0.0.0
 *	\date 07.10.23
 */

#pragma once

#ifndef SCD_SMOOTH_HULL_GENERATOR_VVR_H
#define SCD_SMOOTH_HULL_GENERATOR_VVR_H

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "vector3.h"


/*! \namespace SCD
 *  \brief Smooth Collision Detection Namespace
 */

namespace SCD
{
  typedef vector3<double> vector3d;

  /*! \class VVRSFace
   *	\brief %Class VVRSFace
   *	\author Cochet-Grasset Amelie
   *	\version 0.0.0
   *	\date 07.10.23
   *	\bug None
   *	\warning None
   *
   * A triangular face along with the center of its related sphere
   */
  struct VVRSFace
  {
    int			_point1;
    int			_point2;
    int			_point3;
    vector3d	_center;
  };

  /*! \class VVRsphere
   *	\brief %Class VVRsphere
   *	\author Cochet-Grasset Amelie
   *	\version 0.0.0
   *	\date 07.10.22
   *	\bug None
   *	\warning None
   *
   * A sphere
   */
  struct VVRsphere
  {
    double		_radius;
    vector3d		_center;
  };

  /*! \class VVRtorus
   *	\brief %Class VVRtorus
   *	\author Cochet-Grasset Amelie
   *	\version 0.0.0
   *	\date 07.10.22
   *	\bug None
   *	\warning None
   *
   * A Torus without the inner radius
   */
  struct VVRtorus
  {
    double		_extRadius;
    vector3d	_center;
    vector3d		_normal;
  };

  /*! \class Cone
   *	\brief %Class Cone
   *	\author Cochet-Grasset Amelie
   *	\version 0.0.0
   *	\date 07.07
   *	\bug None
   *	\warning None
   *
   * A sphere
   */
  struct Cone
  {
    double _cosangle;
    vector3d		 _axis;
  };

  typedef std::pair<int, vector3d> VVRplane;
  typedef std::pair<int, Cone> VVRcone;

  struct faceVVR
  {
    VVRplane _plane1;
    VVRplane _plane2;
    VVRplane _plane3;
  };


  /*! \class SmoothHullGeneratorVVR
   *	\brief %Class SmoothHullGeneratorVVR
   *	\author Escande Adrien
   *   \author Cochet-Grasset Amelie
   *	\version 0.0.0
   *	\date 07.10.23
   *	\bug None
   *	\warning None
   *
   * This class read a 3d points cloud and compute its smooth convex hull made of spheres and tori
   */



  class SmoothHullGeneratorVVR
  {
  private:
    struct turnData
    {
      int			_p1;				//we turn around [_p1,_p2]
      int			_p2;
      int			_p3;				//we remember the third point of the previous sphere
      vector3d		_previousCenter;	//center of the sphere we turn (/come) from
    };

  public:
    SmoothHullGeneratorVVR(double r, double R);
    SmoothHullGeneratorVVR(std::vector<vector3d>& points, double r, double R);

  public:
    //WARNING : be sure that there is no double points in the cloud
    void	loadGeometry(const std::string& filename);
    void	compute3DSMaxHull(const std::string& rootPath);
    void	computeSmoothHull(std::vector<vector3d>& outPoints, std::vector<VVRSFace>& outSFaces);
    void	computeVVR(const std::string& filename);
    void	computeVVR_WithPolyhedron(const std::string& filename);
    void	computeVVR_Prime(const std::string& filename);

  private:
    bool	findCenter(int p1, int p2, int p3, vector3d &center);
    bool	isInSphere(vector3d &point, vector3d &center);
    double	distInSphere(vector3d &point, vector3d &center);
    bool	allPointsInSphere(vector3d &center);
    double	distMaxPointsInSphere(vector3d &center);
    double	computeAngle(vector3d &p1, vector3d &p2, vector3d &p3, vector3d &axe);
    int		getKey(int a, int b);
    int		getKey(turnData& td);
    double	getKeyByAngle(turnData& td);
    void	travelCover(VVRSFace s);
    void	printSphere(VVRSFace& s);
    void	cover(void);
    void	readVertex(const std::string& filename);
    void	output(const std::string& rootPath);
    bool	findFirstTriangle(unsigned &i,unsigned &j,unsigned &k,vector3d &c);

  private:
    double							_r;
    double							_R;
    std::vector<vector3d>			_points;
    std::vector<VVRSFace>			_spheres;
    std::set<int,std::less<int> >	_index;
    double							_epsilon;
    bool							_ccw;
  };
}

#endif	//SCD_SMOOTH_HULL_GENERATOR_VVR_H

