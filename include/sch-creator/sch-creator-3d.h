/*! \file sch-creator-3d.h
 *  \brief Declaration file of the Class sch-creator-3d
 *  \author Ana Ruiz
 *  \version 0.0.0
 *	\date 21.07.09
 */

#ifndef SCH_CREATOR_3D_H
#define SCH_CREATOR_3D_H

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <queue>
#include <iomanip>
// #include <memory>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <sch/S_Polyhedron/Polyhedron_algorithms.h>
#include <boost/program_options.hpp>
#include <math.h>
#define PI 3.14159265

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

    struct Face
    {
      std::pair<size_t,Eigen::Vector3d> plane1, plane2, plane3; 
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

    struct SCHvertex
    {
      bool inHull = true;
      Eigen::Vector3d vertex;
      std::vector<size_t> neighbours;
      size_t ssIndex = -1;

      SCHvertex(const Eigen::Vector3d &v,const std::vector<size_t> &n)
      {
        vertex = v;
        neighbours = n;
      }

      void removeFromHull()
      {
        inHull = false;
      }
    };

    struct SCHedge
    {
      bool inHull = true;
      size_t vertex1, vertex2, face1, face2;
      size_t torusIndex = -1;

      SCHedge(size_t p1,size_t p2,size_t f1, size_t f2)
      {
        vertex1 = p1;
        vertex2 = p2;
        face1 = f1;
        face2 = f2;
      }

      void removeFromHull()
      {
        inHull = false;
      }
    };

    struct SCHtriangle
    {
      bool inHull = true;
      size_t p1,p2,p3;
      size_t e1,e2,e3;
      size_t bsIndex = -1;

      SCHtriangle(size_t a,size_t b,size_t c,size_t d,size_t e,size_t f)
      {
        p1 = a;
        p2 = b;
        p3 = c;
        e1 = d;
        e2 = e;
        e3 = f;
      }

      void removeFromHull()
      {
        inHull = false;
      }
    };

    enum type {edge,triangle};

    struct SCHheap
    {
      double radius;
      type t;
      size_t index;

      SCHheap(double r, type e, size_t i)
      {
        radius = r;
        t = e;
        index = i;
      }
      
      bool operator<(const SCHheap &h) const {
        return (radius < h.radius);
      }
    };

    typedef std::pair<Sphere,bool> SCHss;
    typedef std::pair<BigSphere,bool> SCHbs;
    typedef std::pair<Torus,bool> SCHt;
    typedef std::pair<size_t,Cone> SCHcone;
    typedef std::pair<size_t,Eigen::Vector3d> SCHplane;
    typedef std::pair<size_t,size_t> SCHneighbours;

  public:
    SchCreator3D(double r, double R);
    
  private:
    bool findMaxDistance(); 

    void initialize(); 

    void getSmallSpheres();

    void getBigSpheres();

    // get key functions
    size_t getEdgeKey(size_t a, size_t b);
    size_t getEdgeKey(size_t edge);
    size_t getKey(size_t a, size_t b);

    // functions required to get the big spheres
    Sphere findCircumSphere3(size_t a, size_t b, size_t c);
    Sphere findSphereThroughPoints(size_t a, size_t b,
                                   size_t c);  
    Sphere findSphereThroughPoints(size_t a, SchCreator3D::Plane p,
                                   Eigen::Vector2d circleCenter2D);                         
    Sphere findSphereThroughPoints(size_t a, SchCreator3D::Plane p,
                                   Eigen::Vector2d circleCenter2D,
                                   double R);                         
    Sphere findSphereThroughPoints(const SCHtriangle &t);
    Plane findPlaneBase(size_t a, size_t b, size_t c);
    Eigen::Vector2d findCircleThroughPoints(
                        const Eigen::Vector2d &a, 
                        const Eigen::Vector2d &b,
                        const Eigen::Vector2d &c);
    Eigen::Vector2d findCircleThroughPoints(size_t a,
                                            size_t b,
                                            size_t c,
                                            SchCreator3D::Plane p);
    
    void getTorii();
    // functions required to get the torii and plane normals
    Torus getTorus(const Sphere &s, size_t a, size_t b);
    double getCosine(size_t a, size_t b);
    std::pair<SCHcone,SCHcone> getCones(size_t a, size_t b,
                                        const Eigen::Vector3d &n);
    Eigen::Vector3d getPlaneNormal(size_t a,size_t b,Eigen::Vector3d c);
    Eigen::Vector3d getPlaneNormal(size_t a,size_t b,size_t f);
    
    void getVertexNeighbours();

    size_t findCommonFace(size_t e1, size_t e2);

    void getHeap();
    double getEdgeHeap(const SCHedge &e);
    size_t findVertex(size_t f, size_t a, size_t b);
    size_t findVertex(const SCHtriangle &t,size_t a,size_t b);

    Sphere findCircumSphere4(size_t a, size_t b, size_t c, size_t d);
    
    // other functions
    void changeTopology(SCHheap heap);
    void printEdges();
    void printVertexes();
    void printTriangles();
    void printHeap();
    void updateNeighbours(SCHneighbours e, size_t index);
    void updateNeighbours(size_t e, size_t index);
    void updateNeighbours(size_t oldE, size_t f, size_t e);
    void addToVertexNeighbours(size_t e);
    SCHneighbours findEdge(size_t f1, size_t f2, size_t v1, size_t v2, size_t v3);
    SCHneighbours findEdges(size_t f, size_t a, size_t b, size_t e);

    size_t findEdge(size_t v1, size_t v2);
    bool checkSameNeighbour(size_t v1, size_t v2, type t);
    bool checkSameNeighbour(size_t e1, size_t e2, size_t f1);
    void dissapearVertex(SCHheap heap, size_t v1, size_t v2, size_t v3 , size_t v4, 
                         SCHneighbours e12,SCHneighbours e34, size_t e);
    void invertEdge(SCHheap heap, size_t v3, size_t v4,
                    SCHneighbours e12, SCHneighbours e34);
    bool checkPointsInSphere(const Sphere &s);
    bool checkPointsInSphere(size_t a, size_t b, size_t c, size_t d);
    bool getDerivative(size_t v1, size_t v2, size_t v3, size_t v4);
    void updateVertexesIndex();
    void dissapearUnderEdge(SCHheap heap);
    void removeNeighboursFromHull(size_t v);
    size_t findActiveNeighbours(size_t v);
    bool checkHeap();
    double addNoise();
    bool torusOnSphereCheck(size_t a, size_t b, const Eigen::Vector3d &C1, 
                            const Eigen::Vector3d &C2);
    bool checkOrientation(size_t a, size_t b, size_t c);
    double angleBetween(Eigen::Vector3d a, Eigen::Vector3d b);
    void orderTriangle(size_t a, size_t &b, size_t &c);
    void swap(size_t &a, size_t &b);
    bool checkLimitCase(size_t a, size_t b, size_t c, size_t d, double r);
    void makeNewEdge(const SCHedge &e, size_t prevE, size_t newE, size_t t, double maxHeap);
    size_t findVertex(size_t e1, size_t e2);
    size_t findEdge(size_t f, size_t e1, size_t e2);
    void removeTriangle(size_t t);
    bool torusThicknessCheck(size_t v1, size_t v2, size_t f1, size_t f2);
    bool checkTorii(size_t e, size_t e1, size_t e2);
  public:
    void computeSCH(const std::string &filename);
    void writeToFile(const std::string &filename);


    Polyhedron_algorithms poly = Polyhedron_algorithms();

  private:
    double _r, _R, _alpha, _desiredAlpha, _epsilon;
    double noise = 100, maxBodyDistance = 0;

    size_t _activeVertexes;
    std::vector<SCHvertex> _SCHvertexes;
    std::vector<SCHedge> _SCHedges;
    std::vector<SCHtriangle> _SCHtriangles;
    std::map<size_t,size_t> toriiKeys;

    size_t _numberOfVertexes;
    std::vector<Eigen::Vector3d> _vertexes;
    std::vector<std::vector<SCHcone>> _vertexNeighbours;

    std::vector<Sphere> _smallSpheres;
    std::vector<BigSphere> _bigSpheres;
    std::map<size_t,Face> _bigSphereNormals;

    std::map<size_t,size_t> toriiKey;
    std::vector<Torus> _torii;
    std::map<size_t,std::pair<SCHcone,SCHcone>> _toriiCones;
    std::map<size_t,std::pair<SCHplane,SCHplane>> _toriiPlanes;

    std::priority_queue<SCHheap> _heap;
    // std::vector<SCHheap> temp;
    std::set<size_t,std::less<size_t>> difference;
    std::vector<size_t> finalIndexes;

    Eigen::Vector3d initialCenter;
    std::vector<size_t> v;
    
  }; //class SchCreator3D

} // namespace sch

#endif