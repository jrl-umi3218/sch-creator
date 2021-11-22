/*! \file sch-creator-3d.h
 *  \brief Declaration file of the Class sch-creator-3d
 *  \author Ana Ruiz
 *  \version 1.0.0
 *	\date 17.11.21
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
  *	\version 1.0.0
  *	\date 17.11.21
  *	\bug Output files show errors in visualizer
  *	\warning Opening too many files in visualizer crashes vscode
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

  // Structures used to describe the SCH

    /**
     * @brief Describes a regular sphere.
     * 
     */
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

    /**
     * @brief Describes one of the spheres that touch 3 vertexes.
     * 
     */
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

    /**
     * @brief Contains the three plane normals limiting a big sphere.
     * 
     */
    struct Face
    {
      std::pair<size_t,Eigen::Vector3d> plane1, plane2, plane3; 
    };

    /**
     * @brief A torus cone. The axis should be the same as its corresponding
     * torus' normal, which should be the normalized edge.
     * 
     */
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

    /**
     * @brief A torus. The normal corresponds to the "axis" of its cone.
     * 
     */
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

  // Structures used to build the SCH

    /**
     * @brief Contains the state if the current vertex in the hull,
     * a vector indicating its edge neighbours and the index of the final
     * corresponding small sphere.
     * 
     */
    struct SCHvertex
    {
      bool inHull = true;
      std::vector<size_t> neighbours;
      size_t ssIndex = -1;

      SCHvertex(const std::vector<size_t> &n)
      {
        neighbours = n;
      }

      void removeFromHull()
      {
        inHull = false;
      }
    };  

    /**
     * @brief Contains the state of the current edge in the hull,
     * the indexes of the vertexes of the edge, the neighbour faces/triangles
     * to the edge and the final torus index.
     * 
     */
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

    /**
     * @brief Contains the state of the current triangle in the hull, 
     * the three vertexes that make up the triangle, the three edges
     * that are neighbours to the triangle and the final big sphere index.
     * 
     */
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

    /**
     * @brief type is used in the heap to indicate whether it corresponds to
     * an edge or a trianlge
     * 
     */
    enum type {edge,triangle};

    /**
     * @brief Contains the computed circumradius, the type of element that produced
     * this radius and the index of said element.
     * 
     */
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

    /**
     * @brief The first element indicates the index of the cone's initial vertex
     * (according to the direction of its axis), while the second element is the
     * cone itself.
     * 
     */
    typedef std::pair<size_t,Cone> SCHcone;

    /**
     * @brief The first element indicates the index of the triangle to which the
     * plane normal belongs to, while the second element is the plane normal itself.
     * 
     */
    typedef std::pair<size_t,Eigen::Vector3d> SCHplane;

    /**
     * @brief A pair of indexes that will point to the pair of neighbours
     * that fulfill the required conditions.
     * 
     */
    typedef std::pair<size_t,size_t> SCHneighbours;

  public:
    SchCreator3D(double r, double R);
    
  private:
    //  heap functions
    void getHeap();
    double getEdgeHeap(const SCHedge &e);
    bool checkHeap();

    // get key functions
    size_t getEdgeKey(size_t a, size_t b);
    size_t getEdgeKey(size_t edge);
    size_t getKey(size_t a, size_t b);

    // sphere functions
      // circumsphere
    Sphere findCircumSphere3(size_t a, size_t b, size_t c);
    Sphere findCircumSphere4();
    Sphere findCircumSphere4(size_t a, size_t b, size_t c, size_t d);
      // sphere through points
    Sphere findSphereThroughPoints(size_t a, size_t b,
                                   size_t c);  
    Sphere findSphereThroughPoints(size_t a, SchCreator3D::Plane p,
                                   Eigen::Vector2d circleCenter2D);                         
    Sphere findSphereThroughPoints(size_t a, SchCreator3D::Plane p,
                                   Eigen::Vector2d circleCenter2D,
                                   double R);                         
    Sphere findSphereThroughPoints(const SCHtriangle &t);
      // find circle
    Eigen::Vector2d findCircleThroughPoints(
                        const Eigen::Vector2d &a, 
                        const Eigen::Vector2d &b,
                        const Eigen::Vector2d &c);
    Eigen::Vector2d findCircleThroughPoints(size_t a,
                                            size_t b,
                                            size_t c,
                                            SchCreator3D::Plane p);
      //  get plane
    Plane findPlaneBase(size_t a, size_t b, size_t c);
      //  check in sphere
    bool checkPointsInSphere(const Sphere &s);
    bool checkPointsInSphere(size_t a, size_t b, size_t c, size_t d);
    

    //  change of topology functions
    void changeTopology(SCHheap heap);
    void disappearVertex(SCHheap heap, size_t v1, size_t v2, size_t v3 , size_t v4, 
                         SCHneighbours e12,SCHneighbours e34, size_t e);
    void invertEdge(SCHheap heap, size_t v3, size_t v4,
                    SCHneighbours e12, SCHneighbours e34);
    void disappearUnderEdge(SCHheap heap);
    void makeNewEdge(const SCHedge &e, size_t prevE, size_t newE, size_t t, double maxHeap);

    // update functions
    void updateNeighbours(SCHneighbours e, size_t index);
    void updateNeighbours(size_t e, size_t index);
    void updateNeighbours(size_t oldE, size_t f, size_t e);

    // build hull functions
    void getSmallSpheres();
    void getBigSpheres();
    void getTorii();
    void getVertexNeighbours();
        
    // functions required to get the torii and plane normals
    Torus getTorus(const Sphere &s, size_t a, size_t b);
    double getCosine(size_t a, size_t b);
    std::pair<SCHcone,SCHcone> getCones(size_t a, size_t b,
                                        const Eigen::Vector3d &n);
    Eigen::Vector3d getPlaneNormal(size_t a,size_t b,Eigen::Vector3d c);
    Eigen::Vector3d getPlaneNormal(size_t a,size_t b,size_t f);

    // print functions
    void printEdges();
    void printVertexes();
    void printTriangles();
    void printHeap();

    // find functions
    SCHneighbours findEdge(size_t f1, size_t f2, size_t v1, size_t v2, size_t v3);
    SCHneighbours findEdges(size_t f, size_t a, size_t b, size_t e);
    size_t findEdge(size_t v1, size_t v2);
    size_t findEdge(size_t f, size_t e1, size_t e2);

    size_t findCommonFace(size_t e1, size_t e2);

    size_t findVertex(size_t e1, size_t e2);
    size_t findVertex(size_t f, size_t a, size_t b);
    size_t findVertex(const SCHtriangle &t,size_t a,size_t b);

    
    // swap
    void swap(size_t &a, size_t &b);
    void swap(SCHneighbours &a, SCHneighbours &b);
    
    // check neighbours
    bool checkSameNeighbour(size_t v1, size_t v2, type t);
    bool checkSameNeighbour(size_t e1, size_t e2, size_t f1);

    void removeNeighboursFromHull(size_t v);

    size_t findActiveNeighbours(size_t v);

    void addToVertexNeighbours(size_t e);

    // torus checks
    bool torusThicknessCheck(size_t v1, size_t v2, size_t f1, size_t f2);
    bool checkTorii(size_t e1, size_t e2);

    // other
    bool findMaxDistance(); 

    void initialize(); 

    bool getDerivative(size_t v1, size_t v2, size_t v3, size_t v4);
    
    double addNoise();

    bool checkOrientation(size_t a, size_t b, size_t c);

    void orderTriangle(size_t a, size_t &b, size_t &c);

    bool checkLimitCase(size_t a, size_t b, size_t c, size_t d, double r);
  public:
    void computeSCH(const std::string &filename);
    void writeToFile(const std::string &filename);


    Polyhedron_algorithms poly = Polyhedron_algorithms();

  private:
    double _r, _R, _alpha, _desiredAlpha, _epsilon;
    double noise = 100, maxBodyDistance = 0;
    bool _limitCase = false;

    size_t _activeVertexes;
    std::vector<SCHvertex> _SCHvertexes;
    std::vector<SCHedge> _SCHedges;
    std::vector<SCHtriangle> _SCHtriangles;

    size_t _numberOfVertexes;
    std::vector<Eigen::Vector3d> _vertexes;
    std::vector<std::vector<SCHcone>> _vertexNeighbours;

    std::vector<Sphere> _smallSpheres;
    
    std::vector<BigSphere> _bigSpheres;
    std::map<size_t,Face> _bigSphereNormals;

    std::vector<Torus> _torii;
    std::map<size_t,std::pair<SCHcone,SCHcone>> _toriiCones;
    std::map<size_t,std::pair<SCHplane,SCHplane>> _toriiPlanes;

    std::priority_queue<SCHheap> _heap;

    Eigen::Vector3d initialCenter;
    
  }; //class SchCreator3D

} // namespace sch

#endif