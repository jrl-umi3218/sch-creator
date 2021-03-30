#include <Eigen/Core>
#include <iostream>
#include <list>
#include <math.h>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

struct Triangle
{
  Vector2d a, b, c;
  double d; // circumcircle diameter
  bool inHeap = true;

  Triangle(Vector2d A, Vector2d B, Vector2d C)
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

struct Point
{
  Vector2d point;
  bool inHull = true;

  Point(Vector2d p)
  {
    point = p;
  }

  void removeFromHull()
  {
    inHull = false;
  }
};

bool checkHull(const vector<Vector2d> & points,
               vector<Point> originalPoints,
               list<Triangle> triangles,
               const vector<Radius> & heap);

void printPoints(const vector<Point> & points)
{
  int n = 0;
  for(auto i = points.begin(); i != points.end(); i++)
  {
    if((*i).inHull)
    {
      cout << "Point " << n << endl;
      cout << (*i).point << endl;
      cout << "Length: " << (*i).point.norm() << endl;
      cout << "inHull: " << (*i).inHull << endl;
      cout << '\n' << endl;
    }

    n++;
  }
}

void printTriangles(const list<Triangle> & triangles)
{
  int n = 0;
  for(auto i = triangles.begin(); i != triangles.end(); i++)
  {
    if((*i).inHeap)
    {
      cout << "Triangle " << n << ": " << endl;
      cout << "a\t\tb\t\tc" << endl;
      cout << (*i).a.norm() << "\t\t" << (*i).b.norm() << "\t\t" << (*i).c.norm() << endl;
      cout << "Circum circle: " << (*i).d << endl;
      cout << "inHeap: " << (*i).inHeap << '\n' << endl;
    }
    n++;
  }
}

void printHeap(const vector<Radius> heap)
{
  cout << "\nHeap of cc Radius" << endl;
  for(auto i = heap.begin(); i != heap.end(); i++)
  {
    cout << "midpointIndex\ttriangleIndex\tcircumCircle Radius" << endl;
    cout << (*i).midpointIndex << "\t\t" << (*i).triangleIndex << "\t\t" << (*i).radius << endl;
  }
}

void swap(Radius * a, Radius * b)
{
  Radius temp = *a;
  *a = *b;
  *b = temp;
}

void insertElementToHeap(vector<Radius> & heap, Radius a)
{
  heap.push_back(a);
  if(heap.size() > 1)
  {
    int childIndex = heap.size();
    int parentIndex = int(childIndex / 2);

    while(heap[parentIndex - 1].radius < a.radius)
    {
      // cout << heap[parentIndex - 1].r << " < " << heap[childIndex - 1].r <<endl;
      swap(heap[parentIndex - 1], heap[childIndex - 1]);
      if(parentIndex - 2 < 0) break;
      childIndex = parentIndex;
      parentIndex = int(childIndex / 2);
    }
  }
}

// void compareToChild

list<Triangle> listTriangles(vector<Point> points, vector<Radius> & heap)
{
  list<Triangle> triangles;
  int n = points.size();

  for(int i = 0; i < n; i++)
  {
    Triangle t = Triangle(points[i % n].point, points[(i + 1) % n].point, points[(i + 2) % n].point);
    Radius r = Radius(i % n, (i + 1) % n, (i + 2) % n, i, t.d);
    triangles.push_back(t);
    insertElementToHeap(heap, r);
  }

  return triangles;
}

void removePointFromHull(vector<Point> & points, list<Triangle> & triangles, const Radius & heap)
{
  list<Triangle>::iterator it;
  int trianglesSize = triangles.size();
  int trianglesIndex = trianglesSize + heap.triangleIndex;
  for(int i = -1; i <= 1; i++)
  {
    it = triangles.begin();
    cout << "\nREMOVE FROM HULL" << endl;
    cout << "(trianglesIndex + i) >= 2 * points.size() && heap.triangleIndex < points.size()" << endl;
    cout << (trianglesIndex + i) << " >= " << 2 * points.size() << " && " << heap.triangleIndex << " < "
         << points.size() << endl;
    if((trianglesIndex + i) >= 2 * points.size() && heap.triangleIndex < points.size())
    {
      advance(it, (points.size() + heap.triangleIndex + i) % points.size());
      cout << "Triangle " << (points.size() + heap.triangleIndex + i) % points.size() << endl;
    }
    else
    {
      cout << "Triangle " << (trianglesIndex + i) % trianglesSize << endl;
      advance(it, (trianglesIndex + i) % trianglesSize);
    }

    (*it).removeFromHeap();
  }
  points[heap.midpointIndex].removeFromHull();
}

void compareToChild(vector<Radius> & heap)
{
  int parentIndex = 1;
  int childIndex = 2 * parentIndex;
  while(childIndex < heap.size() && heap[parentIndex - 1].radius < heap[childIndex - 1].radius
        || heap[parentIndex - 1].radius < heap[childIndex].radius)
  {
    if(heap[childIndex - 1].radius < heap[childIndex].radius) childIndex++;
    while(heap[parentIndex - 1].radius < heap[childIndex - 1].radius)
    {
      swap(heap[parentIndex - 1], heap[childIndex - 1]);
      parentIndex = childIndex;
      childIndex = 2 * parentIndex;
    }
  }
}

void removeHeap(vector<Radius> & heap)
{
  swap(heap.front(), heap.back());
  heap.pop_back();
  int parentIndex = 1;
  int childIndex = 2 * parentIndex;
  while(childIndex + 1 < heap.size()
        && (heap[parentIndex - 1].radius < heap[childIndex - 1].radius
            || heap[parentIndex - 1].radius < heap[childIndex].radius))
  {
    if(heap[childIndex - 1].radius < heap[childIndex].radius) childIndex++;
    if(heap[parentIndex - 1].radius < heap[childIndex - 1].radius)
    {
      swap(heap[parentIndex - 1], heap[childIndex - 1]);
      parentIndex = childIndex;
      childIndex = 2 * parentIndex;
    }
  }
}

bool checkIfHeapExists(vector<Radius> & heap, vector<Point> points, list<Triangle> & triangles)
{
  list<Triangle>::iterator it = triangles.begin();
  advance(it, heap[0].triangleIndex);
  // cout << "Triangle " << heap[0].triangleIndex << endl;
  // cout << "Circumcircle radius: " << (*it).d;
  // cout << "\tIn Heap: " << (*it).inHeap << endl;
  bool trianglePointsExist = points[heap[0].frontpointIndex].inHull
                             && (points[heap[0].midpointIndex].inHull && points[heap[0].endpointIndex].inHull);
  return (*it).inHeap && trianglePointsExist;
}

int findPreviousPoint(const vector<Point> & points, int pointIndex)
{
  while(!points[pointIndex % points.size()].inHull)
  {
    pointIndex--;
  }
  return pointIndex;
}

int findNextPoint(const vector<Point> & points, int pointIndex)
{
  while(!points[pointIndex % points.size()].inHull)
  {
    pointIndex++;
  }
  return pointIndex;
}

void makeTriangles(const vector<Point> & points,
                   list<Triangle> & triangles,
                   vector<Radius> & heap,
                   int previousMidpoint)
{
  int n = points.size();

  int newMidpoint = findPreviousPoint(points, n + previousMidpoint - 1);

  int newFrontpoint = findPreviousPoint(points, newMidpoint - 1);

  int newEndpoint = findNextPoint(points, n + previousMidpoint + 1);

  Triangle newTriangle1 =
      Triangle(points[newFrontpoint % n].point, points[newMidpoint % n].point, points[newEndpoint % n].point);
  Radius newRadius1 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, triangles.size(), newTriangle1.d);

  newFrontpoint = newMidpoint;
  newMidpoint = newEndpoint;
  newEndpoint = findNextPoint(points, newMidpoint + 1);

  Triangle newTriangle2 =
      Triangle(points[newFrontpoint % n].point, points[newMidpoint % n].point, points[newEndpoint % n].point);
  Radius newRadius2 = Radius(newFrontpoint % n, newMidpoint % n, newEndpoint % n, triangles.size() + 1, newTriangle2.d);

  triangles.push_back(newTriangle1);
  triangles.push_back(newTriangle2);
  insertElementToHeap(heap, newRadius1);
  insertElementToHeap(heap, newRadius2);
}

void updateTriangles(const vector<Point> & points, list<Triangle> & triangles, vector<Radius> & heap)
{
  int eliminatedPointIndex = heap.front().midpointIndex;
  removeHeap(heap);
  makeTriangles(points, triangles, heap, eliminatedPointIndex);
  cout << "\nUPDATE TRIANGLES\n" << endl;
  printTriangles(triangles);
  printHeap(heap);
}

vector<Vector2d> sch(vector<Point> & points)
{
  if(points.size() < 3)
  {
    cout << "You need at least 3 points.\n" << endl;
    return {};
  }
  double alpha = 5.2;
  vector<Vector2d> strictlyConvexHull;

  vector<Radius> radiusHeap;
  list<Triangle> triangles = listTriangles(points, radiusHeap);

  printPoints(points);
  printTriangles(triangles);
  printHeap(radiusHeap);

  while(radiusHeap[0].radius > alpha && checkIfHeapExists(radiusHeap, points, triangles))
  {
    cout << "\nRemove Point " << radiusHeap[0].midpointIndex << endl;
    removePointFromHull(points, triangles, radiusHeap[0]);
    updateTriangles(points, triangles, radiusHeap);

    while(!checkIfHeapExists(radiusHeap, points, triangles))
    {
      removeHeap(radiusHeap);
      cout << "\nHEAP DOESN'T EXIST. NEW HEAP" << endl;
      printHeap(radiusHeap);
    }
  }
  for(auto i = points.begin(); i != points.end(); i++)
  {
    if((*i).inHull)
    {
      strictlyConvexHull.push_back((*i).point);
    }
  }

  if(checkHull(strictlyConvexHull, points, triangles, radiusHeap))
    cout << "Convex" << endl;
  else
    cout << "Not convex" << endl;

  return strictlyConvexHull;
}

bool checkHull(const vector<Vector2d> & points,
               vector<Point> originalPoints,
               list<Triangle> triangles,
               const vector<Radius> & heap)
{
  cout << "\nCHECK IF HULL IS STRICTLY CONVEX\n" << endl;
  list<Triangle>::iterator it = triangles.begin();
  vector<Radius> realHeap;
  for(int i = 0; i < heap.size(); i++)
  {
    advance(it, heap[i].triangleIndex);
    if((*it).inHeap)
    {
      realHeap.push_back(heap[i]);
    }
    it = triangles.begin();
  }

  printHeap(realHeap);

  for(int i = 0; i < realHeap.size(); i++)
  {
    cout << "\nHull Point: " << endl;
    cout << originalPoints[realHeap[i].midpointIndex].point << '\n' << endl;
    for(int j = 0; j < originalPoints.size(); j++)
    {
      double d = (originalPoints[j].point - originalPoints[realHeap[i].midpointIndex].point).norm();
      cout << d << " > " << 2 * realHeap[i].radius << endl;
      if(d > 2 * realHeap[i].radius)
      {
        cout << "not in Hull" << endl;
        return false;
      }
      else
        cout << "in Hull" << endl;
    }
    cout << '\n' << endl;
  }
  return true;
}

int main()
{
  vector<Point> points = {Point(Vector2d(0.5, 3.25)), Point(Vector2d(7.54, 0.63)), Point(Vector2d(4.14, 4.84)),
                          Point(Vector2d(1.14, 6.36)), Point(Vector2d(1.12, 4.45))};

  vector<Vector2d> schPoints = sch(points);

  cout << "\nPoints in strictly convex Hull" << endl;
  for(auto i = schPoints.begin(); i != schPoints.end(); i++) cout << *i << '\n' << endl;

  return 0;
}