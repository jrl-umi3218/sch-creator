#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <sch/CD/CD_Scene.h>
#include <sch/STP-BV/STP_BV.h>
#include <sch/S_Object/S_ObjectNormalized.h>
#include <string>
#include <vector>

using namespace sch;

double randomNormal(double x)
{
  return exp(-0.5*x*x)/sqrt(2*M_PI);
}

int main(int argc, char ** argv)
{
  Point3 v1, v2;
  double maxDistance, maxDifference = 0;
  std::vector<Vector3> v(std::stoi(argv[1]));

  if(argc == 4)
  {
    for(size_t i = 1; i < std::stoi(argv[1]); i++)
    {
      double a = (rand() / double(RAND_MAX));
      double b = (rand() / double(RAND_MAX));
      double c = (rand() / double(RAND_MAX));

      v[i] = Vector3(randomNormal(a), 
                     randomNormal(b), 
                     randomNormal(c)).normalize();
    }

    S_ObjectNormalized * obj1 = new STP_BV();
    S_ObjectNormalized * obj2 = new STP_BV();
    obj1->constructFromFile(argv[2]);
    obj2->constructFromFile(argv[3]);

    for(size_t i = 0; i < v.size(); i++)
    {
      v1 = obj1->support(v[i]);
      v2 = obj2->support(v[i]);

      if(fabs(v1*v[i]-v2*v[i]) > maxDifference)
      {
        maxDistance = (v1 - v2).norm();
        maxDifference = fabs(v1*v[i]-v2*v[i]);
      }
    }

    std::cout << "\nMax difference: " << maxDifference 
              << " Max Distance: " << maxDistance 
              << '\n' << std::endl;
  }
  else
    std::cout << "Please provide a sample size and two filenames" << std::endl;

  return 0;
}
