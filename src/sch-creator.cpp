#include <sch-creator/SmoothHullGeneratorVVR.h>
#include <iostream>
#include <vector>
#include <string>

#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
  double r, R;

  po::options_description desc("Allowed options");
  desc.add_options()
  ("help,h", "produce help message")
  ("r,r", po::value<double>(&r)->default_value(.2), "small sphere radius")
  ("R,R", po::value<double>(&R)->default_value(300.), "big sphere radius")
  ("input-file", po::value<string>(), "input file")
  ("output-file", po::value<string>(), "output file")
  ("poly", po::value<bool>()->default_value(false), "generate the polyhedron of the STP-BV");

  po::positional_options_description pos;
  pos.add("input-file", 1);
  pos.add("output-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    cout << desc << endl;
    return 1;
  }

  cout << "\n STP-BV parameters: r = " << r << ", R = " << R << std::endl << std::endl;

  if (vm.count("input-file") && vm.count("output-file"))
  {
    SCD::SmoothHullGeneratorVVR sg(r, R);
    fstream testfile;
    string input = vm["input-file"].as<string>();
    string output = vm["output-file"].as<string>();

    testfile.open(input.c_str());

    if (testfile.is_open())
    {
      cout << "Opening "<< input << endl;

      testfile.close();

      sg.loadGeometry(input);
      if(vm["poly"].as<bool>())
        sg.computeVVR_WithPolyhedron(output);
      else
        sg.computeVVR_Prime(output);

      cout << "STP-BV Created, output file "<< output << endl;
      cout << "Successfully finished" << endl;
    }
    else
      cout << "Failed to open " << input << endl;
  }

  return 0;
}
