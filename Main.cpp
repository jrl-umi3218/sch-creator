#include "SmoothHullGeneratorVVR.h"

#include <iostream>
#include <string.h>

#ifdef LINUX
struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '/';
    }
};
#else
#	ifdef APPLE
struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '/';
    }
};
#	else //WIN32
struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};
#	endif
#endif

std::string
basename( std::string const& pathname )
{
    return std::string( 
        std::find_if( pathname.rbegin(), pathname.rend(),
                      MatchPathSeparator() ).base(),
        pathname.end() );
}

std::string
removeExtension( std::string const& filename )
{
    std::string::const_reverse_iterator
                        pivot
            = std::find( filename.rbegin(), filename.rend(), '.' );
    return pivot == filename.rend()
        ? filename
        : std::string( filename.begin(), pivot.base() - 1 );
}

int main(int n_params, char *params[])
{
	std::string base_name;
    std::string output_dir;
    std::string r_str; //small radius
    std::string R_str; //big radius
	if (n_params<4) {
	    std::cout << "3 mandatory and 3 optional parameters needed:" << std::endl;
		std::cout << " 1. source file"<<std::endl;
		std::cout << " 2. radius of small spheres" << std::endl;
		std::cout << " 3. radius of big spheres" << std::endl;
	    std::cout << " 4. (optional) output directory" << std::endl;
		std::cout << " 5. (optional) output filename" << std::endl;
	    std::cout << " 6. (optional) polyhedron : if underlying polyhedron of STP-BV must be also generated (default : no)" << std::endl;
	    return -1;
	}

    r_str = params[2];
    R_str = params[3];
	std::string source_file(params[1]);

	

	if (n_params>4) 
		output_dir=std::string(params[4])+"/";
	else
		output_dir="";

	if (n_params>5) 
		base_name=std::string(params[5]);
	else
	{
		base_name=basename(source_file);
		base_name=removeExtension(base_name);
	}
	
    std::string output_file(output_dir+base_name+"_r"+r_str+"_R"+R_str+".txt");
    double r = strtod(r_str.c_str(), NULL);
	double R = strtod(R_str.c_str(), NULL);

    bool gen_polyhedron = false;
    if (n_params>6) {
        if (strcmp("polyhedron",params[6])==0)
            gen_polyhedron = true;
        else {
            std::cout << "ERROR, 6th parameter must be ""polyhedron"" to generate it" << std::endl;
            return -1;
        }
    }

    SCD::SmoothHullGeneratorVVR STPBV_Generator(r,R);
    STPBV_Generator.loadGeometry(source_file);
    if (gen_polyhedron)
        STPBV_Generator.computeVVR_WithPolyhedron(output_file);
    else
        STPBV_Generator.computeVVR_Prime(output_file);

}
