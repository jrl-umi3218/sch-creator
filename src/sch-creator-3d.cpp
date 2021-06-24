#include <Eigen/Dense> 
#include <iostream>
#include <vector>


struct Sphere
{
    Eigen::Vector3d center;
    double radius;

    Sphere(const Eigen::Vector3d &c, double r)
    {
        center = c;
        radius = r;
    }

    friend std::ostream& operator<<(std::ostream& os, const Sphere &s);
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

std::ostream& operator<<(std::ostream& os, const Sphere &s)
{
    os << "----- Sphere -----\n";
    os << "Center: \n" << s.center;
    os << "\nRadius: " << s.radius << '\n';

    return os;
}

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

// ***************** global variables ***************************
std::vector<SphereCenter> sphereCenters;

// **************************************************************


/*
*   Given 3 points, the function returns the circumsphere
*/
Sphere findCircumSphere3(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c)
{
	// find vectors a to b and a to c
	Eigen::Vector3d ab, ac;
	ab = b - a;
	ac = c- a;

	// compute dot products
	double v11 = ab.dot(ab);
	double v22 = ac.dot(ac);
	double v12 = ab.dot(ac);

	double num = v11 * v22 - pow(v12,2);
	double k1 =  0.5 * v22 * (v11 - v12) / num;
	double k2 =  0.5 * v11 * (v22 - v12) / num;

    Eigen::Vector3d center = a + k1 * ab + k2 * ac;
    double radius = (k1 * ab + k2 * ac).norm();

	return Sphere(center, radius);
}

/*
*   Given 3 points, the function finds the plane that touches them
*   and returns a matrix containing the plane's base vectors 
*/
Plane findPlaneBase(const Eigen::Vector3d &a, 
                              const Eigen::Vector3d &b,
		                      const Eigen::Vector3d &c)
{
    // find the normal to the plane
    Eigen::Vector3d ab = b - a, bc = c - b;
	Eigen::Vector3d n = ab.cross(bc).normalized();

	// find base vector ex and ey
	Eigen::Vector3d ex = ab.cross(n).normalized();
	Eigen::Vector3d ey = n.cross(ex); 
	Eigen::MatrixXd base(2,3);
	base << ex.transpose(),
		    ey.transpose();
    
    return Plane(base, n);
}

/*
*   Given 3 points, the function returns the center of the circumcircle.
*
*   NOTE: Points must be 2D. 
*/
 Eigen::Vector2d findCircleThroughPoints(const Eigen::Vector2d &a, 
                                         const Eigen::Vector2d &b,
							             const Eigen::Vector2d &c)
{
    // Using crammer's rule, find the base matrix
    Eigen::MatrixXd A(3,4);
    A << pow(a.norm(),2), a[0], a[1], 1,
         pow(b.norm(),2), b[0], b[1], 1,
         pow(c.norm(),2), c[0], c[1], 1;
    
    // From the base matrix get M1n
    Eigen::Matrix3d M11 = A.block(0,1,3,3), M12, M13, M14;
    M12 << A.leftCols(1) , A.rightCols(2); 
    M13 << A.leftCols(2) , A.rightCols(1); 
    M14 << A.leftCols(3);

    // compute the center's x and y coordinates
    Eigen::Vector2d center;
    center(0) = 0.5 * (M12.determinant() / M11.determinant());
    center(1) = -0.5 * (M13.determinant() / M11.determinant());

    return center;
}

/*
*   Given 3 points and a radius, the function returns the Sphere 
*   on whose surface lay the points. 
*/
Sphere findSphereTrhoughPoints(const Eigen::Vector3d &a, 
                              const Eigen::Vector3d &b,
							  const Eigen::Vector3d &c,
                              const double radius)
{
    // get plane base and normal
    Plane p = findPlaneBase(a,b,c);
    Eigen::MatrixXd base = p.base;
    Eigen::Vector3d n = p.normal;

	// find the 2D projection of the points using the base
	Eigen::Vector3d  ab = b - a, ac = c - a;
	Eigen::Vector2d a2d(0,0), b2d = base * ab, c2d = base * ac;

    // get the center of the circle (in 2D)
    Eigen::Vector2d circleCenter2D = findCircleThroughPoints(a2d,b2d,c2d);

    // convert coordinates to 3D
    Eigen::Vector3d circleCenter3D;
    circleCenter3D(0) = a(0) + circleCenter2D.dot(base.col(0));
    circleCenter3D(1) = a(1) + circleCenter2D.dot(base.col(1));
    circleCenter3D(2) = a(2) + circleCenter2D.dot(base.col(2));

    // get sphere and circle radius
    double circleRadius = circleCenter2D.norm();
    double sphereRadius = radius;

    // check that sphere's radius is larger than the circle's,
    // else, increase the sphere radius 
    while(sphereRadius < circleRadius) sphereRadius += 0.125;

    // get the distance from the center of the sphere
    double distanceFromCenter = sqrt(pow(sphereRadius,2) - pow(circleRadius,2));

    // find the center of the sphere
    Eigen::Vector3d sphereCenter = circleCenter3D + distanceFromCenter * n;

    // add center values to vector
    sphereCenters.push_back(SphereCenter(circleRadius, circleCenter3D, n));

    return Sphere(sphereCenter, sphereRadius);
}

/*
*   Given 4 points, the function returns the Sphere 
*   on whose surface lay the points. 
*/
Sphere findCircumSphere4(const Eigen::Vector3d a, const Eigen::Vector3d b,
                        const Eigen::Vector3d c, const Eigen::Vector3d d)
{
    // ussing crammer's rule find base matrix
    Eigen::Matrix4d T;
    T << a[0], a[1], a[2], 1,
         b[0], b[1], b[2], 1,
         c[0], c[1], c[2], 1,
         d[0], d[1], d[2], 1;
    
    // get squared norm values
    Eigen::Vector4d t;
    t << -pow(a.norm(),2), -pow(b.norm(),2), 
         -pow(c.norm(),2), -pow(d.norm(),2);

    // from base matrix get Mn
    Eigen::Matrix4d M1, M2, M3, M4;
    M1 << t, T.rightCols(3);
    M2 << T.col(0), t, T.rightCols(2);
    M3 << T.leftCols(2), t, T.col(3);
    M4 << T.leftCols(3), t;

    // get the determinants
    double Tdet = T.determinant(),
           D = M1.determinant() / Tdet,
           E = M2.determinant() / Tdet,
           F = M3.determinant() / Tdet,
           G = M4.determinant() / Tdet;

    return Sphere(Eigen::Vector3d(-D/2, -E/2, -F/2),
                  0.5*sqrt(pow(D,2)+pow(E,2)+pow(F,2)-4*G));
}


bool getDerivative(const SphereCenter &currentSphereCenter, 
                   const Eigen::Vector3d B, double radius)
{
    double circleRadius = currentSphereCenter.circleRadius;
    double distanceToCenter = sqrt(pow(radius,2) - pow(circleRadius,2));

    Eigen::Vector3d n = currentSphereCenter.planeNormal;
    Eigen::Vector3d c = currentSphereCenter.circleCenter;

    double diff = n[0] * ((c[0] + distanceToCenter * n[0]) - B[0])
                + n[1] * ((c[1] + distanceToCenter * n[1]) - B[1])
                + n[2] * ((c[2] + distanceToCenter * n[2]) - B[2]);
    
    return ((2 * radius / distanceToCenter) * diff - 2 * radius) > 0;
}

// Run:
// "/home/amrf/balloon-inflating/sch-creator/build/src/sch-creator-3d"

int main()
{
    Eigen::Vector3d a(0,5,1), c(5.25,-2,3), d(10,4.25,6.21), b(5,2,4);
    double radius = 5;
    std::cout << findSphereTrhoughPoints(a,b,c,radius);
    std::cout << findSphereTrhoughPoints(a,d,b,radius);
    std::cout << findCircumSphere3(a,b,c);
    Sphere circumSphere = findCircumSphere4(a,b,c,d);
    std::cout << circumSphere;

    return 0;
}