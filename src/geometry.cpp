// Contains functions to compute geometric quantities (e.g. area, normal,
// centroid, spacing functions)

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>
#include "transformations.h"

/******************************************************************************/
//
// Computes area of triangle from three nodes given in counterclockwise
// order
//
/******************************************************************************/
double tri_area ( const double & x0, const double & y0, const double & z0,
                  const double & x1, const double & y1, const double & z1,
                  const double & x2, const double & y2, const double & z2 )
{
  double dxa, dya, dza, dxb, dyb, dzb;
  double sx, sy, sz;

  // Edge vectors

  dxa = x1 - x0;
  dya = y1 - y0;
  dza = z1 - z0;

  dxb = x2 - x0;
  dyb = y2 - y0;
  dzb = z2 - z0;

  // Face vector

  sx = dya*dzb - dza*dyb;
  sy = dza*dxb - dxa*dzb;
  sz = dxa*dyb - dya*dxb;

  // Area (1/2 magnitude of face vector)

  return 0.5*std::sqrt(sx*sx + sy*sy + sz*sz);
}

/******************************************************************************/
//
// Computes area of quadrilateral from four nodes given in counterclockwise
// order
//
/******************************************************************************/
double quad_area ( const double & x0, const double & y0, const double & z0,
                   const double & x1, const double & y1, const double & z1,
                   const double & x2, const double & y2, const double & z2,
                   const double & x3, const double & y3, const double & z3 )
{
  double dxa, dya, dza, dxb, dyb, dzb;
  double sx, sy, sz;

  // Vectors crossing the face

  dxa = x2 - x0;
  dya = y2 - y0;
  dza = z2 - z0;

  dxb = x3 - x1;
  dyb = y3 - y1;
  dzb = z3 - z1;

  // Face vector

  sx = dya*dzb - dza*dyb;
  sy = dza*dxb - dxa*dzb;
  sz = dxa*dyb - dya*dxb;

  // Area (1/2 magnitude of face vector)

  return 0.5*std::sqrt(sx*sx + sy*sy + sz*sz);
}

/******************************************************************************/
//
// Computes normal for triangle from three nodes given in counterclockwise 
//   order
//
/******************************************************************************/
Eigen::Vector3d tri_normal ( 
                       const double & x0, const double & y0, const double & z0,
                       const double & x1, const double & y1, const double & z1,
                       const double & x2, const double & y2, const double & z2 )
{
  double dxa, dya, dza, dxb, dyb, dzb;
  double length, sx, sy, sz;
  Eigen::Vector3d norm;

  // Edge vectors

  dxa = x1 - x0;
  dya = y1 - y0;
  dza = z1 - z0;

  dxb = x2 - x0;
  dyb = y2 - y0;
  dzb = z2 - z0;

  // Face vector

  sx = dya*dzb - dza*dyb;
  sy = dza*dxb - dxa*dzb;
  sz = dxa*dyb - dya*dxb;

  // Magnitude of face vector

  length = sqrt(sx*sx + sy*sy + sz*sz);

  // Normal vector

  norm(0) = sx/length;
  norm(1) = sy/length;
  norm(2) = sz/length;

  return norm;
}

/******************************************************************************/
//
// Computes normal for quadrilateral from four nodes given in counterclockwise 
//   order
//
/******************************************************************************/
Eigen::Vector3d quad_normal ( 
                       const double & x0, const double & y0, const double & z0,
                       const double & x1, const double & y1, const double & z1,
                       const double & x2, const double & y2, const double & z2,
                       const double & x3, const double & y3, const double & z3 )
{
  double dxa, dya, dza, dxb, dyb, dzb;
  double length, sx, sy, sz;
  Eigen::Vector3d norm;

  // Vectors crossing the face

  dxa = x2 - x0;
  dya = y2 - y0;
  dza = z2 - z0;

  dxb = x3 - x1;
  dyb = y3 - y1;
  dzb = z3 - z1;

  // Face vector

  sx = dya*dzb - dza*dyb;
  sy = dza*dxb - dxa*dzb;
  sz = dxa*dyb - dya*dxb;

  // Magnitude of face vector

  length = sqrt(sx*sx + sy*sy + sz*sz);

  // Normal vector

  norm(0) = sx/length;
  norm(1) = sy/length;
  norm(2) = sz/length;

  return norm;
}

/******************************************************************************/
//
// Computes centroid of triangle from three nodes given in counterclockwise
// order via intersection of vertex-midpoint lines
//
/******************************************************************************/
Eigen::Vector3d tri_centroid (
                       const double & x0, const double & y0, const double & z0,
                       const double & x1, const double & y1, const double & z1,
                       const double & x2, const double & y2, const double & z2 )
{
  unsigned int i, mincompa, mincompb, mincomp;
  double lanorm, lbnorm, mina, minb, r0, r1, ta;
  bool eqswitch;
  Eigen::Vector3d p0, p1, p2, pa, pb, la, lb, latemp, lbtemp;
  Eigen::Vector2d rhsvec, tvec;
  Eigen::Matrix2d lhsmat;

  // Construct node vectors

  p0(0) = x0;
  p0(1) = y0;
  p0(2) = z0;
  p1(0) = x1;
  p1(1) = y1;
  p1(2) = z1;
  p2(0) = x2;
  p2(1) = y2;
  p2(2) = z2;

  // Edge midpoints

  pa = 0.5*(p0 + p1);
  pb = 0.5*(p0 + p2);

  // Construct lines between pa, pb and opposite vertices, e.g.,
  // line_a = pa + unit_vector(p2 - pa)*ta

  la = p2 - pa;
  lanorm = la.norm();
  la = la/lanorm;

  lb = p1 - pb;
  lbnorm = lb.norm();
  lb = lb/lbnorm;

  // Get intersection point between lines.  The system has three equations and
  // two unknowns (ta and tb).  Thus discard the one with the smallest component
  // on the vectors la and lb to ensure the system is not ill-conditioned.

  mina = 1000.;
  minb = 1000.;
  for ( i = 0; i < 3; i++ )
  {
    latemp(i) = std::abs(la(i));
    if (latemp(i) < mina) 
    {
      mina = latemp(i); 
      mincompa = i; 
    }
    lbtemp(i) = std::abs(lb(i));
    if (lbtemp(i) < minb) 
    {
      minb = lbtemp(i); 
      mincompb = i; 
    }
  }

  if (latemp(mincompa) < lbtemp(mincompb)) { mincomp = mincompa; }
  else { mincomp = mincompb; }

  // Set up system of equations

  if (mincomp == 0)
  {
    lhsmat(0,0) = la(1); lhsmat(0,1) = -lb(1); rhsvec(0) = pb(1) - pa(1);
    lhsmat(1,0) = la(2); lhsmat(1,1) = -lb(2); rhsvec(1) = pb(2) - pa(2);
  }
  else if (mincomp == 1)
  {
    lhsmat(0,0) = la(0); lhsmat(0,1) = -lb(0); rhsvec(0) = pb(0) - pa(0);
    lhsmat(1,0) = la(2); lhsmat(1,1) = -lb(2); rhsvec(1) = pb(2) - pa(2);
  }
  else
  {
    lhsmat(0,0) = la(0); lhsmat(0,1) = -lb(0); rhsvec(0) = pb(0) - pa(0);
    lhsmat(1,0) = la(1); lhsmat(1,1) = -lb(1); rhsvec(1) = pb(1) - pa(1);
  }

  // At this point it is still possible to have a singular matrix if the two
  // equations are not linearly independent.  If this is the case, switch out
  // the first equation with the unused one.

  r0 = lhsmat(1,0) / lhsmat(0,0);
  r1 = lhsmat(1,1) / lhsmat(0,1);
  eqswitch = false;
  if (r0 == r1) { eqswitch = true; }
  if (std::abs(r0 - r1) < 1.E-08)
  {
    if ( (std::abs(la(mincomp)) > 1.E-08) && std::abs(lb(mincomp) > 1.E-08) )
    {
      eqswitch = true;
    }
  }

  if (eqswitch)
  {
    lhsmat(0,0) = la(mincomp);
    lhsmat(0,1) = -lb(mincomp);
    rhsvec(0) = pb(mincomp) - pa(mincomp);
  }

  // Solve system of equations for ta and tb

  tvec = lhsmat.lu().solve(rhsvec);

  // Get centroid using intersection on line a

  ta = tvec(0);
  return pa + la*ta;
}

/******************************************************************************/
//
// Computes centroid of quadrialteral from four nodes given in counterclockwise
// order
//
/******************************************************************************/
Eigen::Vector3d quad_centroid (
                       const double & x0, const double & y0, const double & z0,
                       const double & x1, const double & y1, const double & z1,
                       const double & x2, const double & y2, const double & z2,
                       const double & x3, const double & y3, const double & z3 )
{
  double area1, area2;
  Eigen::Vector3d cen1, cen2;

  // Compute area and centroid of triangle 1

  area1 = tri_area(x0, y0, z0, x1, y1, z1, x2, y2, z2);
  cen1 = tri_centroid(x0, y0, z0, x1, y1, z1, x2, y2, z2);

  // Compute area and centroid of triangle 2

  area2 = tri_area(x0, y0, z0, x2, y2, z2, x3, y3, z3);
  cen2 = tri_centroid(x0, y0, z0, x2, y2, z2, x3, y3, z3);

  // Quadrilateral centroid

  return (cen1*area1 + cen2*area2) / (area1 + area2);
}

/******************************************************************************/
//
// Computes inertial->panel transformation from normal vector
//
/******************************************************************************/
Eigen::Matrix3d transform_from_normal ( const double & nx, const double & ny,
                                        const double & nz )
{
  double term, base, a1, a2;

  // Special case where normal is in +y or -y direction

  if (ny == 1.0)
  {
    a1 = 0.0;
    a2 = M_PI/2.0;
  }
  else if (ny == -1.0)
  {
    a1 = 0.0;
    a2 = -M_PI/2.0;
  }

  else
  {
    // Angle between z axis and projection of normal on x-z plane
  
    term = nz / std::sqrt(nx*nx + nz*nz);
    if (term > 1.0) { term = 1.0; }         // Rounding error fix
    else if (term < -1.0) { term = -1.0; }

    base = std::acos(term);
    if (nx <= 0.0) { a1 = base; }
    else { a1 = -base; }
  
    // Angle between normal and projection of normal on x-z plane
  
    term = std::sqrt(nx*nx + nz*nz);
    if (term > 1.0) { term = 1.0; }         // Rounding error fix
    else if (term < -1.0) { term = -1.0; }

    base = std::acos(term);
    if (ny <= 0.0) { a2 = base; }
    else { a2 = -base; }
  }

  // Transformation from inertial frame to local frame

  return euler_rotation(a2*180./M_PI, -a1*180./M_PI, 0.0);
}

/******************************************************************************/
//
// Computes equation of a plane from a point and normal.
// The plane has the equation:
//   ax + by + cz = d
// a, b, c, and d are computed using the relation:
//   normal dot plane_vector = 0, where plane_vector is (x - px, y - py, z - pz)
//
/******************************************************************************/
void compute_plane ( const Eigen::Vector3d & point,
                     const Eigen::Vector3d & norm, double & a, double & b,
                     double & c, double & d )
{
	a = norm(0);
	b = norm(1);
	c = norm(2);
	d = norm(0)*point(0) + norm(1)*point(1) + norm(2)*point(2);
}

/******************************************************************************/
//
// Intersection between a line and a plane. The line is given by a point and
// a vector along the line. The plane is given by a, b, c, d where:
//   ax + by + cz = d is the equation of the plane.
//
// The line is given by:
//   p = p0 + t*l, where p0 is the initial point, l is the vector along the
//   line, and t is the parametric coordinate.

// Note: no checking is done for lines lying on the plane or parallel to the
// plane (i.e., infinite solutions or no solutions). If you try to use it for
// such a case, the result will probably be a divide by 0.
//
/******************************************************************************/
Eigen::Vector3d line_plane_intersection ( const Eigen::Vector3d & p0,
                                    const Eigen::Vector3d & l, const double & a,
                                    const double & b, const double & c,
                                    const double & d )
{
	double den, num, t;
	Eigen::Vector3d p;

	num = d - a*p0(0) - b*p0(1) - c*p0(2);
	den = a*l(0) + b*l(1) + c*l(2);
	t = num/den;
	p = p0 + t*l;

	return p;
}
