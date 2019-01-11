#include <vector>
#include <string>
#include <Eigen/Core>
#include <cmath>
#include "util.h"
#include "algorithms.h"
#include "transformations.h"
#include "sectional_object.h"
#include "vertex.h"
#include "airfoil.h"
#include "section.h"

/******************************************************************************/
//
// Section class. Defines a wing section.
// 
/******************************************************************************/

/******************************************************************************/
//
// Interpolates BL data from airfoil points to section vertices and sets in
// vertex data. scale is optional and defaults to 1.
//
/******************************************************************************/
void Section::setVertexBLData ( const std::vector<double> & bldata,
                                unsigned int dataidx, double scale )
{
    unsigned int i, j1, j2;
    double interpval;

#ifdef DEBUG
    if ( int(bldata.size()) != _foil.nSmoothed() )
        conditional_stop(1, "Section::setVertexBLData",
                         "Wrong size input vector.");
#endif

    for ( i = 0; i < _nverts; i++ )
    {
        j1 = _foilinterp[i].point1;
        j2 = _foilinterp[i].point2;
        interpval = bldata[j1]*_foilinterp[i].weight1
                  + bldata[j2]*_foilinterp[i].weight2;
        _verts[i].setData(dataidx, interpval*scale);
    }
}

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Section::Section ()
{
    _xle = 0.;
    _zle = 0.;
    _chord = 0.;
    _twist = 0.;
    _roll = 0.;
    _nverts = 0;
    _nwake = 0;
    _verts.resize(0); 
    _uverts.resize(0); 
    _wverts.resize(0);
    _re = 0.;
    _cl2dprev = -1.E+06;
    _cl2dguessprev = -1.E+06;
    _cl2dguessprev = -1.E+06;
    _converged = false;
    _unconverged_count = 0;
    _reinitialized = false;
    _foilinterp.resize(0);
}
Section::~Section () {};

/******************************************************************************/
//
// Set or access position, orientation, and scale
//
/******************************************************************************/
void Section::setGeometry ( const double & xle, const double & y,
                            const double & zle, const double & chord,
                            const double & twist, const double & roll )
{
    _xle = xle;
    setY(y);
    _foil.setY(y);
    _zle = zle;
    _chord = chord;
    _twist = twist;
    _roll = roll;
}

void Section::setRoll ( const double & roll ) { _roll = roll; }
const double & Section::xle () const { return _xle; }
const double & Section::zle () const { return _zle; }
const double & Section::chord () const { return _chord; }
const double & Section::twist () const { return _twist; }
const double & Section::roll () const { return _roll; }

/******************************************************************************/
//
// Access vertex
//
/******************************************************************************/
unsigned int Section::nVerts () const { return _nverts; }
Vertex & Section::vert ( unsigned int idx )
{
#ifdef DEBUG
    if ( idx >= _nverts )
        conditional_stop(1, "Section::vert", "Index out of range.");
#endif

    return _verts[idx];
} 

/******************************************************************************/
//
// Set vertices from spacing distribution
//   nchord: number of points along chord
//   lesprat: spacing at leading edge as fraction of uniform spacing
//   tesprat: spacing at leading edge as fraction of uniform spacing
//
/******************************************************************************/
void Section::setVertices ( unsigned int nchord, const double & lesprat,
                            const double & tesprat )
{
    Airfoil foil;
    double slen, sle, unisp, lesp, tesp;
    double a4top, a5top, a4bot, a5bot;
    double sscale, svs;
    std::vector<double> sv;
    std::vector<double> xf, zf;         // Vertices in foil coordinates
    unsigned int i, j, nsmoothed;
    Eigen::Matrix3d rotation;

#ifdef DEBUG
    if (_foil.nBuffer() == 0)
    {
        conditional_stop(1, "Section::setVertices", "Airfoil not loaded yet.");
    }

    if (! _foil.splined())
    {
        conditional_stop(1, "Section::setVertices",
                         "Must call splineFit first.");
    }
#endif

    // Save a copy of the section's airfoil and remove TE gap if necessary. By
    // working with a copy, the gap is only removed for 3D panel calculations, but
    // it is preserved for Xfoil BL calculations.

    foil = _foil;
    if (foil.teGap() > 1.E-14)
    {
        foil.modifyTEGap(0.0, 0.9);
        foil.smoothPaneling();
        foil.splineFit();
    }

    // Get spacings 
    
    slen = foil.sLen();
    sle = foil.sLE();
    unisp = slen / float(2*nchord-2);   // Uniform spacing
    lesp = unisp * lesprat;             // LE spacing
    tesp = unisp * tesprat;             // TE spacing
    
    // Optimize tanh spacing coefficients to minimize stretching
    
    opt_tanh_spacing(nchord, sle, tesp, lesp, a4top, a5top);
    opt_tanh_spacing(nchord, slen-sle, lesp, tesp, a4bot, a5bot);
    
    // Set spacing vector
    
    _nverts = 2*nchord - 1;
    _verts.resize(_nverts);
    _uverts.resize(_nverts);
    sv.resize(_nverts);
    sv[0] = 0.;
    for ( i = 1; i < nchord; i++ )
    {
        sv[i] = sv[i-1] +
                tanh_spacing(i-1, a4top, a5top, nchord, sle, tesp, lesp);
    }
    for ( i = 1; i < nchord; i++ )
    {
        sv[i+nchord-1] = sv[i+nchord-2] +
                  tanh_spacing(i-1, a4bot, a5bot, nchord, slen-sle, lesp, tesp);
    }

    // Get vertices in foil coordinate system (unit chord, 0 <= x <= 1)
    
    xf.resize(_nverts);
    zf.resize(_nverts);
    for ( i = 0; i < _nverts; i++ )
    {
        foil.splineInterp(sv[i], xf[i], zf[i]);
    }

    // Transform to section coordinates. Also store non-rotated, non-translated
    // version for calculating sectional loads.
    
    rotation = inverse_euler_rotation(_roll, _twist, 0.0, "123");
    for ( i = 0; i < _nverts; i++ )
    {
        _verts[i].setCoordinates(xf[i], 0.0, zf[i]);
        _verts[i].translate(-0.25, 0., 0.);
        _verts[i].rotate(rotation);
        _verts[i].translate(0.25, 0., 0.);
        _verts[i].scale(_chord);
        _verts[i].translate(_xle, _y, _zle);

        _uverts[i].setCoordinates(xf[i], 0.0, zf[i]);
        _uverts[i].scale(_chord);
    }

    // Finds interpolation points on airfoil for section vertices

    _foilinterp.resize(_nverts);
    nsmoothed = _foil.nSmoothed();
    
    // Scale section to airfoil length, leaving a little room for roundoff error
    sscale = _foil.sSmoothed(nsmoothed-1) / (sv[_nverts-1] + 1.E-12);

    for ( i = 0; i < _nverts; i++ )
    {
        _foilinterp[i].weight1 = 0.;
        _foilinterp[i].weight2 = 0.;
        svs = sv[i] * sscale;
        for ( j = 0; j < nsmoothed-1; j++ )
        {
            if ( (svs >= _foil.sSmoothed(j)) && (svs <= _foil.sSmoothed(j+1)) )
            {
                _foilinterp[i].point1 = j;
                _foilinterp[i].point2 = j+1;
                _foilinterp[i].weight2 = (svs - _foil.sSmoothed(j))
                                  / (_foil.sSmoothed(j+1) - _foil.sSmoothed(j));
                _foilinterp[i].weight1 = 1. - _foilinterp[i].weight2;
                break;
            }
        }
#ifdef DEBUG
        if ( (_foilinterp[i].weight1 == 0.) && (_foilinterp[i].weight2 == 0.) )
          print_warning("Section::setVertices",
                        "Could not find airfoil->section interpolants." );
#endif
    }
} 

/******************************************************************************/
//
// Prandtl-Glauert geometric transformation (xinc = x/beta)
//
/******************************************************************************/
void Section::transformPrandtlGlauert ( const double & minf )
{
    double beta, x, y, z;
    unsigned int i;

    beta = std::sqrt(1. - std::pow(minf, 2.));
    for ( i = 0; i < _nverts; i++ )
    {
        x = _verts[i].x();
        y = _verts[i].y();
        z = _verts[i].z();
        _verts[i].setIncompressibleCoordinates(x/beta, y, z);
    }
}

/******************************************************************************/
//
// Access airfoil
//
/******************************************************************************/
Airfoil & Section::airfoil () { return _foil; }

/******************************************************************************/
//
// Computes Reynolds number based on local chord. Also sets is for airfoil.
//
/******************************************************************************/
void Section::computeReynoldsNumber ( const double & rhoinf,
                                     const double & uinf, const double & muinf )
{
    _re = rhoinf * uinf * _chord / muinf;
    _foil.setReynoldsNumber(_re);
} 

const double & Section::reynoldsNumber () const { return _re; }

/******************************************************************************/
//
// Sets Mach number for airfoil
//
/******************************************************************************/
void Section::setMachNumber ( const double & minf )
{
    _foil.setMachNumber(minf); 
}

/******************************************************************************/
//
// Boundary layer calculations with Xfoil
//
/******************************************************************************/
void Section::computeBL ( const Eigen::Vector3d & uinfvec,
                          const double & rhoinf, const double & pinf,
                          const double & alpha, int reinit_freq )
{
    Eigen::Vector3d uinfvec_p;
    double qinfp, uinf, uinfp, cl2d, dcl2d, cl2dguessnew;
    double minf, beta, x, y, z;
    Eigen::Matrix3d inertial2section, section2inertial;
    std::vector<double> bldata;
    std::vector<double> xw, zw, dstarw, uedgew;
    int stat;
    unsigned int i;

    uinf = uinfvec.norm();

    // Sectional lift must be computed as an input to Xfoil. Sectional forces
    // and moments will be recomputed after running Xfoil for the purpose of
    // writing data to the sectional output files.

    computeForceMoment(alpha, uinf, rhoinf, pinf, true);

    /** To get 2D Cl:
        1. Transform uinfvec to section frame -> uinfvec_p
        2. 2D angle of attack is atan(uinfvec_p[2]/uinfvec_p[0])
        3. 2D lift is -_fa*sin(alpha2d) + _fn*cos(alpha2d)
    **/

    inertial2section = euler_rotation(_roll, _twist, 0.0, "123");
    uinfvec_p = inertial2section.transpose() * uinfvec;
    uinfp = uinfvec_p.norm();
    qinfp = 0.5*rhoinf*std::pow(uinfp, 2.);
    cl2d = -_fa*uinfvec_p[2]/uinfp + _fn*uinfvec_p[0]/uinfp;
    cl2d /= qinfp*_chord;

    // Approximation of Cl for next iteration
    // Uses 1st order Taylor series approximation for the nonlinear equation
    //  cl = f(x, cl) about clguess. A 0th-order approximation would result
    //  in clguessnew = f(x, clguess). The 1st-order approximation includes
    //  the first derivative of f and has lower error and better convergence
    //  properties.

    if (_cl2dprev > -1.E+06)
    {
        dcl2d = (cl2d - _cl2dprev) / (_cl2dguess - _cl2dguessprev);
        cl2dguessnew = (cl2d - _cl2dguess*dcl2d) / (1. - dcl2d);
    }
    else
        cl2dguessnew = cl2d;
    _cl2dprev = cl2d;
    _cl2dguessprev = _cl2dguess;
    _cl2dguess = cl2dguessnew;

    // Run xfoil at 2D Cl

    _reinitialized = false;
    if (_foil.runXfoil(cl2dguessnew) != 0)
    {
        _converged = false;
        _unconverged_count += 1;
        if (int(_unconverged_count) == reinit_freq)
        {
            _foil.reinitializeBL();
            _reinitialized = true;
            _unconverged_count = 0;
        }
    }
    else
    {
        _converged = true;
        _unconverged_count = 0;
        _reinitialized = false;
    }

    // Interpolate BL quantities to vertices. These will be overwritten for
    // unconverged sections if interpolation/extrapolation is possible. Perform
    // scaling as needed.

    bldata = _foil.blData("cf", stat);
    setVertexBLData(bldata, 9);
    bldata = _foil.blData("deltastar", stat);
    setVertexBLData(bldata, 10, _chord);
    bldata = _foil.blData("ampl", stat);
    setVertexBLData(bldata, 11);
    bldata = _foil.blData("uedge", stat);
    setVertexBLData(bldata, 12, uinf);
    bldata = _foil.blData("cp2d", stat);
    setVertexBLData(bldata, 13);

    // Get BL wake data

    if (_nwake == 0)
    {
        _nwake = _foil.nWake();
        _wverts.resize(_nwake);
    }
    _foil.wakeCoordinates(_nwake, xw, zw);
    dstarw = _foil.wakeDeltastar(_nwake);
    uedgew = _foil.wakeUedge(_nwake);

    // Set wake vertex positions and scaled data

    section2inertial = inverse_euler_rotation(_roll, _twist, 0.0, "123");
    minf = uinf / std::sqrt(1.4*pinf/rhoinf);
    beta = std::sqrt(1. - std::pow(minf, 2.));
    for ( i = 0; i < _nwake; i++ )
    {
        _wverts[i].setCoordinates(xw[i], 0.0, zw[i]);
        _wverts[i].translate(-0.25, 0., 0.);
        _wverts[i].rotate(section2inertial);
        _wverts[i].translate(0.25, 0., 0.);
        _wverts[i].scale(_chord);
        _wverts[i].translate(_xle, _y, _zle);
        _wverts[i].setData(10, dstarw[i]*_chord);
        _wverts[i].setData(12, uedgew[i]*uinf);

        // Prandtl-Glauert transformation

        x = _wverts[i].x();
        y = _wverts[i].y();
        z = _wverts[i].z();
        _wverts[i].setIncompressibleCoordinates(x/beta, y, z);
    }
}

bool Section::blConverged () const { return _converged; }
bool Section::blReinitialized () const { return _reinitialized; }

/******************************************************************************/
//
// Computes force and moment / unit span
//
/******************************************************************************/
void Section::computeForceMoment ( const double & alpha, const double & uinf,
                                   const double & rhoinf, const double & pinf,
                                   bool viscous )
{
    unsigned int i;
    double nx, nz, tx, tz, Vx, Vz, cenx, cenz, pave, cfave, qinf;
    double fap, fav, fnp, fnv, mp, mv, liftp, liftv, dragp, dragv;
    double dfap, dfav, dfnp, dfnv;
    double uinfp, qinfp;
    Eigen::Vector3d forcep, forcev, momentp, momentv, uinfvec, uinfvec_p, fdrag;
    Eigen::Matrix3d section2inertial;
    
    fap = 0.;
    fav = 0.;
    fnp = 0.;
    fnv = 0.;
    mp = 0.;
    mv = 0.;
    qinf = 0.5*rhoinf*std::pow(uinf,2.);
    for ( i = 1; i < _nverts; i++ )
    {
        // Dimensional edge normal vector and edge center
        
        nx = _uverts[i].z() - _uverts[i-1].z();
        nz = _uverts[i-1].x() - _uverts[i].x();
        cenx = 0.5*(_uverts[i].x() + _uverts[i-1].x());
        cenz = 0.5*(_uverts[i].z() + _uverts[i-1].z());
        
        // Integrate pressure
        
        pave = 0.5*(_verts[i].data(5) + _verts[i-1].data(5)) - pinf;
        dfnp = -pave*nz;
        dfap = -pave*nx;
        fnp += dfnp;
        fap += dfap;
        mp += dfnp * (0.25*_chord - cenx) + dfap*cenz;

        // Viscous contribution

        if (viscous)
        {
            Vx = 0.5*(_verts[i].data(2) + _verts[i-1].data(2));
            Vz = 0.5*(_verts[i].data(4) + _verts[i-1].data(4));

            // Tangent direction is positive in local flow direction

            tx = _uverts[i].x() - _uverts[i-1].x();
            tz = _uverts[i].z() - _uverts[i-1].z();
            if (tx*Vx + tz*Vz < 0.)
            {
                tx *= -1.;
                tz *= -1.;
            }

            cfave = 0.5*(_verts[i].data(9) + _verts[i-1].data(9));
            dfnv = cfave*qinf*tz;
            dfav = cfave*qinf*tx;
            fnv += dfnv;
            fav += dfav;
            mv += dfnv * (0.25*_chord - cenx) + dfav*cenz;
        }
    }
    _fn = fnp + fnv;
    _fa = fap + fav;
    
    // Rotate forces and moments to inertial frame
    
    section2inertial = inverse_euler_rotation(_roll, _twist, 0.0, "123");
    forcep << fap, 0., fnp;
    forcep = section2inertial*forcep;
    forcev << fav, 0., fnv;
    forcev = section2inertial*forcev;
    momentp << 0., mp, 0.;
    momentp = section2inertial*momentp;
    momentv << 0., mv, 0.;
    momentv = section2inertial*momentv;
    
    // Sectional lift, drag, and moment coefficients. Note pressure drag
    // comes from Xfoil solution, because integrated pressure drag is not
    // accurate (Xfoil uses BL solution in far wake to compute drag).
    
    liftp = -forcep(0)*sin(alpha*M_PI/180.) + forcep(2)*cos(alpha*M_PI/180.);
    liftv = -forcev(0)*sin(alpha*M_PI/180.) + forcev(2)*cos(alpha*M_PI/180.);
    dragv =  forcev(0)*cos(alpha*M_PI/180.) + forcev(2)*sin(alpha*M_PI/180.);

    // Pressure drag from Xfoil solution, because integrated pressure drag is
    // not accurate (Xfoil uses BL solution in far wake to compute drag).

    uinfvec << uinf*cos(alpha*M_PI/180.), 0., uinf*sin(alpha*M_PI/180.);
    uinfvec_p = section2inertial * uinfvec;
    uinfp = uinfvec_p.norm();
    qinfp = 0.5*rhoinf*std::pow(uinfp, 2.);
    fdrag << _foil.dragCoefficient()*qinfp*_chord, 0., 0.;  // Section frame
    fdrag = section2inertial * fdrag;                       // Inertial frame
    dragp = fdrag(0)*cos(alpha*M_PI/180.) + fdrag(2)*sin(alpha*M_PI/180.)
          - dragv;

    _clp = liftp/(qinf*_chord);
    _clv = liftv/(qinf*_chord);
    _cdp = dragp/(qinf*_chord);
    _cdv = dragv/(qinf*_chord);
    _cmp = momentp(1)/(qinf*_chord*_chord);
    _cmv = momentv(1)/(qinf*_chord*_chord);
}

/******************************************************************************/
//
// Sectional force and moment coefficients
//
/******************************************************************************/
double Section::liftCoefficient () const { return _clp + _clv; }
const double & Section::pressureLiftCoefficient () const { return _clp; }
const double & Section::viscousLiftCoefficient () const { return _clv; }

double Section::dragCoefficient () const { return _cdp + _cdv; }
const double & Section::pressureDragCoefficient () const { return _cdp; }
const double & Section::viscousDragCoefficient () const { return _cdv; }

double Section::pitchingMomentCoefficient () const { return _cmp + _cmv; }
const double & Section::pressurePitchingMomentCoefficient () const
{
    return _cmp;
}
const double & Section::viscousPitchingMomentCoefficient () const
{
    return _cmv;
}

/******************************************************************************/
//
// Viscous wake vertices, scaled and transformed to inertial frame
//
/******************************************************************************/
unsigned int Section::nWake () const { return _nwake; }
Vertex & Section::wakeVert ( unsigned int idx )
{
#ifdef DEBUG
    if (idx >= _nwake)
        conditional_stop(1, "Section::WakeVert", "Index out of range.");
#endif

    return _wverts[idx];
}

/******************************************************************************/
//
// Gets BL data from two other sections. This is needed when Xfoil fails to
// converge.
//
/******************************************************************************/
void Section::interpolateBL ( Section & sec1, Section & sec2,
                              const double & weight1, const double & weight2 )
{
    unsigned int i;
    double cf, dstar, ampl, uedge, cp2d;

#ifdef DEBUG
    if ( (_nverts != sec1.nVerts()) || (_nverts != sec2.nVerts()) )
        conditional_stop(1, "Section::interpolateBL",
                         "Mismatched number of vertices.");
    
    if ( (_nwake != sec1.nWake()) || (_nwake != sec2.nWake()) )
        conditional_stop(1, "Section::interpolateBL",
                         "Mismatched number of wake vertices.");
#endif

    for ( i = 0; i < _nverts; i++ )
    {
        cf = sec1.vert(i).data(9)*weight1 + sec2.vert(i).data(9)*weight2;
        dstar = sec1.vert(i).data(10)*weight1 + sec2.vert(i).data(10)*weight2;
        ampl = sec1.vert(i).data(11)*weight1 + sec2.vert(i).data(11)*weight2;
        uedge = sec1.vert(i).data(12)*weight1 + sec2.vert(i).data(12)*weight2;
        cp2d = sec1.vert(i).data(13)*weight1 + sec2.vert(i).data(13)*weight2;
        _verts[i].setData(9, cf);
        _verts[i].setData(10, dstar);
        _verts[i].setData(11, ampl);
        _verts[i].setData(12, uedge);
        _verts[i].setData(13, cp2d);
    }

    for ( i = 0; i < _nwake; i++ )
    {
        dstar = sec1.wakeVert(i).data(10)*weight1
              + sec2.wakeVert(i).data(10)*weight2;
        uedge = sec1.wakeVert(i).data(12)*weight1
              + sec2.wakeVert(i).data(12)*weight2;
        _wverts[i].setData(10, dstar);
        _wverts[i].setData(12, uedge);
    }
}
