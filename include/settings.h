// Stores global settings

#ifndef SETTINGS_H
#define SETTINGS_H

// Xfoil settings

extern double ncrit;
extern double xtript, xtripb;
extern int maxit;
extern double vaccel;
extern bool fix_unconverged;
extern bool reinitialize;

extern int npan;
extern double cvpar;
extern double cterat, ctrrat;

// Algorithm settings

extern double farfield_distance_factor;

#endif
