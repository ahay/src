#ifndef __IWAVE_SGNPARS__
#define __IWAVE_SGNPARS__

#include "utils.h"

/*
 * coefficient arrays for finite difference computations
 * not all used (perhaps) - retained just in case
 * [from sgcoeffs.h by I. Terentyev, early versions of IWAVE]
 */

static const ireal COEFF1[] = {                   -1.0e0}; /* 2-2 */
static const ireal COEFF2[] = {             -9.0e0/8.0e0,           1.0e0/24.0e0}; /* 2-4 */
static const ireal COEFF3[] = {           -75.0e0/64.0e0,         25.0e0/384.0e0,           -3.0e0/640.0e0}; /* 2-6 */
static const ireal COEFF4[] = {       -1225.0e0/1024.0e0,       245.0e0/3072.0e0,         -49.0e0/5120.0e0,         5.0e0/7168.0e0}; /* 2-8 */
static const ireal COEFF5[] = {     -19845.0e0/16384.0e0,       735.0e0/8192.0e0,       -567.0e0/40960.0e0,     405.0e0/229376.0e0,      -35.0e0/294912.0e0}; /* 2-10*/
static const ireal COEFF6[] = {   -160083.0e0/131072.0e0,   12705.0e0/131072.0e0,   -22869.0e0/1310720.0e0,   5445.0e0/1835008.0e0,    -847.0e0/2359296.0e0,    63.0e0/2883584.0e0};  /* 2-12*/
static const ireal COEFF7[] = { -1288287.0e0/1048576.0e0, 429429.0e0/4194304.0e0, -429429.0e0/20971520.0e0, 61347.0e0/14680064.0e0, -13013.0e0/18874368.0e0, 3549.0e0/46137344.0e0, -231.0e0/54525952.0e0}; /* 2-14*/

/*----------------------------------------------------------------------------
 * Parameters for time step function.
 *----------------------------------------------------------------------------
 */
typedef struct {
    ireal dt;      /* time step - copied from IMODEL.tsinfo */
    RPNT lam;      /* courant params */
    int k;         /* scheme order */
    int ndim;      /* dimension, copied from IMODEL.grid */
    IPNT lbc;      /* flag left boundary conditions */
    IPNT rbc;      /* flag right boundary conditions */
    ireal * ep0_p;     /* precomputed pml arrays */
    ireal * ep0_pp;    /* p = (1-eta*dt^2) */
    ireal * ev0_p;     /* pp = p/(1+eta*dt^2) */
    ireal * ev0_pp;
    ireal * ep1_p;
    ireal * ep1_pp;   
    ireal * ev1_p;
    ireal * ev1_pp;
    ireal * ep2_p;
    ireal * ep2_pp;   
    ireal * ev2_p;
    ireal * ev2_pp;
    /* slightly redundant dim info for aux arrays */
    int nep0;
    int nev0;
    int nep1;
    int nev1;
    int nep2;
    int nev2;
    /* coefficient arrays for FD schemes - set in readschemeinfo as they
    // are data-dependent
    // encoding (as RPNT): c[diff index][half-order]
    // thus for 3D, c23[0] is the coefficient of u[.][.][.+2]-u[.][.][.-1]  
    // or u[.][.][.+1]-u[.][.][.-2] in a sixth-order difference formula */
    RPNT c11;
    RPNT c12, c22;
    RPNT c14, c24, c34, c44;

} SGN_TS_PARS;  

#endif
