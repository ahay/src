#ifndef __SEAM_DEFAULTS__
#define __SEAM_DEFAULTS__

/** Default values for physical constants and units */

#define FPEAK_DEF 0.01                /* default peak frequency = 10 Hz       */
#define CMAX_DEF  5.00                /* default max velocity   = 5 m/ms      */
#define CMIN_DEF  1.00                /* default min velocity   = 1 m/ms      */
#define DMAX_DEF  3.00                /* default max density    = 3 g/cm^3    */
#define DMIN_DEF  5.00                /* default min density    = 0.5 g/cm^3  */
#define CFL_DEF   0.80                /* default cfl fraction   = 0.8         */
#define CREF_DEF  1.50                /* default ref vel        = 1.5 m/ms    */
#define RREF_DEF  1000.0              /* default ref dist       = 1000 m      */ 
#define RHO_DEF   1.0                 /* default ref density    = 1 g/cm^3    */
#define KAPPA_DEF 2.25                /* default bulk modulus   = 2.25 GPa    */

#define KMAX_DEF  CMAX_DEF*CMAX_DEF*DMAX_DEF
                                      /* default max bulk modulus             */
#define KMIN_DEF  CMIN_DEF*CMIN_DEF*DMIN_DEF
                                      /* default min bulk modulus             */
#define BMAX_DEF  1.0/DMIN_DEF
                                      /* default max buoyancy                 */
#define BMIN_DEF  1.0/DMAX_DEF        
                                      /* default min buoyancy                 */

/** m, ms, g (and g/cm^3 for density) - this is one unit 
    system in which all acoustic params are O(1) */
#define KAPPA_MMSGCM3 2.25
#define RHO_MMSGCM3   1.00
/** default unit choice: m, ms, kg */
#define SEAM_UNITS    "MMSGCM3"

#endif
