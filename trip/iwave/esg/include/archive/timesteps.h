/*
Time step functions.
*/
int sgn_ts1d_22(RDOM *dom, int iarr, int it, void *pars);
int sgn_ts1d_24(RDOM *dom, int iarr, int it, void *pars);
int sgn_ts1d_210(RDOM *dom, int iarr, int it, void *pars);
/*int sgn_ts1d_2K(RDOM *dom, int iarr, int it, void *pars);*/

int sgn_ts2d_22(RDOM *dom, int iarr, int it, void *pars);
int sgn_ts2d_24(RDOM *dom, int iarr, int it, void *pars);

/* NOTE: calling this function returns error code 01.11.09 */
int sgn_ts2d_210(RDOM *dom, int iarr, int it, void *pars);

int sgn_ts2d_2K(RDOM *dom, int iarr, int it, void *pars);

int sgn_ts3d_22(RDOM *dom, int iarr, int it, void *pars);
int sgn_ts3d_24(RDOM *dom, int iarr, int it, void *pars);
int sgn_ts3d_210(RDOM *dom, int iarr, int it, void *pars);
int sgn_ts3d_2K(RDOM *dom, int iarr, int it, void *pars);

/*
Model timestep function and posttimestep function
*/
int sgn_modelts(RDOM *dom, int iarr, int it, void *pars);
/*D.S. 06.08.09 : add in two arguments (RDOM *dom, *rdom) and cut off one (int it) 
  WWS 19.06.09: rename, restore original */
int sgn_modelgts(RDOM *dom, RDOM *rdom, RDOM *cdom, int iarr, void *pars); 
int sgn_modelpostts(RDOM *dom, int iarr, int it, void *pars);
/*----------------------------------------------------------------------------*/
/*
Auxilary functions.
Center pressure functions call corresponding sg pressure function
if single pressure array.
*/
int sgn_ts2d_22p(RDOM *dom, void *pars);              /* p0,p1: no etas       */
int sgn_ts2d_22p0(RDOM *dom, void *pars);             /* p0,p1: etax          */
int sgn_ts2d_22p1(RDOM *dom, void *pars);             /* p0,p1: etay          */
int sgn_ts2d_22p01(RDOM *dom, void *pars);            /* p0,p1: etax, etay    */
int sgn_ts2d_22v0(RDOM *dom, void *pars);             /* v0: etax             */
int sgn_ts2d_22v1(RDOM *dom, void *pars);             /* v1: etay             */

int sgn_ts3d_22p(RDOM *dom, void *pars);              /* p0,p1,p2: no etas    */
int sgn_ts3d_22p0(RDOM *dom, void *pars);             /* p0,p1,p2: etax       */
int sgn_ts3d_22pj(RDOM *dom, void *pars, int ind);    /* p0,p1,p2: etaj       */
int sgn_ts3d_22p0j(RDOM *dom, void *pars, int ind);   /* p0,p1,p2: etax, etaj */
int sgn_ts3d_22p12(RDOM *dom, void *pars);            /* p0,p1,p2: etay, etaz */
int sgn_ts3d_22p012(RDOM *dom, void *pars);           /* p0,p1,p2: all etas   */
int sgn_ts3d_22v0(RDOM *dom, void *pars);             /* v0:       eta0       */
int sgn_ts3d_22vj(RDOM *dom, void *pars, int ind);    /* vj:       etaj       */
/*----------------------------------------*/
/* replaced by generalized fcns 06.09 DS */

int sgn_ts2d_24p(RDOM *dom, void *pars);
int sgn_ts2d_24p0(RDOM *dom, void *pars);
int sgn_ts2d_24p1(RDOM *dom, void *pars);
int sgn_ts2d_24p01(RDOM *dom, void *pars);
int sgn_ts2d_24v0(RDOM *dom, void *pars);
int sgn_ts2d_24v1(RDOM *dom, void *pars);

int sgn_ts3d_24p(RDOM *dom, void *pars);
int sgn_ts3d_24p0(RDOM *dom, void *pars);
int sgn_ts3d_24pj(RDOM *dom, void *pars, int ind);
int sgn_ts3d_24p0j(RDOM *dom, void *pars, int ind);
int sgn_ts3d_24p12(RDOM *dom, void *pars);
int sgn_ts3d_24p012(RDOM *dom, void *pars);
int sgn_ts3d_24v0(RDOM *dom, void *pars);
int sgn_ts3d_24vj(RDOM *dom, void *pars, int ind);
/*----------------------------------------*/

/* NOTE: NOT IMPLEMENTED 01.11.09
int sgn_ts2d_210p(RDOM *dom, void *pars);
int sgn_ts2d_210p0(RDOM *dom, void *pars);
int sgn_ts2d_210p1(RDOM *dom, void *pars);
int sgn_ts2d_210p01(RDOM *dom, void *pars);
int sgn_ts2d_210v0(RDOM *dom, void *pars);
int sgn_ts2d_210v1(RDOM *dom, void *pars);
*/

int sgn_ts3d_210p(RDOM *dom, void *pars);
int sgn_ts3d_210p0(RDOM *dom, void *pars);
int sgn_ts3d_210pj(RDOM *dom, void *pars, int ind);
int sgn_ts3d_210p0j(RDOM *dom, void *pars, int ind);
int sgn_ts3d_210p12(RDOM *dom, void *pars);
int sgn_ts3d_210p012(RDOM *dom, void *pars);
int sgn_ts3d_210v0(RDOM *dom, void *pars);
int sgn_ts3d_210vj(RDOM *dom, void *pars, int ind);
/*----------------------------------------*/

int sgn_ts2d_2Kp(RDOM *dom, void *pars);
int sgn_ts2d_2Kp0(RDOM *dom, void *pars);
int sgn_ts2d_2Kp1(RDOM *dom, void *pars);
int sgn_ts2d_2Kp01(RDOM *dom, void *pars);
int sgn_ts2d_2Kv0(RDOM *dom, void *pars);
int sgn_ts2d_2Kv1(RDOM *dom, void *pars);

int sgn_ts3d_2Kp(RDOM *dom, void *pars);
int sgn_ts3d_2Kp0(RDOM *dom, void *pars);
int sgn_ts3d_2Kpj(RDOM *dom, void *pars, int ind);
int sgn_ts3d_2Kp0j(RDOM *dom, void *pars, int ind);
int sgn_ts3d_2Kp12(RDOM *dom, void *pars);
int sgn_ts3d_2Kp012(RDOM *dom, void *pars);
int sgn_ts3d_2Kv0(RDOM *dom, void *pars);
int sgn_ts3d_2Kvj(RDOM *dom, void *pars, int ind);

/*----------------------------------------------------------------------------*/
/* D.S. 28.05.09 */ 
/*
linpre timestep function

modified to keep original interface 19.06.09 - WWS
- change name to gts instead of ts
- define ts by calling gts
*/
/*int sgn_bornts(RDOM *dom, RDOM *rdom, RDOM *cdom, int iarr, void *pars);*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
Time step functions - linpre timestep
Add in two arguments (RDOM * rdom, * cdom) and cut off one (int it) 
*/
/*
int sgn_gts2d_24(RDOM *dom,  RDOM *rdom, RDOM *cdom, int iarr, void *pars);
*/
/*
Auxilary functions - linpre timestep.
Center pressure functions call corresponding sg pressure function
if single pressure array.
*/
/*----------------------------------------*/
/*
int sgn_gts2d_24p(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars);
int sgn_gts2d_24p0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars);
int sgn_gts2d_24p1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars);
int sgn_gts2d_24p01(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars);
int sgn_gts2d_24v0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars);
int sgn_gts2d_24v1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars);
*/
/*----------------------------------------*/
