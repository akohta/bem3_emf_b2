/*
 * bem3_emf_b2.h
 *
 *  Created on: Feb 16, 2019
 *      Author: ohta
 */

#ifndef BEM3_EMF_B2_H_
#define BEM3_EMF_B2_H_

#include "d3b1_const.h"
#include "d3b1_elem.h"

#define FNLGTH   256 // maximum lengh of filename
#define ITER_MIN  10 // minimum number of iterative operations
#define ITER_MAX 100 // maximum number of iterative operations
#define PREC_DEF   2 // type setting for coefficient matrix and fields for force analysis, 0:4p-GL, 1:9p-GL, 2:GLNp-GL,3:GHNp-GL
#define CCD   1.0e-6 // for convergence criterion 

typedef struct multi_object_data{
  int N;         // number of objects 

  char **mod_fn; // object datafile name
  double **rs;   // translation vector of each object

  DOMD *md;      // object data 
} MOBJ;

// ---------- for create_cmatrix ( create coefficient matrix of single object ) ------------
// setup2.c
void read_domd2(int argc,char **argv,DOMD *md); // read datafile 
void print_domd2(DOMD *md);                     // print readed data 
void print_domd2_mksa(DOMD *md);                // print readed data in MKSA system of units
void initialize_domd2(DOMD *md);                // memory allocation and initialize the data
void finalize_domd2(DOMD *md);                  // free allocated memory
int domain_id(double *rt,DOMD *md);             // return domain id of point rt, return the main domain id if rt is on boundary
// solve_bieq2.c
void create_matrix(DOMD *md,char *ofn);         // create coefficient matrix data ( file extension is .cmat ) and its inverse matrix ( solution of boundary integral equations )


// ----------- for d3b2_bv_solver -------------------------------------
// mo_setup.c
void mo_initialize(int argc,char **argv,MOBJ *mo);   // read datafile and initalize multi object 
void mo_print_data(MOBJ *mo);                        // print readed data
void mo_print_data_mksa(MOBJ *mo);                   // print readed data in MKSA system of units
void mo_finalize(MOBJ *mo);                          // free allocated memory 
void mo_dat_write(char *fname,MOBJ *mo);             // output analysis result to binary file with specified filename
void mo_dat_read(char *fname,MOBJ *mo);              // read datafile outputed by mo_dat_write()
void mo_output_node_particles(char *fname,MOBJ *mo); // outputs the nodes as point cloud data ( .particles file ) 
// mo_solve.c
void mo_solve_bv(MOBJ *mo);                        // solve boundary integral equation by using iterative solution


// ----------- for analysis -------------------------------------------
// mo_field.c
void mo_object_domain_id(int *oid,int *did,double *rt,MOBJ *mo);
// outputs
// oid : object id ( oid=-1:open region ), did : domain id of the object.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z.
// return the main domain id if rt is on the boundary.
// the smaller domain id side is selected as main domain.
 
// analysis modified electromagnetic potential by using boundary integral equations 
int mo_mEMP_s(double complex *U,double *rt,int type,MOBJ *mo); // internal field and total scattered field  
int mo_mEMP_t(double complex *U,double *rt,int type,MOBJ *mo); // total field 
int mo_mEMP_i(double complex *U,double *rt,int type,MOBJ *mo); // initial incident field 
// outputs
// U[0]=Ux, U[1]=Uy, U[2]=Uz ( vector potential ), U[3]=phi ( scalar potential ), return object id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of MOBJ.
// type=0:4-point GL, type=1:9-pont or 7-point GL, type=2:GLN-point GL, type=3:GHN-point GL, type=4:DE integration.
// the higher the type numbers are, the lower error is, the slower calculation speed is. it is usually setted to 0 or 1. 
// the type must be set higher number when analyse near the boundary ( in the order of a element ).

// analysis electromagnetic field by using boundary integral equations
int mo_EH_s(double complex *E,double complex *H,double *rt,int type,MOBJ *mo); // internal and total scattered field
int mo_EH_t(double complex *E,double complex *H,double *rt,int type,MOBJ *mo); // total field 
int mo_EH_i(double complex *E,double complex *H,double *rt,int type,MOBJ *mo); // initial incident field 
// outputs
// E[0]=Ex, E[1]=Ey, E[2]=Ez, H[0]=Hx, H[1]=Hy, H[2]=Hz, return object id.
// inputs
// rt[0]=x, rt[1]=y, rt[2]=z, type, pointer of MOBJ.
// type is the same as mo_mEMP_*()

// analysis electromagnetic field by using derivative boundary integral equations
int mo_EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,MOBJ *mo); // internal and total scattered field
int mo_EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,MOBJ *mo); // total field
int mo_EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,MOBJ *mo); // initial incident field 
// outputs and inputs are the same as mo_EH_*() functions.
// these fields are calculated by modified electromagnetic potential using derivative boundary integral equations.
// these are slower than EH_*() functions. the error is smaller than EH_*() functions in far-field.

// analysis electromagnetic field on the boundary 
void mo_EH_s_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo); // internal and scattered field of the single object.
void mo_EH_t_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo); // total field of the single object
void mo_EH_i_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo); // incident field of the single object ( include scattered field of others )
void mo_EH_a_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo); // total scattered field on boundary (exclude initial incident field)
// outputs
// E[0]=Ex, E[1]=Ey, E[2]=Ez, H[0]=Hx, H[1]=Hy, H[2]=Hz.
// inputs
// oid:object id, did:domain id of the object, t:element id defined in each domain, zeta_t,eta_t:point parameter on element surface ( main domain side coordinate ), type, pointer of MOBJ.
// type is the same as mo_EH_*().

// mo_force.c
int mo_force_FN(double *F,double *N,double *rc,int oid,int type,MOBJ *mo);
// outputs
// F[0]=Fx, F[1]=Fy, F[2]=Fz, N[0]=Nx, N[1]=Ny, N[2]=Nz.
// inputs
// rc[0]=x, rc[1]=y, rc[2]=z ( center of rotation ), oid:object id, type, pointer of MOBJ.
// type=0:4-point GL, type!=0:9-point or 7-point GL.
// return 
// 0 : abnormal termination, -1 : no boundary elements

// mo_absorb.c 
int mo_absorb_P(double *P,int oid,int type,MOBJ *mo);
// output
// P : absourbed energy 
// return 
// -2 : loss of significance (catastrophic cancellation) occurred during the surface integral of Poynting vector.
//      The absorbed energy is even smaller than the returned value.
// others are the same as mo_force_FN().

#endif /* BEM3_EMF_B2_H_ */
