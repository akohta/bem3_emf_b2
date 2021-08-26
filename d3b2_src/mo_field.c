/*
 * mo_field.c
 *
 *  Created on: Feb 17, 2019
 *      Author: ohta
 */

#include "bem3_emf_b2.h"


void mo_object_domain_id(int *oid,int *did,double *rt,MOBJ *mo)
{
  int o;

  for(o=0;o<mo->N;o++){
    *did=domain_id(rt,&(mo->md[o]));
    if(*did!=0){
      *oid=o;
      return;
    }
  }

  *oid=-1; // opened region
}

int mo_mEMP_s(double complex *U,double *rt,int type,MOBJ *mo)
{
  int mEMP_s(double complex *U,double *rt,int type,DOMD *md); // d3b1_field.c
  
  double complex tu[4];
  int oid,did,m,i;

  mo_object_domain_id(&oid,&did,rt,mo);

  if(oid<0){ // open region
    for(i=0;i<4;i++) U[i]=0.0;
    for(m=0;m<mo->N;m++){
      mEMP_s(tu,rt,type,&(mo->md[m]));
      for(i=0;i<4;i++) U[i]+=tu[i];
    }
  }
  else {
    mEMP_s(U,rt,type,&(mo->md[oid]));
  }

  return oid;
}

int mo_mEMP_t(double complex *U,double *rt,int type,MOBJ *mo)
{
  int mEMP_s(double complex *U,double *rt,int type,DOMD *md); // d3b1_field.c
  
  double complex tu[4],te[3],th[3];
  int oid,did,m,i;

  mo_object_domain_id(&oid,&did,rt,mo);

  if(oid<0){ // open region
    calc_mfb_EH(te,th,rt,&(mo->md[0].mw));
    for(i=0;i<3;i++) U[i]=te[i];
    U[3]=0.0;
    for(m=0;m<mo->N;m++){
      mEMP_s(tu,rt,type,&(mo->md[m]));
      for(i=0;i<4;i++) U[i]+=tu[i];
    }
  }
  else {
    mEMP_s(U,rt,type,&(mo->md[oid]));
  }

  return oid;
}

int mo_mEMP_i(double complex *U,double *rt,int type,MOBJ *mo)
{
  double complex e[3],h[3];
  int oid,did,i;

  mo_object_domain_id(&oid,&did,rt,mo);

  calc_mfb_EH(e,h,rt,&(mo->md[0].mw));

  for(i=0;i<3;i++) U[i]=e[i];
  U[3]=0.0;

  return oid;
}

int mo_EH_s(double complex *E,double complex *H,double *rt,int type,MOBJ *mo)
{
  int EH_s(double complex *E,double complex *H,double *rt,int type,DOMD *md); // d3b1_field.c
  
  double complex te[3],th[3];
  int oid,did,i,m;

  mo_object_domain_id(&oid,&did,rt,mo);

  if(oid<0){ // open region
    for(i=0;i<3;i++){
      E[i]=0.0;
      H[i]=0.0;
    }
    for(m=0;m<mo->N;m++){
      EH_s(te,th,rt,type,&(mo->md[m]));
      for(i=0;i<3;i++){
        E[i]+=te[i];
        H[i]+=th[i];
      }
    }
  }
  else {
    EH_s(E,H,rt,type,&(mo->md[oid]));
  }

  return oid;
}

int mo_EH_t(double complex *E,double complex *H,double *rt,int type,MOBJ *mo)
{
  int EH_s(double complex *E,double complex *H,double *rt,int type,DOMD *md);
  
  double complex te[3],th[3];
  int oid,did,i,m;

  mo_object_domain_id(&oid,&did,rt,mo);

  if(oid<0){ // open region
    calc_mfb_EH(E,H,rt,&(mo->md[0].mw));
    for(m=0;m<mo->N;m++){
      EH_s(te,th,rt,type,&(mo->md[m]));
      for(i=0;i<3;i++){
        E[i]+=te[i];
        H[i]+=th[i];
      }
    }
  }
  else {
    EH_s(E,H,rt,type,&(mo->md[oid]));
  }

  return oid;
}

int mo_EH_i(double complex *E,double complex *H,double *rt,int type,MOBJ *mo)
{
  int oid,did;

  mo_object_domain_id(&oid,&did,rt,mo);
  calc_mfb_EH(E,H,rt,&(mo->md[0].mw));

  return oid;
}

int mo_EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,MOBJ *mo)
{
  int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DOMD *md); // d3b1_field.c
  
  double complex te[3],th[3];
  int oid,did,i,m;

  mo_object_domain_id(&oid,&did,rt,mo);

  if(oid<0){ // open region
    for(i=0;i<3;i++){
      E[i]=0.0;
      H[i]=0.0;
    }
    for(m=0;m<mo->N;m++){
      EH_mEMP_s(te,th,rt,type,&(mo->md[m]));
      for(i=0;i<3;i++){
        E[i]+=te[i];
        H[i]+=th[i];
      }
    }
  }
  else {
    EH_mEMP_s(E,H,rt,type,&(mo->md[oid]));
  }

  return oid; 
}

int mo_EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,MOBJ *mo)
{
  int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DOMD *md);
  
  double complex te[3],th[3];
  int oid,did,i,m;

  mo_object_domain_id(&oid,&did,rt,mo);

  if(oid<0){ // open region
    calc_mfb_EH(E,H,rt,&(mo->md[0].mw));
    for(m=0;m<mo->N;m++){
      EH_mEMP_s(te,th,rt,type,&(mo->md[m]));
      for(i=0;i<3;i++){
        E[i]+=te[i];
        H[i]+=th[i];
      }
    }
  }
  else {
    EH_mEMP_s(E,H,rt,type,&(mo->md[oid]));
  }

  return oid; 
}

int mo_EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,MOBJ *mo)
{
  int oid,did;

  mo_object_domain_id(&oid,&did,rt,mo);
  calc_mfb_EH(E,H,rt,&(mo->md[0].mw));

  return oid;
}


void mo_EH_s_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo)
{
  void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md);
  
  EH_s_bd(E,H,did,t,zeta_t,eta_t,type,&(mo->md[oid]));
}

void mo_EH_t_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo)
{
  void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md); // d3b1_field.c
  int EH_s(double complex *E,double complex *H,double *rt,int type,DOMD *md); // d3b1_field.c
  
  double complex te[3],th[3];
  double rt[3];
  int m,i;

  EH_s_bd(E,H,did,t,zeta_t,eta_t,type,&(mo->md[oid]));

  if(did==0){
    r_bd(rt,mo->md[oid].bd.sb[did].sid[t],zeta_t,eta_t,&(mo->md[oid].bd));
    for(m=0;m<mo->N;m++){
      if(m==oid) continue;
      EH_s(te,th,rt,type,&(mo->md[m]));
      for(i=0;i<3;i++){
        E[i]+=te[i];
        H[i]+=th[i];
      }
    }
    calc_mfb_EH(te,th,rt,&(mo->md[0].mw));
    for(i=0;i<3;i++){
      E[i]+=te[i];
      H[i]+=th[i];
    }
  }
}

void mo_EH_i_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo)
{
  int EH_s(double complex *E,double complex *H,double *rt,int type,DOMD *md); // d3b1_field.c
  
  double complex te[3],th[3];
  double rt[3];
  int i,m;

  if(did==0){
    r_bd(rt,mo->md[oid].bd.sb[did].sid[t],zeta_t,eta_t,&(mo->md[oid].bd));
    calc_mfb_EH(E,H,rt,&(mo->md[0].mw));
    for(m=0;m<mo->N;m++){
      if(m==oid) continue;
      EH_s(te,th,rt,type,&(mo->md[m]));
      for(i=0;i<3;i++){
        E[i]+=te[i];
        H[i]+=th[i];
      }
    }
  }
  else {
    for(i=0;i<3;i++){
      E[i]=0.0;
      H[i]=0.0;
    }
  }
}

void mo_EH_a_bd(double complex *E,double complex *H,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo)
{
  void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md);
  int EH_s(double complex *E,double complex *H,double *rt,int type,DOMD *md);
  
  double complex te[3],th[3];
  double rt[3];
  int m,i;

  EH_s_bd(E,H,did,t,zeta_t,eta_t,type,&(mo->md[oid]));

  if(did==0){
    r_bd(rt,mo->md[oid].bd.sb[did].sid[t],zeta_t,eta_t,&(mo->md[oid].bd));
    for(m=0;m<mo->N;m++){
      if(m==oid) continue;
      EH_s(te,th,rt,type,&(mo->md[m]));
      for(i=0;i<3;i++){
        E[i]+=te[i];
        H[i]+=th[i];
      }
    }
  }
}


