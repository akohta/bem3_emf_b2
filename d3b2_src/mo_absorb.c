/*
 * mo_force.c
 *
 *  Created on: Nov 22, 2021
 *      Author: ohta
 */
#include "bem3_emf_b2.h"

int mo_absorb_P(double *P,int oid,int type,MOBJ *mo)
{
  void mo_bil_absorb_4p(double *P,int oid,int s,MOBJ *mo);
  void mo_bil_absorb_9p(double *P,int oid,int s,MOBJ *mo);
  void mo_lit_absorb_4p(double *P,int oid,int s,MOBJ *mo);
  void mo_lit_absorb_7p(double *P,int oid,int s,MOBJ *mo);
  
  double tp,p;
  int t,td,nc;

  p=0.0;

  nc=0;
  #pragma omp parallel for schedule(dynamic) reduction(+:p,nc) private(td,tp)
  for(t=1;t<=mo->md[oid].bd.sb[0].Ne;t++){
    td=mo->md[oid].bd.sb[0].sid[t];
    nc++;
    if( ELT4==check_element_type(td,&(mo->md[oid].bd)) ){
      if(type==0) mo_bil_absorb_4p(&tp,oid,t,mo);
      else mo_bil_absorb_9p(&tp,oid,t,mo);
    }
    else {
      if(type==0) mo_lit_absorb_4p(&tp,oid,t,mo);
      else mo_lit_absorb_7p(&tp,oid,t,mo);
    }
    p+=tp;
  }
  *P=-p;
  if(*P<0.0){
    *P*=-1.0;
    return -2;
  }
  
  if(nc!=0) return 0;
  else return -1;
}

///////////////////////////////////////////////////////////////////////
void mo_bil_absorb_4p(double *P,int oid,int s,MOBJ *mo)
{
  double complex te[3],th[3],e[3],h[3];
  double vp[3],w[3];
  int i,j,asd;

  asd=abs(mo->md[oid].bd.sb[0].sid[s]);
  *P=0.0;

  for(i=0;i<4;i++){
    for(j=0;j<3;j++){
      e[j]=mo->md[oid].bd.sb[0].E[s][i][j];
      h[j]=mo->md[oid].bd.sb[0].H[s][i][j];
      w[j]=mo->md[oid].bd.wen[asd][i][j];
    }
    mo_EH_i_bd(te,th,oid,0,s,mo->md[oid].bd.zt_44[i],mo->md[oid].bd.et_44[i],PREC_DEF,mo);
    for(j=0;j<3;j++){
      e[j]+=te[j];
      h[j]+=th[j];
    }
    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*mo->md[oid].bd.wt_44[i];
  }

  *P*=-0.5;
}

void mo_bil_absorb_9p(double *P,int oid,int s,MOBJ *mo)
{
  double complex e[3],h[3];
  double vp[3],r[3],w[3],cr[3][4],cw[3][3];
  int i,asd;

  asd=abs(mo->md[oid].bd.sb[0].sid[s]);
  bil_copy_elem_const_rw(cr,cw,asd,&(mo->md[oid].bd));

  *P=0.0;

  for(i=0;i<9;i++){
    bil_rw_zeta_eta(r,w,mo->md[oid].bd.zt_49[i],mo->md[oid].bd.et_49[i],cr,cw);
    mo_EH_t_bd(e,h,oid,0,s,mo->md[oid].bd.zt_49[i],mo->md[oid].bd.et_49[i],PREC_DEF,mo);

    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*mo->md[oid].bd.wt_49[i];
  }
  *P*=-0.5;
}

void mo_lit_absorb_4p(double *P,int oid,int s,MOBJ *mo)
{
  double complex te[3],th[3],e[3],h[3];
  double vp[3],r[3],w[3],cr[3][4],cw[3][3];
  int i,j,asd;

  asd=abs(mo->md[oid].bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(mo->md[oid].bd));

  *P=0.0;

  for(i=0;i<4;i++){
    if(i<3){
      for(j=0;j<3;j++){
        e[j]=mo->md[oid].bd.sb[0].E[s][i][j];
        h[j]=mo->md[oid].bd.sb[0].H[s][i][j];
        r[j]=mo->md[oid].bd.ren[asd][i][j];
        w[j]=mo->md[oid].bd.wen[asd][i][j];
      }
      mo_EH_i_bd(te,th,oid,0,s,mo->md[oid].bd.zt_34[i],mo->md[oid].bd.et_34[i],PREC_DEF,mo);
      for(j=0;j<3;j++){
        e[j]+=te[j];
        h[j]+=th[j];
      }
    }
    else {
      lit_rw_zeta_eta(r,w,mo->md[oid].bd.zt_34[i],mo->md[oid].bd.et_34[i],cr,cw);
      mo_EH_t_bd(e,h,oid,0,s,mo->md[oid].bd.zt_34[i],mo->md[oid].bd.et_34[i],PREC_DEF,mo);
    }

    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*mo->md[oid].bd.wt_34[i];
  }
  *P*=-0.25;
}

void mo_lit_absorb_7p(double *P,int oid,int s,MOBJ *mo)
{
  double complex e[3],h[3];
  double vp[3],r[3],w[3],cr[3][4],cw[3][3];
  int i,asd;

  asd=abs(mo->md[oid].bd.sb[0].sid[s]);
  lit_copy_elem_const_rw(cr,cw,asd,&(mo->md[oid].bd));

  *P=0.0;

  for(i=0;i<7;i++){
    lit_rw_zeta_eta(r,w,mo->md[oid].bd.zt_37[i],mo->md[oid].bd.et_37[i],cr,cw);
    mo_EH_t_bd(e,h,oid,0,s,mo->md[oid].bd.zt_37[i],mo->md[oid].bd.et_37[i],PREC_DEF,mo);

    // poynting vector
    vp[0]=creal(e[1]*conj(h[2])-e[2]*conj(h[1]));
    vp[1]=creal(e[2]*conj(h[0])-e[0]*conj(h[2]));
    vp[2]=creal(e[0]*conj(h[1])-e[1]*conj(h[0]));
    *P+=(vp[0]*w[0]+vp[1]*w[1]+vp[2]*w[2])*mo->md[oid].bd.wt_37[i];
  }
  *P*=-0.25;
}
