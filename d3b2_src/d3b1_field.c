/*
 * d3b1_field.c
 *
 *  Created on: Jan 21, 2019
 *      Author: ohta
 */

#include "bem3_emf_b2.h"



int mEMP_s(double complex *U,double *rt,int type,DOMD *md)
{
  double complex CC[9],kc;
  double F;
  int did,s,sd,l,n;

  did=domain_id(rt,md);

  for(l=0;l<4;l++) U[l]=0.0;
  F=0.0;
  kc=md->kn[did];

  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];
    coef_rt(CC,rt,sd,kc,type,&(md->bd));

    for(n=0;n<4;n++)
      for(l=0;l<4;l++) U[l]+=CC[n+0]*md->bd.sb[did].dU[s][n][l]-CC[n+4]*md->bd.sb[did].U[s][n][l];

    F+=creal(CC[8]);
  }

  if(did==0) for(l=0;l<4;l++) U[l]/=1.0+F;
  else for(l=0;l<4;l++) U[l]/=F;

  return did;
}

int mEMP_t(double complex *U,double *rt,int type,DOMD *md)
{
  double complex e[3],h[3];
  int did,l;

  did=mEMP_s(U,rt,type,md);
  if(did==0){
    calc_mfb_EH(e,h,rt,&(md->mw));
    for(l=0;l<3;l++) U[l]+=e[l];
  }

  return did;
}

int mEMP_i(double complex *U,double *rt,int type,DOMD *md)
{
  double complex e[3],h[3];
  int did,i;

  did=domain_id(rt,md);
  calc_mfb_EH(e,h,rt,&(md->mw));

  for(i=0;i<3;i++) U[i]=e[i];
  U[3]=0.0;

  return did;
}

int EH_mEMP_s(double complex *E,double complex *H,double *rt,int type,DOMD *md)
{
  double complex CC[9],dCC[9],kc,U[4],dU[3][4],ch;
  double v[3][3],F,dF[3];
  int did,l,sd,s,i,n;

  did=domain_id(rt,md);
  kc=md->kn[did];
  ch=md->mw.lambda_0/(2.0*M_PI*I);

  F=0.0;
  for(i=0;i<4;i++) U[i]=0.0;
  for(i=0;i<3;i++){
    for(l=0;l<4;l++) dU[i][l]=0.0;

    for(l=0;l<3;l++) {
      if(i==l) v[i][l]=1.0;
      else v[i][l]=0.0;
    }
    dF[i]=0.0;
  }

  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];

    for(i=0;i<3;i++){
      dcoef_rt(CC,dCC,rt,&(v[i][0]),sd,kc,type,&(md->bd));

      for(n=0;n<4;n++)
        for(l=0;l<4;l++) dU[i][l]+=dCC[n+0]*md->bd.sb[did].dU[s][n][l]-dCC[n+4]*md->bd.sb[did].U[s][n][l];
      dF[i]+=creal(dCC[8]);
    }

    for(n=0;n<4;n++)
      for(l=0;l<4;l++) U[l]+=CC[n+0]*md->bd.sb[did].dU[s][n][l]-CC[n+4]*md->bd.sb[did].U[s][n][l];
    F+=creal(CC[8]);
  }

  if(did==0){
    for(l=0;l<4;l++){
      U[l]/=1.0+F;
      for(i=0;i<3;i++){
        dU[i][l]=(dU[i][l]-U[l]*dF[i])/(1.0+F);
      }
    }
  }
  else {
    for(l=0;l<4;l++){
      U[l]/=F;
      for(i=0;i<3;i++){
        dU[i][l]=(dU[i][l]-U[l]*dF[i])/F;
      }
    }
  }

  if(fabs(dF[0])>CBD_CDF || fabs(dF[1])>CBD_CDF || fabs(dF[2])>CBD_CDF){
    for(i=0;i<3;i++){
      E[i]=0.0;
      H[i]=0.0;
    }
    return -did;
  }

  for(i=0;i<3;i++){
    E[i]=-dU[i][3]+U[i];
  }
  H[0]=ch*(dU[1][2]-dU[2][1]);
  H[1]=ch*(dU[2][0]-dU[0][2]);
  H[2]=ch*(dU[0][1]-dU[1][0]);

  return did;
}

int EH_mEMP_t(double complex *E,double complex *H,double *rt,int type,DOMD *md)
{
  double complex e[3],h[3];
  int did,l;

  did=EH_mEMP_s(E,H,rt,type,md);
  if(did<0) return did;

  if(did==0){
    calc_mfb_EH(e,h,rt,&(md->mw));
    for(l=0;l<3;l++){
      E[l]+=e[l];
      H[l]+=h[l];
    }
  }

  return did;
}

int EH_mEMP_i(double complex *E,double complex *H,double *rt,int type,DOMD *md)
{
  int did;

  did=domain_id(rt,md);
  calc_mfb_EH(E,H,rt,&(md->mw));

  return did;
}

int EH_s(double complex *E,double complex *H,double *rt,int type,DOMD *md)
{
  double complex CC[9];
  double F;
  int did,s,sd,l,n;

  did=domain_id(rt,md);

  for(l=0;l<3;l++){
    E[l]=0.0;
    H[l]=0.0;
  }
  F=0.0;
  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];
    coef_rt(CC,rt,sd,md->kn[did],type,&(md->bd));

    for(n=0;n<4;n++){
      for(l=0;l<3;l++){
        E[l]+=CC[n+0]*md->bd.sb[did].dE[s][n][l]-CC[n+4]*md->bd.sb[did].E[s][n][l];
        H[l]+=CC[n+0]*md->bd.sb[did].dH[s][n][l]-CC[n+4]*md->bd.sb[did].H[s][n][l];
      }
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    for(l=0;l<3;l++){
      E[l]/=1.0+F;
      H[l]/=1.0+F;
    }
  }
  else{
    for(l=0;l<3;l++){
      E[l]/=F;
      H[l]/=F;
    }
  }
  return did;
}

int EH_t(double complex *E,double complex *H,double *rt,int type,DOMD *md)
{
  double complex e[3],h[3];
  int did,l;

  did=EH_s(E,H,rt,type,md);
  if(did==0){
    calc_mfb_EH(e,h,rt,&(md->mw));
    for(l=0;l<3;l++){
      E[l]+=e[l];
      H[l]+=h[l];
    }
  }

  return did;
}

int EH_i(double complex *E,double complex *H,double *rt,int type,DOMD *md)
{
  int did;

  did=domain_id(rt,md);
  calc_mfb_EH(E,H,rt,&(md->mw));

  return did;
}

void EH_s_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md)
{
  double complex CC[9];
  double F,rt[3];
  int s,sd,l,n,td;

  for(l=0;l<3;l++){
    E[l]=0.0;
    H[l]=0.0;
  }
  F=0.0;

  td=md->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(md->bd));

  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];

    coef_bd(CC,rt,td,zeta_t,eta_t,sd,md->kn[did],type,&(md->bd));

    for(n=0;n<4;n++){
      for(l=0;l<3;l++){
        E[l]+=CC[n+0]*md->bd.sb[did].dE[s][n][l]-CC[n+4]*md->bd.sb[did].E[s][n][l];
        H[l]+=CC[n+0]*md->bd.sb[did].dH[s][n][l]-CC[n+4]*md->bd.sb[did].H[s][n][l];
      }
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    for(l=0;l<3;l++){
      E[l]/=1.0+F;
      H[l]/=1.0+F;
    }
  }
  else{
    for(l=0;l<3;l++){
      E[l]/=F;
      H[l]/=F;
    }
  }
}

void EH_t_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md)
{
  double complex CC[9],e[3],h[3];
  double F,rt[3];
  int s,sd,l,n,td;

  for(l=0;l<3;l++){
    E[l]=0.0;
    H[l]=0.0;
  }
  F=0.0;

  td=md->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(md->bd));

  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];

    coef_bd(CC,rt,td,zeta_t,eta_t,sd,md->kn[did],type,&(md->bd));

    for(n=0;n<4;n++){
      for(l=0;l<3;l++){
        E[l]+=CC[n+0]*md->bd.sb[did].dE[s][n][l]-CC[n+4]*md->bd.sb[did].E[s][n][l];
        H[l]+=CC[n+0]*md->bd.sb[did].dH[s][n][l]-CC[n+4]*md->bd.sb[did].H[s][n][l];
      }
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    calc_mfb_EH(e,h,rt,&(md->mw));
    for(l=0;l<3;l++){
      E[l]=E[l]/(1.0+F)+e[l];
      H[l]=H[l]/(1.0+F)+h[l];
    }
  }
  else{
    for(l=0;l<3;l++){
      E[l]/=F;
      H[l]/=F;
    }
  }
}

void EH_i_bd(double complex *E,double complex *H,int did,int t,double zeta_t,double eta_t,int type,DOMD *md)
{
  int td;
  double rt[3];
  int i;

  td=md->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(md->bd));

  if(did==0){
    calc_mfb_EH(E,H,rt,&(md->mw));
  }
  else {
    for(i=0;i<3;i++){
      E[i]=0.0;
      H[i]=0.0;
    }
  }
}

int dmEMP_dt_node(double complex *U,double complex *dUdtz,double complex *dUdte,double *vtz,double *vte,
    int did,int t,int nid,int type,DOMD *md)
{
  double complex CC[9],dCz[9],dCe[9];
  double zt,et,rt[3],F,dFz,dFe;
  size_t i,l,s,n;
  int sd,td,atd;

  td=md->bd.sb[did].sid[t];
  atd=abs(td);

  if(nid>3) return -1;
  if(ELT4==check_element_type(td,&(md->bd))){
    zt=md->bd.zt_44[nid];
    et=md->bd.et_44[nid];
  }
  else {
    if(nid>2) return -1;
    zt=md->bd.zt_34[nid];
    et=md->bd.et_34[nid];
  }

  for(i=0;i<3;i++) rt[i]=md->bd.ren[atd][nid][i];
  tz_te_bd_node(vtz,vte,td,nid,&(md->bd));

  for(l=0;l<4;l++){
    U[l]=0.0;
    dUdtz[l]=0.0;
    dUdte[l]=0.0;
  }
  F=0.0;
  dFz=0.0;
  dFe=0.0;

  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];

    dcoef_bd_t2(CC,dCz,dCe,rt,vtz,vte,td,zt,et,sd,md->kn[did],type,&(md->bd));

    F+=creal(CC[8]);
    dFz+=creal(dCz[8]);
    dFe+=creal(dCe[8]);
    for(n=0;n<4;n++){
      for(l=0;l<4;l++){
        U [l]   += CC[n+0]*md->bd.sb[did].dU[s][n][l]- CC[n+4]*md->bd.sb[did].U[s][n][l];
        dUdtz[l]+=dCz[n+0]*md->bd.sb[did].dU[s][n][l]-dCz[n+4]*md->bd.sb[did].U[s][n][l];
        dUdte[l]+=dCe[n+0]*md->bd.sb[did].dU[s][n][l]-dCe[n+4]*md->bd.sb[did].U[s][n][l];
      }
    }
  }
  if(did==0){
    for(l=0;l<4;l++){
      U[l]/=(1.0+F);
      dUdtz[l]=(dUdtz[l]-U[l]*dFz)/(1.0+F);
      dUdte[l]=(dUdte[l]-U[l]*dFe)/(1.0+F);
    }
  }
  else {
    for(l=0;l<4;l++){
      U[l]/=F;
      dUdtz[l]=(dUdtz[l]-U[l]*dFz)/F;
      dUdte[l]=(dUdte[l]-U[l]*dFe)/F;
    }
  }

  return 0;
}

