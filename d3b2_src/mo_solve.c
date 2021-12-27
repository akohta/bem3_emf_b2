/*
 * mo_solve.c
 *
 *  Created on: Feb 17, 2019
 *      Author: ohta
 */

#include "bem3_emf_b2.h"

void mo_solve_bv(MOBJ *mo)
{
  void solve_coefficient(int flg,DOMD *md);
  void solve_coefficient_mEMP(int flg,DOMD *md);
  void reset_incident_bd(DOMD *md);
  void add_scattered_mEMP(DOMD *dst,DOMD *src);
  double ccd_f(MOBJ *mo);
  void solve_coefficient_EH(int flg,DOMD *md);
  
  time_t start,end,ms,me;
  double cc;
  int iter,m,d,itermax;

  itermax=ITER_MAX;

  if(mo->N==1){ // single object
    printf("\nsingle object\n");
    solve_coefficient(1,&(mo->md[0]));
  }
  else { // multi object
    printf("\ninitialize boundary value\n");
    for(m=0;m<mo->N;m++){ // calc boundary value each object
      printf("object % 2d\n",m);
      solve_coefficient_mEMP(1,&(mo->md[m]));
    }

    printf("\nthe first iterative operation\n");
    time(&start);
    for(m=0;m<mo->N;m++){ // reset incidient field
      reset_incident_bd(&(mo->md[m]));
    }
    for(m=0;m<mo->N;m++){ // add scattered field to initial incidient field
      printf("object % 2d\n",m);
      printf("  add scattered field       "); fflush(stdout);
      time(&ms);
      for(d=0;d<mo->N;d++) {
        if(d==m) continue;
        add_scattered_mEMP(&(mo->md[d]),&(mo->md[m]));
      }
      time(&me);
      printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
    }
    for(m=0;m<mo->N;m++){
      printf("object % 2d\n",m);
      solve_coefficient_mEMP(1,&(mo->md[m]));
    }
    time(&end);
    printf("Estimated time per one iteration : %5g (sec)\n\n",difftime(end,start));

    printf("iterative operation start. convergence criterion cc < %g and count > %d\n",CCD,ITER_MIN);
    time(&start);
    for(iter=2;iter<itermax;iter++){
      printf("count =% 3d",iter);
      for(m=0;m<mo->N;m++){ // reset incidient field
        reset_incident_bd(&(mo->md[m]));
      }
      for(m=0;m<mo->N;m++){ // add scattered field to incidient field
        for(d=0;d<mo->N;d++) {
          if(d==m) continue;
          add_scattered_mEMP(&(mo->md[d]),&(mo->md[m]));
        }
      }
      for(m=0;m<mo->N;m++){
        solve_coefficient_mEMP(0,&(mo->md[m]));
      }

      cc=ccd_f(mo);
      printf(", cc = %g\n",cc);
      if(cc<CCD && iter>ITER_MIN) break;
    }
    if(iter==itermax){
      printf("The maximum number of iterations has been reached (The result has not converged).\n");
    }
    for(m=0;m<mo->N;m++) solve_coefficient_EH(0,&(mo->md[m]));

    time(&end);
    printf("finished. Elapsed time : %5g (sec)\n",difftime(end,start));
  }
}


///////////////////////////////////////////////////////////////////////
void solve_coefficient(int flg,DOMD *md)
{
  void create_Bmatrix(double complex *B,DOMD *md);
  void solve_ABmatrix(double complex *Bc,DOMD *md);
  void tmatrix_bd_store(double complex *X,DOMD *md); 
  void solve_eh_bv(CMD *cm,DOMD *md); 
  void solve_deh_bv(CMD *cm,DOMD *md); 

  time_t ms,me;
  double complex *B;

  if(flg==1){ printf("  solve mEMP boundary value "); fflush(stdout); }
  time(&ms);
  B=(double complex *)m_alloc2(md->cm.na,sizeof(double complex),"mo_solve.c, solve_coefficient(),B");
  create_Bmatrix(B,md);
  solve_ABmatrix(B,md);
  tmatrix_bd_store(B,md);
  free(B);

  time(&me);
  if(flg==1) printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  if(flg==1){ printf("  solve EH boundary value   "); fflush(stdout); }
  time(&ms);
  solve_eh_bv(&(md->cm),md);
  solve_deh_bv(&(md->cm),md);
  time(&me);
  if(flg==1) printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
}

void solve_coefficient_EH(int flg,DOMD *md)
{
  void solve_eh_bv(CMD *cm,DOMD *md); 
  void solve_deh_bv(CMD *cm,DOMD *md); 

  time_t ms,me;

  if(flg==1){ printf("  solve EH boundary value   "); fflush(stdout); }
  time(&ms);
  solve_eh_bv(&(md->cm),md);
  solve_deh_bv(&(md->cm),md);
  time(&me);
  if(flg==1) printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
}

void create_Bmatrix(double complex *B,DOMD *md)
{
  void create_Bmatrix_dac(size_t *j,size_t did,size_t cid,double complex *B,DOMD *md);
  
  size_t did,cid,j;

  j=0;

  for(did=0;did<=md->MN;did++)
    for(cid=0;cid<4;cid++) create_Bmatrix_dac(&j,did,cid,B,md);

}

void create_Bmatrix_dac(size_t *j,size_t did,size_t cid,double complex *B,DOMD *md)
{
  FILE *fg,*fh;

  double complex *tG,*tH,k1,k2,delt,beps,sigm;
  double vn[3];
  size_t Ne,t,tn,s,asd,mdid,sdid,etype,l;
  int td,sd;

  if((fg=fopen(md->cm.tgfn[did],"rb"))==NULL){    printf("mo_solve.c, create_Bmatrix_dac(),*fg. Failed to open %s file.\n",md->cm.tgfn[did]);    exit(1);  }
  if((fh=fopen(md->cm.thfn[did],"rb"))==NULL){    printf("mo_solve.c, create_Bmatrix_dac(),*fh. Failed to open %s file.\n",md->cm.thfn[did]);    exit(1);  }

  Ne=(size_t)md->bd.sb[did].Ne;
  tG=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"mo_solve.c, create_Bmatrix_dac(),tG"); // malloc
  tH=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"mo_solve.c, create_Bmatrix_dac(),tH"); // malloc

  for(t=1;t<=Ne;t++){
    td=md->bd.sb[did].sid[t];

    for(tn=0;tn<4;tn++){
      if(fread(tG,sizeof(double complex),Ne*4,fg)!=Ne*4){
        printf("mo_solve.c, create_Bmatrix_dac(), failed to read the tG. exit...\n");
        exit(1);
      }
      if(fread(tH,sizeof(double complex),Ne*4,fh)!=Ne*4){
        printf("mo_solve.c, create_Bmatrix_dac(), failed to read the tH. exit...\n");
        exit(1);
      }
      if( tn==3 && ELT3==check_element_type(td,&(md->bd)) )  continue;

      B[*j]=0.0;
      for(s=1;s<=Ne;s++){
        sd=md->bd.sb[did].sid[s]; // signed element id
        asd=(size_t)abs(sd);
        mdid=(size_t)md->bd.md[asd]; // main domain id
        sdid=(size_t)md->bd.sd[asd]; // sub domain id
        etype=check_element_type(sd,&(md->bd));        

        if(did!=mdid){
          k1=md->kn[mdid];
          k2=md->kn[sdid];
          delt=k1*k1-k2*k2;
          beps=k1*k1/(k2*k2);
          sigm=1.0-beps;

          if(cid<3){ // U
            if(etype==ELT3){ // linear-triangular element
              for(l=0;l<3;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                if(mdid==0) B[*j]+=-tH[(s-1)*4+l]*md->bd. Ui[asd][l][cid]
                                   -tG[(s-1)*4+l]*md->bd.dUi[asd][l][cid]
                                   -tG[(s-1)*4+l]*md->bd. Ui[asd][l][  3]*delt*vn[cid];
              }
            }
            else { // bi-linear element
              for(l=0;l<4;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                if(mdid==0) B[*j]+=-tH[(s-1)*4+l]*md->bd. Ui[asd][l][cid]
                                   -tG[(s-1)*4+l]*md->bd.dUi[asd][l][cid]
                                   -tG[(s-1)*4+l]*md->bd. Ui[asd][l][  3]*delt*vn[cid];
              }
            }
          } // end U
          else { // varphi
            if(etype==ELT3){ // linear-triangular element
              for(l=0;l<3;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                if(mdid==0) B[*j]+=-tH[(s-1)*4+l]*md->bd.Ui[asd][l][3]
                                   -sigm*tG[(s-1)*4+l]*(vn[0]*md->bd.Ui[asd][l][0]+vn[1]*md->bd.Ui[asd][l][1]+vn[2]*md->bd.Ui[asd][l][2])
                                   -beps*tG[(s-1)*4+l]*md->bd.dUi[asd][l][3];
              }
            }
            else { // bi-linear element
              for(l=0;l<4;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                if(mdid==0) B[*j]+=-tH[(s-1)*4+l]*md->bd.Ui[asd][l][3]
                                  -sigm*tG[(s-1)*4+l]*(vn[0]*md->bd.Ui[asd][l][0]+vn[1]*md->bd.Ui[asd][l][1]+vn[2]*md->bd.Ui[asd][l][2])
                                  -beps*tG[(s-1)*4+l]*md->bd.dUi[asd][l][3];
              }
            }
          } // end varphi
        } // end sub domain
      } // end for s

      *j+=1;

    } // end for tn

  } // end for t

  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);
}

void solve_ABmatrix(double complex *Bc,DOMD *md)
{
  FILE *fa;
  double complex *tA,*X;
  size_t N,i,j;

  N=md->cm.na;

  tA=(double complex *)m_alloc2(N,sizeof(double complex),"mo_solve.c, solve_ABmatrix(),tA");
  X=(double complex *)m_alloc2(N,sizeof(double complex),"mo_solve.c, solve_ABmatrix(),X");

  if((fa=fopen(md->cm.lupfn,"rb"))==NULL){     printf("mo_solve.c, solve_ABmatrix(), Failed to open the %s file.\n",md->cm.lupfn);    exit(1); }

  for(i=0;i<N;i++){
    if(fread(tA,sizeof(double complex),N,fa)!=N){
      printf("mo_solve.c, solve_ABmatrix(), failed to read the tA. exit...\n");
      exit(1);
    }
    X[i]=0.0;
    for(j=0;j<N;j++) X[i]+=tA[j]*Bc[j];
  }

  for(i=0;i<N;i++){
    Bc[i]=X[i];
  }

  free(tA);
  free(X);
}

void solve_coefficient_mEMP(int flg,DOMD *md)
{
  void tmatrix_bd_store(double complex *X,DOMD *md); 

  time_t ms,me;
  double complex *B;

  if(flg==1){ printf("  solve mEMP boundary value "); fflush(stdout); }
  time(&ms);
  B=(double complex *)m_alloc2(md->cm.na,sizeof(double complex),"mo_solve.c, solve_coefficient_mEMP(),B");
  create_Bmatrix(B,md);

  solve_ABmatrix(B,md);
  tmatrix_bd_store(B,md);
  free(B);

  time(&me);
  if(flg==1) printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
}

void add_scattered_mEMP(DOMD *dst,DOMD *src)
{
  int mEMP2_s(double complex *U,double complex *dU,double *rt,double *v,int type,DOMD *md);
  
  double complex U[4],dU[4];
  double r[3],w[3];
  int i,j,l,did;

  int type=PREC_DEF;

  #pragma omp parallel for schedule(dynamic) private(j,l,r,w,U,dU,did)
  for(i=1;i<=dst->bd.Ne;i++){
    if(dst->bd.ed[i][3]==0){ // linear triangular element
      for(j=0;j<4;j++){
        if(dst->bd.md[i]==0){
          for(l=0;l<3;l++){
            r[l]=dst->bd.ren[i][j][l];
            w[l]=dst->bd.wen[i][j][l];
          }
          vuni_d(w);
          did=mEMP2_s(U,dU,r,w,type,src);
          if(did!=0){
            printf("mo_solve.c, add_scattered_mEMP(), scattered field domain id error!(overlapped). domain id=%d. Exit...\n",did);
            exit(1);
          }
          for(l=0;l<4;l++){
            dst->bd. Ui[i][j][l]+=U[l];
            dst->bd.dUi[i][j][l]+=dU[l];
          }
        }
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        if(dst->bd.md[i]==0){
          for(l=0;l<3;l++){
            r[l]=dst->bd.ren[i][j][l];
            w[l]=dst->bd.wen[i][j][l];
          }
          vuni_d(w);
          did=mEMP2_s(U,dU,r,w,type,src);
          if(did!=0){
            printf("mo_solve.c, add_scattered_mEMP(), scattered field domain id error!(overlapped). domain id=%d. Exit...\n",did);
            exit(1);
          }
          for(l=0;l<4;l++){
            dst->bd. Ui[i][j][l]+=U[l];
            dst->bd.dUi[i][j][l]+=dU[l];
          }
        }
      }
    }
  }
}

int mEMP2_s(double complex *U,double complex *dU,double *rt,double *v,int type,DOMD *md)
{
  double complex CC[9],dCC[9],kc;
  double F,dF;
  int did,s,sd,l,n;

  did=domain_id(rt,md);

  for(l=0;l<4;l++){
    U[l]=0.0;
    dU[l]=0.0;
  }
  F=0.0;
  dF=0.0;
  kc=md->kn[did];

  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];
    dcoef_rt(CC,dCC,rt,v,sd,kc,type,&(md->bd));

    for(n=0;n<4;n++){
      for(l=0;l<4;l++){
        U [l]+= CC[n+0]*md->bd.sb[did].dU[s][n][l]- CC[n+4]*md->bd.sb[did].U[s][n][l];
        dU[l]+=dCC[n+0]*md->bd.sb[did].dU[s][n][l]-dCC[n+4]*md->bd.sb[did].U[s][n][l];
      }
    }

    F+=creal(CC[8]);
    dF+=creal(dCC[8]);
  }

  if(did==0){
    for(l=0;l<4;l++){
      U[l]/=1.0+F;
      dU[l]=(dU[l]-U[l]*dF)/(1.0+F);
    }
  }
  else{
    for(l=0;l<4;l++){
      U[l]/=F;
      dU[l]=(dU[l]-U[l]*dF)/F;
    }
  }

  return did;
}

void reset_incident_bd(DOMD *md)
{
  int i,j,l;

  for(i=1;i<=md->bd.Ne;i++){
    for(j=0;j<4;j++){
      for(l=0;l<4;l++){
        md->bd.Uip[i][j][l]=md->bd. Ui [i][j][l];
        md->bd. Ui[i][j][l]=md->bd. Ui0[i][j][l];
        md->bd.dUi[i][j][l]=md->bd.dUi0[i][j][l];
      }
    }
  }
}

double ccd_f(MOBJ *mo)
{
  double complex *tc,*tp,dv[4];
  double usum,psum,ac,ap,bc,bp,dcp,nrm,db,ds;
  int i,s,p,sc;

  usum=0.0;
  psum=0.0;

  sc=0;

  for(i=0;i<mo->N;i++){
    for(s=1;s<=mo->md[i].bd.Ne;s++){
      if(mo->md[i].bd.md[s]==0){
        for(p=0;p<4;p++){
          tc=mo->md[i].bd.Ui[s][p];
          tp=mo->md[i].bd.Uip[s][p];

          ac=creal(tc[0]*conj(tc[0])+tc[1]*conj(tc[1])+tc[2]*conj(tc[2]));
          ap=creal(tp[0]*conj(tp[0])+tp[1]*conj(tp[1])+tp[2]*conj(tp[2]));
          dv[0]=tc[0]-tp[0];
          dv[1]=tc[1]-tp[1];
          dv[2]=tc[2]-tp[2];
          dv[3]=tc[3]-tp[3];
          dcp=creal(dv[0]*conj(dv[0])+dv[1]*conj(dv[1])+dv[2]*conj(dv[2]));
          nrm=2.0*dcp/(ac+ap);
          usum+=nrm;

          bc=creal(tc[3]*conj(tc[3]));
          bp=creal(tp[3]*conj(tp[3]));
          db=creal(dv[3]*conj(dv[3]));
          ds=2.0*db/(bc+bp);
          psum+=ds;

          sc+=1;
        }
      }
    }
  }

  if(usum>psum) return usum/(double)sc;
  else return psum/(double)sc;
}

void tmatrix_bd_store(double complex *X,DOMD *md)
{
  double complex k1,k2,delt,beps,sigm,tU[4],tUi[4];
  double n[3];
  size_t d,s,l,nn,i;
  int sd,asd,etype,mdid,sdid;

  nn=md->bd.NN;

  for(d=0;d<=md->MN;d++){
    for(s=1;s<=md->bd.sb[d].Ne;s++){
      sd=md->bd.sb[d].sid[s];
      asd=abs(sd);
      mdid=md->bd.md[asd];
      sdid=md->bd.sd[asd];
      etype=check_element_type(sd,&(md->bd));

      if(d==mdid){ // main domain
        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            for(i=0;i<4;i++){ // store U, dU
              md->bd.sb[d]. U[s][l][i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d].dU[s][l][i]=X[nn*(i+4)+md->bd.eni[asd][l]];
            }
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            for(i=0;i<4;i++){ // store U, dU
              md->bd.sb[d]. U[s][l][i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d].dU[s][l][i]=X[nn*(i+4)+md->bd.eni[asd][l]];
            }
          }
        }
      }
      else { // subdomain
        k1=md->kn[mdid];
        k2=md->kn[sdid];
        delt=k1*k1-k2*k2;
        beps=k1*k1/(k2*k2);
        sigm=1.0-beps;

        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            n_bd_node(n,asd,l,&(md->bd));
            tU[3]=X[nn*3+md->bd.eni[asd][l]];
            for(i=0;i<3;i++){ // store U, dU
              tU[i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d]. U[s][l][i]=tU[i];
              md->bd.sb[d].dU[s][l][i]=-(X[nn*(i+4)+md->bd.eni[asd][l]]+delt*tU[3]*n[i]);
              if(mdid==0){
                tUi[i]=md->bd.Ui[asd][l][i];
                md->bd.sb[d]. U[s][l][i]+= tUi[i];
                md->bd.sb[d].dU[s][l][i]+=-md->bd.dUi[asd][l][i];
              }
            }
            // store phi, dphi
            md->bd.sb[d]. U[s][l][3]=tU[3];
            md->bd.sb[d].dU[s][l][3]=-(beps*X[nn*7+md->bd.eni[asd][l]]+sigm*(n[0]*tU[0]+n[1]*tU[1]+n[2]*tU[2]));
            if(mdid==0)  md->bd.sb[d].dU[s][l][3]+=-sigm*(n[0]*tUi[0]+n[1]*tUi[1]+n[2]*tUi[2]);
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            n_bd_node(n,asd,l,&(md->bd));
            tU[3]=X[nn*3+md->bd.eni[asd][l]];
            for(i=0;i<3;i++){ // store U, dU
              tU[i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d]. U[s][l][i]=tU[i];
                md->bd.sb[d].dU[s][l][i]=-(X[nn*(i+4)+md->bd.eni[asd][l]]+delt*tU[3]*n[i]);
              if(mdid==0){
                tUi[i]=md->bd.Ui[asd][l][i];
                md->bd.sb[d]. U[s][l][i]+= tUi[i];
                md->bd.sb[d].dU[s][l][i]+=-md->bd.dUi[asd][l][i];
              }
            }
            // store phi, dphi
            md->bd.sb[d]. U[s][l][3]=tU[3];
            md->bd.sb[d].dU[s][l][3]=-(beps*X[nn*7+md->bd.eni[asd][l]]+sigm*(n[0]*tU[0]+n[1]*tU[1]+n[2]*tU[2]));
            if(mdid==0)  md->bd.sb[d].dU[s][l][3]+=-sigm*(n[0]*tUi[0]+n[1]*tUi[1]+n[2]*tUi[2]);
          }
        }
      } // end subdomain

    }  // end s

  } // end d

}

void solve_eh_bv(CMD *cm,DOMD *md)
{
  void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte);
  void dudt_fcm(double complex *dUz,double complex *dUe,int t,int tn,size_t Ne,double complex ***U,double complex ***dU,double complex *tdG,double complex *tdH,double *tdF);
  
  FILE *fdg,*fdh,*fdf;
  double complex *tdG,*tdH,dUz[4],dUe[4],ch;
  double tdF[3],vtz[3],vte[3],vn[3],a[9],sig;
  size_t Ne,d,t,atd,tn,i;
  int at;

  ch=md->mw.lambda_0/(2.0*M_PI*I);

  for(d=0;d<=md->MN;d++){
    if((fdg=fopen(cm->tdgfn[d],"rb"))==NULL){    printf("mo_solve.c, solve_eh_bv(),*fdg. Failed to open %s file.\n",cm->tdgfn[d]);    exit(1);  }
    if((fdh=fopen(cm->tdhfn[d],"rb"))==NULL){    printf("mo_solve.c, solve_eh_bv(),*fdh. Failed to open %s file.\n",cm->tdhfn[d]);    exit(1);  }
    if((fdf=fopen(cm->tdffn[d],"rb"))==NULL){    printf("mo_solve.c, solve_eh_bv(),*fdf. Failed to open %s file.\n",cm->tdffn[d]);    exit(1);  }

    Ne=md->bd.sb[d].Ne;
    tdG=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"mo_solve.c, solve_eh_bv(),tdG"); // malloc
    tdH=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"mo_solve.c, solve_eh_bv(),tdH"); // malloc

    for(t=1;t<=Ne;t++){
      at=md->bd.sb[d].sid[t];
      if(at<0) sig=-1.0;
      else sig=1.0;
      atd=abs(at);

      for(tn=0;tn<4;tn++){
        if(fread(tdG,sizeof(double complex),2*Ne*4,fdg)!=2*Ne*4){
          printf("mo_solve.c, solve_eh_bv(), failed to read the tdG. exit...\n");
          exit(1);
        }
        if(fread(tdH,sizeof(double complex),2*Ne*4,fdh)!=2*Ne*4){
          printf("mo_solve.c, solve_eh_bv(), failed to read the tdH. exit...\n");
          exit(1);
        }
        if(fread(tdF,sizeof(double),3,fdf)!=3){
          printf("mo_solve.c, solve_eh_bv(), failed to read the tdF. exit...\n");
          exit(1);
        }
        if( tn==3 && ELT3==check_element_type(atd,&(md->bd)) )  continue;

        // tangential vector
        tz_te_bd_node(vtz,vte,atd,tn,&(md->bd));
        // normal vector
        for(i=0;i<3;i++) vn[i]=sig*md->bd.wen[atd][tn][i];
        vuni_d(vn);
        // inv-matrix
        inv_matrix_33(a,vn,vtz,vte);
        // dUdtz,dUdte
        dudt_fcm(dUz,dUe,t,tn,Ne,md->bd.sb[d].U,md->bd.sb[d].dU,tdG,tdH,tdF);
        // electric field
        for(i=0;i<3;i++) md->bd.sb[d].E[t][tn][i]=-(md->bd.sb[d].dU[t][tn][3]*a[i*3+0]+dUz[3]*a[i*3+1]+dUe[3]*a[i*3+2])+md->bd.sb[d].U[t][tn][i];
        // magnetic field
        md->bd.sb[d].H[t][tn][0]=ch*( (md->bd.sb[d].dU[t][tn][2]*a[1*3+0]+dUz[2]*a[1*3+1]+dUe[2]*a[1*3+2])
                                     -(md->bd.sb[d].dU[t][tn][1]*a[2*3+0]+dUz[1]*a[2*3+1]+dUe[1]*a[2*3+2]) );
        md->bd.sb[d].H[t][tn][1]=ch*( (md->bd.sb[d].dU[t][tn][0]*a[2*3+0]+dUz[0]*a[2*3+1]+dUe[0]*a[2*3+2])
                                     -(md->bd.sb[d].dU[t][tn][2]*a[0*3+0]+dUz[2]*a[0*3+1]+dUe[2]*a[0*3+2]) );
        md->bd.sb[d].H[t][tn][2]=ch*( (md->bd.sb[d].dU[t][tn][1]*a[0*3+0]+dUz[1]*a[0*3+1]+dUe[1]*a[0*3+2])
                                     -(md->bd.sb[d].dU[t][tn][0]*a[1*3+0]+dUz[0]*a[1*3+1]+dUe[0]*a[1*3+2]) );
      } // end for tn
    } // end for t
    fclose(fdg);
    fclose(fdh);
    fclose(fdf);
    free(tdG);
    free(tdH);
  } // end for d
}

void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte)
{
  double i_det;
  int i;

  i_det=1.0/(vn[0]*vtz[1]*vte[2]+vn[1]*vtz[2]*vte[0]+vn[2]*vtz[0]*vte[1]
          -( vn[2]*vtz[1]*vte[0]+vn[1]*vtz[0]*vte[2]+vn[0]*vtz[2]*vte[1]));

  ma[0*3+0]= vtz[1]*vte[2]-vtz[2]*vte[1];
  ma[0*3+1]=- vn[1]*vte[2]+ vn[2]*vte[1];
  ma[0*3+2]=  vn[1]*vtz[2]- vn[2]*vtz[1];

  ma[1*3+0]=-vtz[0]*vte[2]+vtz[2]*vte[0];
  ma[1*3+1]=  vn[0]*vte[2]- vn[2]*vte[0];
  ma[1*3+2]=- vn[0]*vtz[2]+ vn[2]*vtz[0];

  ma[2*3+0]= vtz[0]*vte[1]-vtz[1]*vte[0];
  ma[2*3+1]=- vn[0]*vte[1]+ vn[1]*vte[0];
  ma[2*3+2]=  vn[0]*vtz[1]- vn[1]*vtz[0];

  for(i=0;i<9;i++) ma[i]*=i_det;
}

void dudt_fcm(double complex *dUz,double complex *dUe,int t,int tn,size_t Ne,double complex ***U,double complex ***dU,double complex *tdG,double complex *tdH,double *tdF)
{
  size_t i,s,sn;

  for(i=0;i<4;i++){
    dUz[i]=0.0;
    dUe[i]=0.0;
  }
  for(s=1;s<=Ne;s++){
    for(sn=0;sn<4;sn++){
      for(i=0;i<4;i++){
        dUz[i]+=tdG[Ne*4*0+4*(s-1)+sn]*dU[s][sn][i]-tdH[Ne*4*0+4*(s-1)+sn]*U[s][sn][i];
        dUe[i]+=tdG[Ne*4*1+4*(s-1)+sn]*dU[s][sn][i]-tdH[Ne*4*1+4*(s-1)+sn]*U[s][sn][i];
      }
    }
  }
  for(i=0;i<4;i++){
    dUz[i]=(dUz[i]-U[t][tn][i]*tdF[1])/tdF[0];
    dUe[i]=(dUe[i]-U[t][tn][i]*tdF[2])/tdF[0];
  }
}

void solve_deh_bv(CMD *cm,DOMD *md)
{
  FILE *fg,*fh;
  MKL_Complex16 *A,*B;
  double complex *tG,*tH;
  size_t Ne,d,t,atd,tn,nn,nc,nr,s,ns,asd,i;
  MKL_INT *ipiv,nrhs,lda,ldb,info;

  nrhs=6;

  for(d=0;d<=md->MN;d++){
    Ne=md->bd.sb[d].Ne;
    if(Ne==0) continue;

    if((fg=fopen(cm->tgfn[d],"rb"))==NULL){    printf("mo_solve.c, solve_deh_bv(),*fg. Failed to open %s file.\n",cm->tgfn[d]);    exit(1);  }
    if((fh=fopen(cm->thfn[d],"rb"))==NULL){    printf("mo_solve.c, solve_deh_bv(),*fh. Failed to open %s file.\n",cm->thfn[d]);    exit(1);  }

    // matrix size
    nn=0;
    for(s=1;s<=Ne;s++){
      asd=abs(md->bd.sb[d].sid[s]);
      if( ELT3==check_element_type(asd,&(md->bd)) ) nn+=3;
      else nn+=4;
    }
    lda=nn;
    ldb=6;

    A =(MKL_Complex16 *)m_alloc2(nn*nn,sizeof(MKL_Complex16),"mo_solve.c, solve_deh_bv(),A"); // malloc
    B =(MKL_Complex16 *)m_alloc2(nn*6,sizeof(MKL_Complex16),"mo_solve.c, solve_dev_bv(),B"); // malloc
    ipiv=(MKL_INT *)m_alloc2(nn,sizeof(MKL_INT),"mo_solve.c, solve_deh_bv(),ipiv"); // malloc
    tG=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"mo_solve.c, solve_dev_bv(),tG"); // malloc
    tH=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"mo_solve.c, solve_dev_bv(),tH"); // malloc

    nc=0;
    for(t=1;t<=Ne;t++){
      atd=abs(md->bd.sb[d].sid[t]);

      for(tn=0;tn<4;tn++){
        if(fread(tG,sizeof(double complex),Ne*4,fg)!=Ne*4){
          printf("mo_solve.c, solve_deh_bv(), failed to read the tG. exit...\n");
          exit(1);
        }
        if(fread(tH,sizeof(double complex),Ne*4,fh)!=Ne*4){
          printf("mo_solve.c, solve_deh_bv(), failed to read the tH. exit...\n");
          exit(1);
        }
        if( tn==3 && ELT3==check_element_type(atd,&(md->bd)) )  continue;

        for(i=0;i<6;i++){
          B[6*nc+i].real=0.0;
          B[6*nc+i].imag=0.0;
        }
        nr=0;
        for(s=1;s<=Ne;s++){
          asd=abs(md->bd.sb[d].sid[s]);
          if( ELT4==check_element_type(asd,&(md->bd)) ){
            for(ns=0;ns<4;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[6*nc+i+0].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+0].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+3].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
                B[6*nc+i+3].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
              }
            }
          }
          else {
            for(ns=0;ns<3;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[6*nc+i+0].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+0].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+3].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
                B[6*nc+i+3].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
              }
            }
          }
        } // end for s
        nc+=1;
      } // end for tn

    } // end for t

    // solve mkl lapack
    info = LAPACKE_zgesv( LAPACK_ROW_MAJOR, nn, nrhs, A, lda, ipiv, B, ldb );
    if(info!=0){
      printf("mo_solve.c, solve_deh_bv(), LAPACKE_zgesv(). info=%lld error. Exit...\n",info);
      exit(1);
    }

    // store data
    nc=0;
    for(s=1;s<=Ne;s++){
      asd=abs(md->bd.sb[d].sid[s]);
      if( ELT4==check_element_type(asd,&(md->bd)) ){
        for(ns=0;ns<4;ns++){
          for(i=0;i<3;i++) {
            md->bd.sb[d].dE[s][ns][i]=B[6*nc+i+0].real+B[6*nc+i+0].imag*I;
            md->bd.sb[d].dH[s][ns][i]=B[6*nc+i+3].real+B[6*nc+i+3].imag*I;
          }
          nc+=1;
        }
      }
      else {
        for(ns=0;ns<3;ns++){
          for(i=0;i<3;i++) {
            md->bd.sb[d].dE[s][ns][i]=B[6*nc+i+0].real+B[6*nc+i+0].imag*I;
            md->bd.sb[d].dH[s][ns][i]=B[6*nc+i+3].real+B[6*nc+i+3].imag*I;
          }
          nc+=1;
        }
      }
    }

    fclose(fg);
    fclose(fh);
    free(A);
    free(B);
    free(ipiv);
    free(tG);
    free(tH);
  } // end for d
}


