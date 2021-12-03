/*
 * solve_bieq2.c
 *
 *  Created on: Feb 16, 2019
 *      Author: ohta
 */

#include "bem3_emf_b2.h"

void create_matrix(DOMD *md,char *ofn)
{
  void create_cmatrix(CMD *cm,DOMD *md); 
  void initalize_cmd2(CMD *cm);
  void create_Amatrix_csr(DOMD *md);
  void lu_dec_A(CMD *cm);
  void dat_write2(char *fname,DOMD *md);
  void finalize_cmd2(CMD *cm);
  
  time_t start,end,ms,me;

  printf("\ncreate matrix data \n");
  time(&start);

  printf("  coefficient matrix          "); fflush(stdout);
  time(&ms);
  md->cm.type=PREC_DEF; // type setting, 0:4p GL,1:9p or 7p(triangular) GL, 2:GLN p GL, 3: GHN p GL
  md->cm.MN=md->MN;
  initalize_cmd2(&(md->cm));
  create_cmatrix(&(md->cm),md);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  printf("  inverse matrix              "); fflush(stdout);
  time(&ms);
  md->cm.nn=(size_t)md->bd.NN;
  md->cm.na=md->cm.nn*8;
  create_Amatrix_csr(md);
  lu_dec_A(&(md->cm));
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  dat_write2(ofn,md);

  finalize_cmd2(&(md->cm));
  time(&end);
  printf("Total elapsed time : %g (sec)\n",difftime(end,start));
}





//////////////////////////////////////////////////////////////////////////////////////////////////
void initalize_cmd2(CMD *cm)
{
  int i,MN;

  MN=cm->MN;
  cm->tgfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq2.c, initialize_cmd2(), cm->tgfn");
  cm->thfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq2.c, initialize_cmd2(), cm->thfn");
  cm->tdgfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq2.c, initialize_cmd2(), cm->tdgfn");
  cm->tdhfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq2.c, initialize_cmd2(), cm->tdhfn");
  cm->tdffn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq2.c, initialize_cmd2(), cm->tdffn");
  for(i=0;i<=MN;i++){
    cm->tgfn[i]=(char *)m_alloc2(128,sizeof(char ),"solve_bieq2.c, initialize_cmd2(), cm->tgfn[i]");
    cm->thfn[i]=(char *)m_alloc2(128,sizeof(char ),"solve_bieq2.c, initialize_cmd2(), cm->thfn[i]");
    sprintf(cm->tgfn[i],"tmpG_%05d.cmat",i);
    sprintf(cm->thfn[i],"tmpH_%05d.cmat",i);
    cm->tdgfn[i]=(char *)m_alloc2(128,sizeof(char ),"solve_bieq2.c, initialize_cmd2(), cm->tdgfn[i]");
    cm->tdhfn[i]=(char *)m_alloc2(128,sizeof(char ),"solve_bieq2.c, initialize_cmd2(), cm->tdhfn[i]");
    cm->tdffn[i]=(char *)m_alloc2(128,sizeof(char ),"solve_bieq2.c, initialize_cmd2(), cm->tdffn[i]");
    sprintf(cm->tdgfn[i],"tmpdG_%05d.cmat",i);
    sprintf(cm->tdhfn[i],"tmpdH_%05d.cmat",i);
    sprintf(cm->tdffn[i],"tmpdF_%05d.cmat",i);
  }
  cm->lupfn=(char *)m_alloc2(128,sizeof(char ),"solve_bieq2.c, initialize_cmd2(), cm->tgfn[i]");
  sprintf(cm->lupfn,"tmpinvA.cmat");

  cm->aval=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq2.c, initialize_cmd2(), cm->aval");
  cm->aptr=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq2.c, initialize_cmd2(), cm->aptr");
  cm->aidx=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq2.c, initialize_cmd2(), cm->aidx");
  sprintf(cm->aval,"tmpAval.dat");
  sprintf(cm->aptr,"tmpAprt.dat");
  sprintf(cm->aidx,"tmpAidx.dat");
  cm->b=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq2.c, initialize_cmd2(), cm->b");
  sprintf(cm->b,"tmpB.dat");
}

void create_cmatrix(CMD *cm,DOMD *md)
{
  void create_cmatrix_domain(int did,CMD *cm,DOMD *md);
  
  int i;

  for(i=0;i<=md->MN;i++) create_cmatrix_domain(i,cm,md);
}

void create_cmatrix_domain(int did,CMD *cm,DOMD *md)
{
  FILE *fg,*fh,*fdg,*fdh,*fdf;
  double complex *tG,*tH,*tdG,*tdH,CC[9],dCz[9],dCe[9],kc;
  double F,vtz[3],vte[3],dFz,dFe,*tdF;
  size_t Ne,N,t,s,tn,tl,i;
  int td,sd;

  Ne=(size_t)md->bd.sb[did].Ne;
  N=4*Ne;

  tG=(double complex *)m_alloc2(N*N,sizeof(double complex),"solve_bieq2.c, create_cmatrix_domain(),tG"); // malloc
  tH=(double complex *)m_alloc2(N*N,sizeof(double complex),"solve_bieq2.c, create_cmatrix_domain(),tH"); // malloc
  if((fg=fopen(cm->tgfn[did],"wb"))==NULL){    printf("solve_bieq2.c, create_cmatrix_domain(),*fg. Failed to create %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"wb"))==NULL){    printf("solve_bieq2.c, create_cmatrix_domain(),*fh. Failed to create %s file.\n",cm->thfn[did]);    exit(1);  }

  tdG=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"solve_bieq2.c, create_cmatrix_domain(),tdG"); // malloc
  tdH=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"solve_bieq2.c, create_cmatrix_domain(),tdH"); // malloc
  tdF=(double *)m_alloc2(3*N,sizeof(double),"solve_bieq2.c, create_cmatrix_domain(),tdF"); // malloc
  if((fdg=fopen(cm->tdgfn[did],"wb"))==NULL){    printf("solve_bieq2.c, create_cmatrix_domain(),*fdg. Failed to create %s file.\n",cm->tdgfn[did]);    exit(1);  }
  if((fdh=fopen(cm->tdhfn[did],"wb"))==NULL){    printf("solve_bieq2.c, create_cmatrix_domain(),*fdh. Failed to create %s file.\n",cm->tdhfn[did]);    exit(1);  }
  if((fdf=fopen(cm->tdffn[did],"wb"))==NULL){    printf("solve_bieq2.c, create_cmatrix_domain(),*fdf. Failed to create %s file.\n",cm->tdffn[did]);    exit(1);  }


  kc=md->kn[did];
  
  #pragma omp parallel for schedule(dynamic) private(td,tl,tn,vtz,vte,F,dFz,dFe,s,sd,CC,dCz,dCe,i) // omp parallel
  for(t=1;t<=Ne;t++){
    td=md->bd.sb[did].sid[t];
    if( ELT3==check_element_type(td,&(md->bd)) ) tl=3;
    else tl=4;

    for(tn=0;tn<tl;tn++){
      tz_te_bd_node(vtz,vte,td,tn,&(md->bd));

      F=0.0;
      dFz=0.0;
      dFe=0.0;
      for(s=1;s<=Ne;s++){
        sd=md->bd.sb[did].sid[s];

        dcoef_bd_node_t2(CC,dCz,dCe,td,tn,vtz,vte,sd,kc,cm->type,&(md->bd)); // precision check with derivative coefficient

        for(i=0;i<4;i++){
          tG[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+0];
          tH[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+4];

          tdG[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+4];
          tdG[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+4];
        }
        F+=creal(CC[8]);
        dFz+=creal(dCz[8]);
        dFe+=creal(dCe[8]);
      }

      if(did==0){
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=1.0+F;
        tdF[(t-1)*4*3+tn*3+0]=1.0+F;
      }
      else {
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=F;
        tdF[(t-1)*4*3+tn*3+0]=F;
      }
      tdF[(t-1)*4*3+tn*3+1]=dFz;
      tdF[(t-1)*4*3+tn*3+2]=dFe;
    }
  }

  fwrite(tG,sizeof(double complex),N*N,fg);
  fwrite(tH,sizeof(double complex),N*N,fh);
  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);

  fwrite(tdG,sizeof(double complex),2*N*N,fdg);
  fwrite(tdH,sizeof(double complex),2*N*N,fdh);
  fwrite(tdF,sizeof(double),3*N,fdf);
  fclose(fdg);
  fclose(fdh);
  fclose(fdf);
  free(tdG);
  free(tdH);
  free(tdF);
}

void create_Amatrix_csr(DOMD *md)
{
  void create_Amatrix_csr_dac(int did,int cid,FILE *av,FILE *ap,FILE *ai,CMD *cm,DOMD *md);
  
  FILE *av,*ai,*ap;

  size_t did,cid;

  if((av=fopen(md->cm.aval,"wb"))==NULL){    printf("solve_bieq2.c, create_Amatrix_csr(),*av. Failed to create %s file.\n",md->cm.aval);    exit(1);  }
  if((ai=fopen(md->cm.aidx,"wb"))==NULL){    printf("solve_bieq2.c, create_Amatrix_csr(),*ai. Failed to create %s file.\n",md->cm.aidx);    exit(1);  }
  if((ap=fopen(md->cm.aptr,"wb"))==NULL){    printf("solve_bieq2.c, create_Amatrix_csr(),*ap. Failed to create %s file.\n",md->cm.aidx);    exit(1);  }

  md->cm.nnz=0; // initialize nnz
  for(did=0;did<=md->cm.MN;did++)
    for(cid=0;cid<4;cid++) create_Amatrix_csr_dac(did, cid,av,ap,ai,&(md->cm),md);

  fwrite(&(md->cm.nnz),sizeof(size_t),1,ap);

  fclose(av);
  fclose(ai);
  fclose(ap);
}

void create_Amatrix_csr_dac(int did,int cid,FILE *av,FILE *ap,FILE *ai,CMD *cm,DOMD *md)
{
  FILE *fg,*fh;

  double complex *tG,*tH,*tA,k1,k2,delt,beps,sigm;
  double vn[3];
  size_t Ne,t,tn,s,asd,mdid,sdid,etype,l,cc,*ti;
  int td,sd;

  if((fg=fopen(cm->tgfn[did],"rb"))==NULL){    printf("solve_bieq2.c, create_Amatrix_csr_dac(),*fg. Failed to open %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"rb"))==NULL){    printf("solve_bieq2.c, create_Amatrix_csr_dac(),*fh. Failed to open %s file.\n",cm->thfn[did]);    exit(1);  }

  Ne=(size_t)md->bd.sb[did].Ne;
  tG=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"solve_bieq2.c, create_Amatrix_csr_dac(),tG"); // malloc
  tH=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"solve_bieq2.c, create_Amatrix_csr_dac(),tH"); // malloc

  tA=(double complex *)m_alloc2(cm->na,sizeof(double complex),"solve_bieq2.c, create_Amatrix_csr_dac(),tA"); // malloc
  ti=(size_t *)m_alloc2(cm->na,sizeof(size_t),"solve_bieq2.c, create_Amatrix_csr_dac(),ti"); // malloc

  for(t=1;t<=Ne;t++){
    td=md->bd.sb[did].sid[t];

    for(tn=0;tn<4;tn++){
      if(fread(tG,sizeof(double complex),Ne*4,fg)!=Ne*4){
        printf("solve_bieq2.c, create_Amatrix_csr_dac(), failed to read the tG. exit...\n");
        exit(1);
      }
      if(fread(tH,sizeof(double complex),Ne*4,fh)!=Ne*4){
        printf("solve_bieq2.c, create_Amatrix_csr_dac(), failed to read the tH. exit...\n");
        exit(1);
      }
      if( tn==3 && ELT3==check_element_type(td,&(md->bd)) )  continue;

      if(fwrite(&(cm->nnz),sizeof(size_t),1,ap)!=1){ // write A pointer
        printf("solve_bieq2.c, create_Amatrix_csr_dac(), failed to write the nnz. exit...\n");
        exit(1);
      }
      for(l=0;l<cm->na;l++) tA[l]=0.0;

      for(s=1;s<=Ne;s++){
        sd=md->bd.sb[did].sid[s]; // signed element id
        asd=(size_t)abs(sd);
        mdid=(size_t)md->bd.md[asd]; // main domain id
        sdid=(size_t)md->bd.sd[asd]; // sub domain id
        etype=check_element_type(sd,&(md->bd));        

        if(did==mdid){ // main domain
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++) {
              tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++) {
              tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
        } // end main domain
        else { // sub domain
          k1=md->kn[mdid];
          k2=md->kn[sdid];
          delt=k1*k1-k2*k2;
          beps=k1*k1/(k2*k2);
          sigm=1.0-beps;

          if(cid<3){ // U
            if(etype==ELT3){ // linear-triangular element
              for(l=0;l<3;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*(    3) + md->bd.eni[asd][l] ]=delt*vn[cid]*tG[(s-1)*4+l];
                tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]= tG[(s-1)*4+l];
              }
            }
            else { // bi-linear element
              for(l=0;l<4;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*(    3) + md->bd.eni[asd][l] ]=delt*vn[cid]*tG[(s-1)*4+l];
                tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]= tG[(s-1)*4+l];
              }
            }
          } // end U
          else { // varphi
            if(etype==ELT3){ // linear-triangular element
              for(l=0;l<3;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*0 + md->bd.eni[asd][l] ]= sigm*vn[0]*tG[(s-1)*4+l];
                tA[ cm->nn*1 + md->bd.eni[asd][l] ]= sigm*vn[1]*tG[(s-1)*4+l];
                tA[ cm->nn*2 + md->bd.eni[asd][l] ]= sigm*vn[2]*tG[(s-1)*4+l];
                tA[ cm->nn*3 + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*7 + md->bd.eni[asd][l] ]= beps*tG[(s-1)*4+l];
              }
            }
            else { // bi-linear element
              for(l=0;l<4;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*0 + md->bd.eni[asd][l] ]= sigm*vn[0]*tG[(s-1)*4+l];
                tA[ cm->nn*1 + md->bd.eni[asd][l] ]= sigm*vn[1]*tG[(s-1)*4+l];
                tA[ cm->nn*2 + md->bd.eni[asd][l] ]= sigm*vn[2]*tG[(s-1)*4+l];
                tA[ cm->nn*3 + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*7 + md->bd.eni[asd][l] ]= beps*tG[(s-1)*4+l];
              }
            }
          } // end varphi
        } // end sub domain
      } // end for s

      // compress and store data
      cc=0;
      for(l=0;l<cm->na;l++){
        if( creal(tA[l])==0.0 && cimag(tA[l])==0.0) continue;
        tA[cc]=tA[l];
        ti[cc]=l;
        cc+=1;
      }
      fwrite(tA,sizeof(double complex),cc,av);
      fwrite(ti,sizeof(size_t),cc,ai);
      cm->nnz+=cc;

    } // end for tn

  } // end for t

  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);
  free(tA);
  free(ti);
}

void lu_dec_A(CMD *cm)
{
  FILE *fa,*fxa,*fas;
  MKL_Complex16 *A,tc;
  MKL_INT is,ie,i,j,p;

  A=(MKL_Complex16 *)m_alloc2(cm->na*cm->na,sizeof(MKL_Complex16),"create_matrix.c, lu_dec_A(),A");
  // read matrix A
  if((fa =fopen(cm->aval,"rb"))==NULL){     printf("solve_bieq2.c, lu_dec_A(), Failed to open the %s file.\n",cm->aval);    exit(1); }
  if((fxa=fopen(cm->aptr,"rb"))==NULL){     printf("solve_bieq2.c, lu_dec_A(), Failed to open the %s file.\n",cm->aptr);    exit(1); }
  if((fas=fopen(cm->aidx,"rb"))==NULL){     printf("solve_bieq2.c, lu_dec_A(), Failed to open the %s file.\n",cm->aidx);    exit(1); }

  if(fread(&is,sizeof(MKL_INT),1,fxa)!=1){
    printf("solve_bieq2.c, lu_dec_A(), failed to read the is. exit...\n");
    exit(1);
  }
 
  for(j=0;j<cm->na;j++){
    if(fread(&ie,sizeof(MKL_INT),1,fxa)!=1){
      printf("solve_bieq2.c, lu_dec_A(), failed to read the ie1. exit...\n");
      exit(1);
    }
    for(p=is;p<ie;p++){
      if(fread(&i,sizeof(MKL_INT),1,fas)!=1){
        printf("solve_bieq2.c, lu_dec_A(), failed to read the i. exit...\n");
        exit(1);
      }
      if(fread(&tc,sizeof(MKL_Complex16),1,fa)!=1){
        printf("solve_bieq2.c, lu_dec_A(), failed to read the tc. exit...\n");
        exit(1);
      }
      A[j*cm->na+i]=tc;
    }
    is=ie;
  }
  fclose(fa);
  fclose(fxa);
  fclose(fas);

  MKL_INT N,lda,info,*ipiv;

  N=(MKL_INT)cm->na;
  lda=N;
  ipiv=(MKL_INT*)m_alloc2(cm->na,sizeof(MKL_INT),"solve_bieq2.c, lu_dec_A(),ipiv");

  // LU factorization
  info=LAPACKE_zgetrf(LAPACK_ROW_MAJOR , N , N , A , lda , ipiv );
  if(info!=0){
    printf("solve_bieq2.c, lu_dec_A(), LAPCKE_zgetrf() error! info=%d. Exit...\n",(int)info);
    exit(1);
  }
  // inverse matrix
  info=LAPACKE_zgetri(LAPACK_ROW_MAJOR , N , A , lda , ipiv );
  if(info!=0){
    printf("solve_bieq2.c, lu_dec_A(), LAPCKE_zgetri() error! info=%d. Exit...\n",(int)info);
    exit(1);
  }

  // output
  if((fa=fopen(cm->lupfn,"wb"))==NULL){     printf("solve_bieq2.c, lu_dec_A(), Failed to open the %s file.\n",cm->lupfn);    exit(1); }
  fwrite(A,sizeof(MKL_Complex16),cm->na*cm->na,fa); // matrix A
  fclose(fa);

  free(A);
  free(ipiv);
}

void finalize_cmd2(CMD *cm)
{
  int i;

  remove(cm->aval);
  remove(cm->aptr);
  remove(cm->aidx);
  remove(cm->b);

  // free memory
  for(i=0;i<=cm->MN;i++){
    free(cm->tgfn[i]);    free(cm->thfn[i]);
    free(cm->tdgfn[i]);    free(cm->tdhfn[i]);  free(cm->tdffn[i]);
  }
  free(cm->tgfn);  free(cm->thfn);
  free(cm->tdgfn);  free(cm->tdhfn); free(cm->tdffn);
  free(cm->lupfn);

  cm->MN=0;

  free(cm->aval);
  free(cm->aptr);
  free(cm->aidx);
  cm->nn=0;
  cm->na=0;
  cm->nnz=0;
  free(cm->b);
}

void dat_write2(char *fname,DOMD *md)
{
  FILE *fp;
  char tfn[256],pfn[128];
  int i,j,d;

  if((fp=fopen(fname,"wb"))==NULL){    printf("solve_bieq2.c, dat_write2(), Failed to create the %s file.\n",fname);    exit(1);  }

  // fname
  fwrite(md->med_fn,sizeof(char),128,fp);
  fwrite(md->msh_fn,sizeof(char),128,fp);
  // rotation and translation data
  fwrite(md->rv,sizeof(double),3,fp);
  fwrite(&(md->th),sizeof(double),1,fp);
  fwrite(md->tv,sizeof(double),3,fp);
  // material def
  fwrite(&(md->MN),sizeof(int),1,fp);
  fwrite(md->n,sizeof(double complex),md->MN+1,fp);
  fwrite(md->kn,sizeof(double complex),md->MN+1,fp);
  // multi_fbeam
  fwrite(&(md->mw),sizeof(Bobj),1,fp);
  fwrite(md->mw.bd.ipw,sizeof(Ipw),md->mw.n_ipw,fp);
  fwrite(md->mw.bd.fpw,sizeof(Fpw),md->mw.n_fpw,fp);
  fwrite(md->mw.bd.lgb,sizeof(LGb),md->mw.n_lgb,fp);
  fwrite(md->mw.bd.bsb,sizeof(Bsb),md->mw.n_bsb,fp);
  fwrite(md->mw.bd.blg,sizeof(BsLGb),md->mw.n_blg,fp);
  fwrite(md->mw.bd.rab,sizeof(RAb),md->mw.n_rab,fp);
  // BOUD
  fwrite(&(md->bd.Nn),sizeof(int),1,fp);
  fwrite(&(md->bd.Ne),sizeof(int),1,fp);
  for(i=0;i<=md->bd.Nn;i++) fwrite(md->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.ed[i],sizeof(int),4,fp);
  fwrite(&(md->bd.NN),sizeof(int),1,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.eni[i],sizeof(int),4,fp);
  fwrite(md->bd.md,sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.sd,sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.gd,sizeof(int),md->bd.Ne+1,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.wen[i][j],sizeof(double),3,fp);
  // sub domain data
  for(d=0;d<=md->MN;d++){
    fwrite(&(md->bd.sb[d].Ne),sizeof(int),1,fp);
    fwrite(md->bd.sb[d].sid,sizeof(int),md->bd.sb[d].Ne+1,fp);
  }

  // matrix data
  for(i=0;i<128;i++){
    if(fname[i]=='.'){
      pfn[i]='\0';
      break;
    }
    pfn[i]=fname[i];
  }
  for(i=0;i<=md->MN;i++){
    sprintf(tfn,"%s_%s",pfn,md->cm. tgfn[i]+3);  rename(md->cm. tgfn[i],tfn);  fwrite(tfn,sizeof(char),128,fp); // G
    sprintf(tfn,"%s_%s",pfn,md->cm. thfn[i]+3);  rename(md->cm. thfn[i],tfn);  fwrite(tfn,sizeof(char),128,fp); // H
    sprintf(tfn,"%s_%s",pfn,md->cm.tdgfn[i]+3);  rename(md->cm.tdgfn[i],tfn);  fwrite(tfn,sizeof(char),128,fp); // dG
    sprintf(tfn,"%s_%s",pfn,md->cm.tdhfn[i]+3);  rename(md->cm.tdhfn[i],tfn);  fwrite(tfn,sizeof(char),128,fp); // dH
    sprintf(tfn,"%s_%s",pfn,md->cm.tdffn[i]+3);  rename(md->cm.tdffn[i],tfn);  fwrite(tfn,sizeof(char),128,fp); // dF
  }
  sprintf(tfn,"%s_%s",pfn,md->cm.lupfn+3);  rename(md->cm.lupfn,tfn);  fwrite(tfn,sizeof(char),128,fp); // LU + pivot

  fclose(fp);
}

