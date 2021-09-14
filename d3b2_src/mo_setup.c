/*
 * mo_setup.c
 *
 *  Created on: Feb 17, 2019
 *      Author: ohta
 */

#include "bem3_emf_b2.h"

void mo_initialize(int argc,char **argv,MOBJ *mo)
{ 
  void rotation_translation_obj(double *rv,double th,double *tv,DOMD *md); // setup2.c
  void init_elem_const(BOUD *bd); // setup2.c
  void init_boundary_data2(DOMD *md); // setup2.c
  void init_boundary_data_u2(DOMD *md);
  void mo_malloc(MOBJ *mo);
  void dat_read2(char *fname,DOMD *md);
  void renew_beam_data(DOMD *md);
  
  FILE *fp;
  int ti,i,flg;
  char tmp[FNLGTH];

  if(argc!=3){
    printf("This program needs command line argument as follows.\n");
    printf("%s object_setting_file output_data_file_name\n",argv[0]);
    exit(0);
  }

  if((fp=fopen(argv[1],"rt"))==NULL){    printf("mo_read_data(),Can not open the %s file.\n",argv[1]);    exit(1);  }

  fgets(tmp,FNLGTH,fp);

  fscanf(fp,"%d %d\n",&ti,&flg); // number of objects and incident field renew flag
  mo->N=ti;  
  mo_malloc(mo); // malloc
  
  fgets(tmp,256,fp);

  for(i=0;i<mo->N;i++){
    fscanf(fp,"%s %lf %lf %lf",mo->mod_fn[i],&(mo->rs[i][0]),&(mo->rs[i][1]),&(mo->rs[i][2]));
    dat_read2(mo->mod_fn[i],&(mo->md[i]));

    if(flg==1) renew_beam_data(&(mo->md[i]));

    rotation_translation_obj(mo->md[i].rv,0.0,mo->rs[i],&(mo->md[i])); 
    init_elem_const(&(mo->md[i].bd));
    init_boundary_data2(&(mo->md[i]));
    init_boundary_data_u2(&(mo->md[i]));
  }

  fclose(fp);
}

void mo_print_data(MOBJ *mo)
{
  void print_medium_data(DOMD *md); // setup2.c
  void print_mesh_data(DOMD *md); // setup2.c
  
  int i;

  printf("-- multi_fbeam --\n");
  print_data_mfb(&(mo->md[0].mw));

  for(i=0;i<mo->N;i++){
    printf("********** object %d **********\n",i);
    printf("data file name                : %s\n",mo->mod_fn[i]);
    printf("additional translation vector : (% 8.7g,% 8.7g,% 8.7g)\n",mo->rs[i][0],mo->rs[i][1],mo->rs[i][2]);
    print_medium_data(&(mo->md[i]));
    print_mesh_data(&(mo->md[i]));
    if(mo->md[i].th!=0.0 || vabs_d(mo->md[i].tv)!=0.0){
      printf("-- rotation and translation settings --\n");
      if(mo->md[i].th!=0.0){
        printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",mo->md[i].rv[0],mo->md[i].rv[1],mo->md[i].rv[2]);
        printf("rotation angle           [rad]: % 8.7g\n",mo->md[i].th);
      }
      if(vabs_d(mo->md[i].tv)!=0.0){
        printf("translation vector            :(% 8.7g,% 8.7g,% 8.7g)\n",mo->md[i].tv[0],mo->md[i].tv[1],mo->md[i].tv[2]);
      }
    }
  }
  printf("\n");
}

void mo_print_data_mksa(MOBJ *mo)
{
  void print_medium_data(DOMD *md); // setup2.c
  void print_mesh_data(DOMD *md); // setup2.c

  int i;

  printf("-- multi_fbeam --\n");
  print_data_mfb_mksa(&(mo->md[0].mw));

  for(i=0;i<mo->N;i++){
    printf("********** object %d **********\n",i);
    printf("data file name                   : %s\n",mo->mod_fn[i]);
    printf("additional translation vector [m]: (% 8.7g,% 8.7g,% 8.7g)\n",OSUtoMKSA_length(mo->rs[i][0]),OSUtoMKSA_length(mo->rs[i][1]),OSUtoMKSA_length(mo->rs[i][2]));
    print_medium_data(&(mo->md[i]));
    print_mesh_data(&(mo->md[i]));
    if(mo->md[i].th!=0.0 || vabs_d(mo->md[i].tv)!=0.0){
      printf("-- rotation and translation settings --\n");
      if(mo->md[i].th!=0.0){
        printf("vector defining rotation axis    :(% 8.7g,% 8.7g,% 8.7g)\n",mo->md[i].rv[0],mo->md[i].rv[1],mo->md[i].rv[2]);
        printf("rotation angle              [rad]: % 8.7g\n",mo->md[i].th);
      }
      if(vabs_d(mo->md[i].tv)!=0.0){
        printf("translation vector            [m]:(% 8.7g,% 8.7g,% 8.7g)\n",OSUtoMKSA_length(mo->md[i].tv[0]),OSUtoMKSA_length(mo->md[i].tv[1]),OSUtoMKSA_length(mo->md[i].tv[2]));
      }
    }    
  }
  printf("\n");
}

void mo_finalize(MOBJ *mo)
{
  void finalize_cmd2(CMD *cm); // solve_bieq2.c
  void finalize_domd2(DOMD *md); // setup2.c
  
  int N,i;

  N=mo->N;
  for(i=0;i<N;i++){
    finalize_cmd2(&(mo->md[i].cm));
    finalize_domd2(&(mo->md[i]));
    free(mo->mod_fn[i]);
    free(mo->rs[i]);
  }
  free(mo->md);

  free(mo->mod_fn);
  free(mo->rs);
}

void mo_dat_write(char *fname,MOBJ *mo)
{
  FILE *fp;
  int i,j,d,m;

  if((fp=fopen(fname,"wb"))==NULL){    printf("mo_dat_write(), Failed to create the %s file.\n",fname);    exit(1);  }

  fwrite(&(mo->N),sizeof(int),1,fp);
  for(m=0;m<mo->N;m++){
    fwrite(mo->mod_fn[m],sizeof(char),FNLGTH,fp);
    fwrite(mo->rs[m],sizeof(double),3,fp);
    // fname
    fwrite(mo->md[m].med_fn,sizeof(char),128,fp);
    fwrite(mo->md[m].msh_fn,sizeof(char),128,fp);
    // trans 
    fwrite(mo->md[m].rv,sizeof(double),3,fp);
    fwrite(&(mo->md[m].th),sizeof(double),1,fp);
    fwrite(mo->md[m].tv,sizeof(double),3,fp);
    // medium 
    fwrite(&(mo->md[m].MN),sizeof(int),1,fp);
    fwrite(mo->md[m].n,sizeof(double complex),mo->md[m].MN+1,fp);
    fwrite(mo->md[m].kn,sizeof(double complex),mo->md[m].MN+1,fp);
    // multi_fbeam
    fwrite(&(mo->md[m].mw),sizeof(Bobj),1,fp);
    fwrite(mo->md[m].mw.bd.ipw,sizeof(Ipw),mo->md[m].mw.n_ipw,fp);
    fwrite(mo->md[m].mw.bd.fpw,sizeof(Fpw),mo->md[m].mw.n_fpw,fp);
    fwrite(mo->md[m].mw.bd.lgb,sizeof(LGb),mo->md[m].mw.n_lgb,fp);
    fwrite(mo->md[m].mw.bd.bsb,sizeof(Bsb),mo->md[m].mw.n_bsb,fp);
    fwrite(mo->md[m].mw.bd.blg,sizeof(BsLGb),mo->md[m].mw.n_blg,fp);
    fwrite(mo->md[m].mw.bd.rab,sizeof(RAb),mo->md[m].mw.n_rab,fp);
    // BOUD
    fwrite(&(mo->md[m].bd.Nn),sizeof(int),1,fp);
    fwrite(&(mo->md[m].bd.Ne),sizeof(int),1,fp);
    for(i=0;i<=mo->md[m].bd.Nn;i++) fwrite(mo->md[m].bd.rn[i],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fwrite(mo->md[m].bd.ed[i],sizeof(int),4,fp);
    fwrite(&(mo->md[m].bd.NN),sizeof(int),1,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fwrite(mo->md[m].bd.eni[i],sizeof(int),4,fp);
    fwrite(mo->md[m].bd.md,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fwrite(mo->md[m].bd.sd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fwrite(mo->md[m].bd.gd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.ren[i][j],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.wen[i][j],sizeof(double),3,fp);
    // sub domain data
    for(d=0;d<=mo->md[m].MN;d++){
      fwrite(&(mo->md[m].bd.sb[d].Ne),sizeof(int),1,fp);
      fwrite(mo->md[m].bd.sb[d].sid,sizeof(int),mo->md[m].bd.sb[d].Ne+1,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d]. U[i][j],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d].dU[i][j],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d]. E[i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d].dE[i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d]. H[i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d].dH[i][j],sizeof(double complex),3,fp);
    }
  }

  fclose(fp);
}

void mo_dat_read(char *fname,MOBJ *mo)
{
  void mo_malloc(MOBJ *mo);
  void malloc_node(BOUD *bd); // setup2.c
  void malloc_elem(BOUD *bd); // setup2.c
  void init_elem_const(BOUD *bd); // setup2.c
  void malloc_sub_domain(DOMD *md); // setup2.c
  void initalize_cmd2(CMD *cm); // solve_bieq2.c
  
  FILE *fp;
  int i,j,d,tmp,m;

  if((fp=fopen(fname,"rb"))==NULL){    printf("mo_dat_read(), Failed to open the %s file.\n",fname);    exit(1);  }

  fread(&(mo->N),sizeof(int),1,fp);
  mo_malloc(mo);
  for(m=0;m<mo->N;m++){
    fread(mo->mod_fn[m],sizeof(char),FNLGTH,fp);
    fread(mo->rs[m],sizeof(double),3,fp);
    // DOMD data
    // fname
    fread(mo->md[m].med_fn,sizeof(char),128,fp);
    fread(mo->md[m].msh_fn,sizeof(char),128,fp);
    // trans
    fread(mo->md[m].rv,sizeof(double),3,fp);
    fread(&(mo->md[m].th),sizeof(double),1,fp);
    fread(mo->md[m].tv,sizeof(double),3,fp);
    // medium
    fread(&(mo->md[m].MN),sizeof(int),1,fp);
    mo->md[m].n=(double complex *)m_alloc2(mo->md[m].MN+1,sizeof(double complex),"read_medium_data(), mo->md[m].n"); // malloc
    mo->md[m].kn=(double complex *)m_alloc2(mo->md[m].MN+1,sizeof(double complex),"read_medium_data(),mo->md[m].kn"); // malloc
    fread(mo->md[m].n,sizeof(double complex),mo->md[m].MN+1,fp);
    fread(mo->md[m].kn,sizeof(double complex),mo->md[m].MN+1,fp);
    // multi_fbeam
    fread(&(mo->md[m].mw),sizeof(Bobj),1,fp);
    mo->md[m].mw.bd.ipw=(Ipw *)malloc(sizeof(Ipw)*mo->md[m].mw.n_ipw); // malloc
    mo->md[m].mw.bd.fpw=(Fpw *)malloc(sizeof(Fpw)*mo->md[m].mw.n_fpw); // malloc
    mo->md[m].mw.bd.lgb=(LGb *)malloc(sizeof(LGb)*mo->md[m].mw.n_lgb); // malloc
    mo->md[m].mw.bd.bsb=(Bsb *)malloc(sizeof(Bsb)*mo->md[m].mw.n_bsb); // malloc
    mo->md[m].mw.bd.blg=(BsLGb *)malloc(sizeof(BsLGb)*mo->md[m].mw.n_blg); // malloc
    mo->md[m].mw.bd.rab=(RAb *)malloc(sizeof(RAb)*mo->md[m].mw.n_rab); // malloc
    fread(mo->md[m].mw.bd.ipw,sizeof(Ipw),mo->md[m].mw.n_ipw,fp);
    fread(mo->md[m].mw.bd.fpw,sizeof(Fpw),mo->md[m].mw.n_fpw,fp);
    fread(mo->md[m].mw.bd.lgb,sizeof(LGb),mo->md[m].mw.n_lgb,fp);
    fread(mo->md[m].mw.bd.bsb,sizeof(Bsb),mo->md[m].mw.n_bsb,fp);
    fread(mo->md[m].mw.bd.blg,sizeof(BsLGb),mo->md[m].mw.n_blg,fp);
    fread(mo->md[m].mw.bd.rab,sizeof(RAb),mo->md[m].mw.n_rab,fp);
    setup_mfb(&(mo->md[m].mw));
    // BOUD
    fread(&(mo->md[m].bd.Nn),sizeof(int),1,fp);
    fread(&(mo->md[m].bd.Ne),sizeof(int),1,fp);
    malloc_node(&(mo->md[m].bd)); // malloc
    malloc_elem(&(mo->md[m].bd)); // malloc
    for(i=0;i<=mo->md[m].bd.Nn;i++) fread(mo->md[m].bd.rn[i],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fread(mo->md[m].bd.ed[i],sizeof(int),4,fp);
    fread(&(mo->md[m].bd.NN),sizeof(int),1,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fread(mo->md[m].bd.eni[i],sizeof(int),4,fp);
    fread(mo->md[m].bd.md,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fread(mo->md[m].bd.sd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fread(mo->md[m].bd.gd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.ren[i][j],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.wen[i][j],sizeof(double),3,fp);
    init_elem_const(&(mo->md[m].bd)); // setup
    // sub domain data
    malloc_sub_domain(&(mo->md[m])); // malloc
    for(d=0;d<=mo->md[m].MN;d++){
      fread(&tmp,sizeof(int),1,fp);
      fread(mo->md[m].bd.sb[d].sid,sizeof(int),mo->md[m].bd.sb[d].Ne+1,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d]. U[i][j],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d].dU[i][j],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d]. E[i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d].dE[i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d]. H[i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d].dH[i][j],sizeof(double complex),3,fp);
    }

    mo->md[m].cm.MN=mo->md[m].MN;
    mo->md[m].cm.na=8*mo->md[m].bd.NN;
    initalize_cmd2(&(mo->md[m].cm));
  }

  fclose(fp);
}

void mo_output_node_particles(char *fname,MOBJ *mo)
{
  FILE *fp;
  int s1,s2,oid,i,j;
  char *sd,fo[128]="";

  sd=strrchr(fname,'.');
  if(sd==NULL){ // no file extension
    sprintf(fo,"%s.particles",fname);
  }
  else {
    s1=strlen(fname);
    s2=strlen(sd);
    strncpy(fo,fname,s1-s2);
    sprintf(fo,"%s.particles",fo);
  }
  
  if((fp=fopen(fo,"wt"))==NULL){    printf("Can not open the %s file.\n",fo);    exit(1);  }
  fprintf(fp,"# x y z object_id\n");
  
  for(oid=0;oid<mo->N;oid++){
    for(i=1;i<=mo->md[oid].bd.Ne;i++){
      for(j=0;j<4;j++){
        fprintf(fp,"%15.14e %15.14e %15.14e %d\n",mo->md[oid].bd.ren[i][j][0],mo->md[oid].bd.ren[i][j][1],mo->md[oid].bd.ren[i][j][2],oid);
      }
    }
  }

  fclose(fp);
}

////////////////////////////////////////////////////////////////////
void mo_malloc(MOBJ *mo)
{
  int i,N;

  N=mo->N;
  mo->md=(DOMD *)m_alloc2(N,sizeof(DOMD),"d3_md_c1.c, mo_malloc(), mo->md");

  mo->mod_fn=(char **)m_alloc2(N,sizeof(char *),"d3_md_c1.c, mo_malloc(), mo->mod_fn");
  mo->rs=(double **)m_alloc2(N,sizeof(double *),"d3_md_c1.c, mo_malloc(), mo->rs");
  for(i=0;i<N;i++){
    mo->mod_fn[i]=(char *)m_alloc2(FNLGTH,sizeof(char),"d3_md_c1.c, mo_malloc(), mo->mod_fn[i]");
    mo->rs[i]=(double *)m_alloc2(3,sizeof(double),"d3_md_c1.c, mo_malloc(), mo->rs[i]");
  }

}

void renew_beam_data(DOMD *md)
{
  double lambda0,n_0;

  lambda0=md->mw.lambda_0;
  n_0=md->mw.n_0;

  free_mfb(&(md->mw));
  init_mfb(&(md->mw));
  read_data_mfb(&(md->mw));
  setup_mfb(&(md->mw));

  // check refractive index
  if(lambda0!=md->mw.lambda_0 || n_0!=md->mw.n_0){
    printf("solve_bieq2.c, renew_beam_data(), wave number check error!\n");
    printf("multi wave old data lambda0=%g, n_0=%g\n",lambda0,n_0);
    printf("multi wave new data lambda0=%g, n_0=%g\n",md->mw.lambda_0,md->mw.n_0);
    printf("Exit...\n");
    exit(1);
  }
}

void dat_read2(char *fname,DOMD *md)
{
  void malloc_node(BOUD *bd); // setup2.c
  void malloc_elem(BOUD *bd); // setup2.c
  void init_elem_const(BOUD *bd); // setup2.c
  void malloc_sub_domain(DOMD *md); // setup2.c
  void initalize_cmd2(CMD *cm); // solve_bieq2.c
  
  FILE *fp;
  int i,j,d,tmp;

  if((fp=fopen(fname,"rb"))==NULL){    printf("dat_read(), Failed to open the %s file.\n",fname);    exit(1);  }

  // fname
  fread(md->med_fn,sizeof(char),128,fp);
  fread(md->msh_fn,sizeof(char),128,fp);
  // trans param
  fread(md->rv,sizeof(double),3,fp);
  fread(&(md->th),sizeof(double),1,fp);
  fread(md->tv,sizeof(double),3,fp);
  // medium def
  fread(&(md->MN),sizeof(int),1,fp);
  md->n=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"read_medium_data(), md->n"); // malloc
  md->kn=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"read_medium_data(),md->kn"); // malloc
  fread(md->n,sizeof(double complex),md->MN+1,fp);
  fread(md->kn,sizeof(double complex),md->MN+1,fp);
  // multi_fbeam
  fread(&(md->mw),sizeof(Bobj),1,fp);
  md->mw.bd.ipw=(Ipw *)malloc(sizeof(Ipw)*md->mw.n_ipw); // malloc
  md->mw.bd.fpw=(Fpw *)malloc(sizeof(Fpw)*md->mw.n_fpw); // malloc
  md->mw.bd.lgb=(LGb *)malloc(sizeof(LGb)*md->mw.n_lgb); // malloc
  md->mw.bd.bsb=(Bsb *)malloc(sizeof(Bsb)*md->mw.n_bsb); // malloc
  md->mw.bd.blg=(BsLGb *)malloc(sizeof(BsLGb)*md->mw.n_blg); // malloc
  md->mw.bd.rab=(RAb *)malloc(sizeof(RAb)*md->mw.n_rab); // malloc
  fread(md->mw.bd.ipw,sizeof(Ipw),md->mw.n_ipw,fp);
  fread(md->mw.bd.fpw,sizeof(Fpw),md->mw.n_fpw,fp);
  fread(md->mw.bd.lgb,sizeof(LGb),md->mw.n_lgb,fp);
  fread(md->mw.bd.bsb,sizeof(Bsb),md->mw.n_bsb,fp);
  fread(md->mw.bd.blg,sizeof(BsLGb),md->mw.n_blg,fp);
  fread(md->mw.bd.rab,sizeof(RAb),md->mw.n_rab,fp);
  setup_mfb(&(md->mw));
  // BOUD
  fread(&(md->bd.Nn),sizeof(int),1,fp);
  fread(&(md->bd.Ne),sizeof(int),1,fp);
  malloc_node(&(md->bd)); // malloc
  malloc_elem(&(md->bd)); // malloc
  for(i=0;i<=md->bd.Nn;i++) fread(md->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fread(md->bd.ed[i],sizeof(int),4,fp);
  fread(&(md->bd.NN),sizeof(int),1,fp);
  for(i=0;i<=md->bd.Ne;i++) fread(md->bd.eni[i],sizeof(int),4,fp);
  fread(md->bd.md,sizeof(int),md->bd.Ne+1,fp);
  fread(md->bd.sd,sizeof(int),md->bd.Ne+1,fp);
  fread(md->bd.gd,sizeof(int),md->bd.Ne+1,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd.wen[i][j],sizeof(double),3,fp);
  init_elem_const(&(md->bd)); // setup
  // sub domain data
  malloc_sub_domain(md); // malloc
  for(d=0;d<=md->MN;d++){
    fread(&tmp,sizeof(int),1,fp);
    fread(md->bd.sb[d].sid,sizeof(int),md->bd.sb[d].Ne+1,fp);
  }

  // coefficient matrix data
  md->cm.MN=md->MN;
  md->cm.na=8*md->bd.NN;
  initalize_cmd2(&(md->cm));
  for(i=0;i<=md->MN;i++){
    fread(md->cm.tgfn[i],sizeof(char),128,fp); // G
    fread(md->cm.thfn[i],sizeof(char),128,fp); // H
    fread(md->cm.tdgfn[i],sizeof(char),128,fp); // dG
    fread(md->cm.tdhfn[i],sizeof(char),128,fp); // dH
    fread(md->cm.tdffn[i],sizeof(char),128,fp); // dF
  }
  fread(md->cm.lupfn,sizeof(char),128,fp); // LU + pivot

  fclose(fp);
}

void init_boundary_data_u2(DOMD *md)
{
  double complex h[3],dh[3];
  double cr[3][4],cw[3][3],r[3],w[3];
  int i,j,l,m;

  // element node and incident field data
  for(i=1;i<=md->bd.Ne;i++){
    for(l=0;l<3;l++){
      for(m=0;m<4;m++) cr[l][m]=md->bd.cr[i][l][m];
      for(m=0;m<3;m++) cw[l][m]=md->bd.cw[i][l][m];
    }
    if(md->bd.ed[i][3]==0){ // linear triangular element
      for(j=0;j<4;j++){
        lit_rw_zeta_eta(r,w,md->bd.zt_34[j],md->bd.et_34[j],cr,cw);
        for(l=0;l<3;l++){
          md->bd.ren[i][j][l]=r[l];
          md->bd.wen[i][j][l]=w[l];
        }
        if(md->bd.md[i]==0){
          vuni_d(w);
          calc_mfb_EH_dv(md->bd.Ui[i][j],h,md->bd.dUi[i][j],dh,r,w,&(md->mw));
          md->bd.Ui[i][j][3]=0.0;
          md->bd.dUi[i][j][3]=0.0;
          for(l=0;l<4;l++){
            md->bd. Ui0[i][j][l]=md->bd. Ui[i][j][l];
            md->bd.dUi0[i][j][l]=md->bd.dUi[i][j][l];
            md->bd. Uip[i][j][l]=md->bd. Ui[i][j][l];
          }
        }
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        bil_rw_zeta_eta(r,w,md->bd.zt_44[j],md->bd.et_44[j],cr,cw);
        for(l=0;l<3;l++){
          md->bd.ren[i][j][l]=r[l];
          md->bd.wen[i][j][l]=w[l];
        }
        if(md->bd.md[i]==0){
          vuni_d(w);
          calc_mfb_EH_dv(md->bd.Ui[i][j],h,md->bd.dUi[i][j],dh,r,w,&(md->mw));
          md->bd.Ui[i][j][3]=0.0;
          md->bd.dUi[i][j][3]=0.0;
          for(l=0;l<4;l++){
            md->bd. Ui0[i][j][l]=md->bd. Ui[i][j][l];
            md->bd.dUi0[i][j][l]=md->bd.dUi[i][j][l];
            md->bd. Uip[i][j][l]=md->bd. Ui[i][j][l];
          }
        }
      }
    }
  }
}

