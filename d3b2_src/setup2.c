/*
 * setup2.c
 *
 *  Created on: Feb 16, 2019
 *      Author: ohta
 */

#include "bem3_emf_b2.h"


void read_domd2(int argc,char **argv,DOMD *md)
{
  void filename_chk(int argc,char **argv);
  void read_medium_data(char *med_fn,DOMD *md);
  void read_mesh_data(char *msh_fn,DOMD *md);
  
  filename_chk(argc,argv);
  // multi_fbeam
  init_mfb(&(md->mw));
  read_data_mfb(&(md->mw));
  
  read_medium_data(argv[1],md);
  read_mesh_data(argv[2],md);
  
  if(argc==11){
    md->rv[0]=atof(argv[ 4]);
    md->rv[1]=atof(argv[ 5]);
    md->rv[2]=atof(argv[ 6]);
    md->th   =atof(argv[ 7]);
    md->tv[0]=atof(argv[ 8]);
    md->tv[1]=atof(argv[ 9]);
    md->tv[2]=atof(argv[10]);
  }
  else {
    md->rv[0]=1.0;
    md->rv[1]=0.0;
    md->rv[2]=0.0;
    md->th=0.0;
    md->tv[0]=0.0;
    md->tv[1]=0.0;
    md->tv[2]=0.0;
  }
}

void print_domd2(DOMD *md)
{
  void print_medium_data(DOMD *md);
  void print_mesh_data(DOMD *md);
  
  printf("-- multi_fbeam data --\n");
  print_data_mfb(&(md->mw));

  print_medium_data(md);
  printf("\n");

  print_mesh_data(md);
  printf("\n");
  
  if(md->th!=0.0 || vabs_d(md->tv)!=0.0){
    printf("-- rotation and translation settings --\n");
    if(md->th!=0.0){
      printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",md->rv[0],md->rv[1],md->rv[2]);
      printf("rotation angle           [rad]: %8.7g\n",md->th);
    }
    if(vabs_d(md->tv)!=0.0){
      printf("translation vector            :(%8.7g,%8.7g,%8.7g)\n",md->tv[0],md->tv[1],md->tv[2]);
    }
    printf("\n");
  }
}

void print_domd2_mksa(DOMD *md)
{
  void print_medium_data(DOMD *md);
  void print_mesh_data(DOMD *md);

  printf("-- multi_fbeam data --\n");
  print_data_mfb_mksa(&(md->mw));

  print_medium_data(md);
  printf("\n");

  print_mesh_data(md);
  printf("\n");
  
  if(md->th!=0.0 || vabs_d(md->tv)!=0.0){
    printf("-- rotation and translation settings --\n");
    if(md->th!=0.0){
      printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",md->rv[0],md->rv[1],md->rv[2]);
      printf("rotation angle           [rad]: %8.7g\n",md->th);
    }
    if(vabs_d(md->tv)!=0.0){
      printf("translation vector         [m]:(%8.7g,%8.7g,%8.7g)\n",OSUtoMKSA_length(md->tv[0]),OSUtoMKSA_length(md->tv[1]),OSUtoMKSA_length(md->tv[2]));
    }
    printf("\n");
  }
}

void initialize_domd2(DOMD *md)
{
  void rotation_translation_obj(double *rv,double th,double *tv,DOMD *md); 
  void init_elem_const(BOUD *bd); 
  void malloc_sub_domain(DOMD *md); 
  void init_sub_domain(DOMD *md); 
  void init_boundary_data2(DOMD *md);

  int i;

  // multi_fbeam
  setup_mfb(&(md->mw));
  // medium
  md->n[0]=md->mw.n_0;
  for(i=0;i<=md->MN;i++) md->kn[i]=2.0*M_PI/md->mw.lambda_0*md->n[i];
  // rotation and translation object
  if(md->th!=0.0 || vabs_d(md->tv)!=0.0) rotation_translation_obj(md->rv,md->th,md->tv,md);  
  // element constant
  init_elem_const(&(md->bd));
  init_boundary_data2(md);
  // sub domain
  malloc_sub_domain(md);
  init_sub_domain(md);
}

void finalize_domd2(DOMD *md)
{
  void mfree_sub_domain(DOMD *md);
  void mfree_elem(BOUD *bd);
  void mfree_node(BOUD *bd);
  
  mfree_sub_domain(md);

  mfree_elem(&(md->bd));
  mfree_node(&(md->bd));
  free(md->n);
  free(md->kn);

  free_mfb(&(md->mw));
}


/////////////////////////////////////////////////////////////////////////////////////////////
void filename_chk(int argc,char **argv)
{
  if(argc!=4 && argc!=11){
    printf("This program needs command line arguments as follows.\n");
    printf("%s medium_datafile_name mesh_datafile_name output_datafile_name [rv_x rv_y rv_z theta tr_x tr_y tr_z ](optional)\n",argv[0]);
    printf("rv : vector defining rotation axis, theta : rotation angle ( using Rodrigues' rotation formula ), tr : translation vector\n"); 
    printf("Exit...\n");
    exit(0);
  }
}

void read_medium_data(char *med_fn,DOMD *md)
{
  FILE *fp;
  double td,td2;
  char buf[256]="";
  int i,ti;

  if((fp=fopen(med_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",med_fn);    exit(1);  }
  strcpy(md->med_fn,med_fn);
  if(fgets(buf,256,fp)==NULL){
    printf("setup2.c, read_medium_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("setup2.c, read_medium_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%d\n",&ti)!=1){
    printf("setup2.c, read_medium_data(), failed to read the MN. exit...\n");
    exit(1);
  }
  md->MN=ti;
  md-> n=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"setup2.c, read_medium_data(),md->n");
  md->kn=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"setup2.c, read_medium_data(),md->kn");
  if(fgets(buf,256,fp)==NULL){
    printf("setup2.c, read_medium_data(), failed to read the line. exit...\n");
    exit(1);
  }
  for(i=1;i<=md->MN;i++){
    if(fscanf(fp,"%lf",&td)!=1){
      printf("setup2.c, read_medium_data(), failed to read the real(n). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&td2)!=1){
      printf("setup2.c, read_medium_data(), failed to read the imag(n). exit...\n");
      exit(1);
    }
    md->n[i]=td+td2*I;
  }
  fclose(fp);
}

void print_medium_data(DOMD *md)
{
  int i;
  printf("-- medium data --\n");
  printf("medium data file name                  : %s\n",md->med_fn);
  for(i=1;i<=md->MN;i++){
    printf("medium (domain) id %2d refractive index :%8.7g + %8.7gI\n",i,creal(md->n[i]),cimag(md->n[i]));
  }
}

void read_mesh_data(char *msh_fn,DOMD *md)
{
  void malloc_node(BOUD *bd);
  void malloc_elem(BOUD *bd);
  
  FILE *fp;
  char buf[256]="";
  double td;
  int ti,i,j,ti2,etype,tmpi,nc;

  if((fp=fopen(msh_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",msh_fn);    exit(1);  }
  strcpy(md->msh_fn,msh_fn);
  if(fgets(buf,256,fp)==NULL){
    printf("setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td)!=1){
    printf("setup.c, read_mesh_data(), failed to read the msh version. exit...\n");
    exit(1);
  }
  // check file version
  if(td<MSHVER){
    printf("This program supports mesh file version %g later. Reading file version is %g. Exit...\n",MSHVER,td);
    fclose(fp);
    exit(1);
  }
  if(fscanf(fp,"%d",&ti)!=1){
    printf("setup.c, read_mesh_data(), failed to read the msh formant. exit...\n");
    exit(1);
  }
  //check data format
  if(ti!=MSHASCI){
    printf("This program supports 'ASCII' data format mesh file. Exit...\n");
    fclose(fp);
    exit(1);
  }
  if(fscanf(fp,"%d\n",&ti)!=1){
    printf("setup.c, read_mesh_data(), failed to read the data precision. exit...\n");
    exit(1);
  }
  //check data precision
  if(ti!=MSHPREC){
    printf("This program supports double precision mesh data. Exit...\n");
    fclose(fp);
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }

  if(fscanf(fp,"%d\n",&ti)!=1){
    printf("setup.c, read_mesh_data(), failed to read the Nn. exit...\n");
    exit(1);
  } 
  md->bd.Nn=ti;
  malloc_node(&(md->bd));
  for(i=1;i<=md->bd.Nn;i++){
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the id. exit...\n");
      exit(1);
    }
    if(ti!=i)       printf("bad id %d\n",ti);
    if(fscanf(fp,"%lf",&td)!=1){
      printf("setup.c, read_mesh_data(), failed to read the rn[i][0]. exit...\n");
      exit(1);
    }
    md->bd.rn[i][0]=td;
    if(fscanf(fp,"%lf",&td)!=1){
      printf("setup.c, read_mesh_data(), failed to read the rn[i][1]. exit...\n");
      exit(1);
    }
    md->bd.rn[i][1]=td;
    if(fscanf(fp,"%lf\n",&td)!=1){
      printf("setup.c, read_mesh_data(), failed to read the rn[i][2]. exit...\n");
      exit(1);
    }
    md->bd.rn[i][2]=td;
  }
  if(fgets(buf,256,fp)==NULL){
    printf("setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("setup.c, read_mesh_data(), failed to read the line. exit...\n");
    exit(1);
  }

  if(fscanf(fp,"%d",&ti)!=1){
    printf("setup.c, read_mesh_data(), failed to read the Ne. exit...\n");
    exit(1);
  }
  md->bd.Ne=ti/2; 
  malloc_elem(&(md->bd));

  nc=0;
  for(i=1;i<=md->bd.Ne;i++){
    // element id
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the id. exit...\n");
      exit(1);
    }
    if(ti!=i*2-1){
      printf("bad id :%d. Exit...\n",ti);
      exit(1);
    }
    // element type
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the type. exit...\n");
      exit(1);
    }
    etype=ti;
    if(ti!=ELT3 && ti!=ELT4){
      printf("bad element type. element type must be %d or %d. Exit...\n",ELT3,ELT4);
      exit(1);
    }
    // number of tags
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the number of tags. exit...\n");
      exit(1);
    }
    for(j=0;j<ti;j++){
      if(fscanf(fp,"%d",&ti2)!=1){
        printf("setup.c, read_mesh_data(), failed to read the physical entity. exit...\n");
        exit(1);
      }
      if(j==0){ // domain id ( gmsh physical entity)
        if(ti2==OPENDID) ti2=0;
        if(md->MN>=ti2) md->bd.md[i]=ti2;
        else {
          printf("domain id %d is not defined medium data. check domain and medium data. exit..\n",ti2);
          exit(1);
        }
      }
      else if(j==1){ // group id ( elementary geometrical entity )
        md->bd.gd[i]=ti2;
      }
    }
    // node id
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the ed[i][0]. exit...\n");
      exit(1);
    }
    md->bd.ed[i][0]=ti;
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the ed[i][1]. exit...\n");
      exit(1);
    }
    md->bd.ed[i][1]=ti;
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the ed[i][2]. exit...\n");
      exit(1);
    }
    md->bd.ed[i][2]=ti;
    if(etype==ELT3) md->bd.ed[i][3]=0;
    else {
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][3]. exit...\n");
        exit(1);
      }
      md->bd.ed[i][3]=ti;
    }
    // element node id
    if(etype==ELT3){
      md->bd.eni[i][0]=nc++;
      md->bd.eni[i][1]=nc++;
      md->bd.eni[i][2]=nc++;
      md->bd.eni[i][3]=-1;
    }
    else {
      md->bd.eni[i][0]=nc++;
      md->bd.eni[i][1]=nc++;
      md->bd.eni[i][2]=nc++;
      md->bd.eni[i][3]=nc++;
    }

    // element id
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the element id. exit...\n");
      exit(1);
    }
    if(ti!=i*2){
      printf("bad id :%d. Exit...\n",ti);
      exit(1);
    }
    // element type
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the element type. exit...\n");
      exit(1);
    }
    etype=ti;
    if(ti!=ELT3 && ti!=ELT4){
      printf("bad element type. element type must be %d or %d. Exit...\n",ELT3,ELT4);
      exit(1);
    }
    // number of tags
    if(fscanf(fp,"%d",&ti)!=1){
      printf("setup.c, read_mesh_data(), failed to read the number of tags. exit...\n");
      exit(1);
    }
    for(j=0;j<ti;j++){
      if(fscanf(fp,"%d",&ti2)!=1){
        printf("setup.c, read_mesh_data(), failed to read the domain id. exit...\n");
        exit(1);
      }
      if(j==0){ // domain id
        if(ti2==OPENDID) ti2=0;
        if(md->MN>=ti2) md->bd.sd[i]=ti2;
        else {
          printf("domain id %d is not defined medium data. check domain and medium data! exit..\n",ti2);
          exit(1);
        }
      }
    }
    // check node id
    if(etype==ELT3){
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][0]. exit...\n");
        exit(1);
      }
      if(md->bd.ed[i][0]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][1]. exit...\n");
        exit(1);
      }
      if(md->bd.ed[i][2]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][2]. exit...\n");
        exit(1);
      }
      if(md->bd.ed[i][1]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
    }
    else {
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][0]. exit...\n");
        exit(1);
      }
      if(md->bd.ed[i][0]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][1]. exit...\n");
        exit(1);
      }
      if(md->bd.ed[i][3]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][2]. exit...\n");
        exit(1);
      }
      if(md->bd.ed[i][2]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      if(fscanf(fp,"%d",&ti)!=1){
        printf("setup.c, read_mesh_data(), failed to read the ed[i][3]. exit...\n");
        exit(1);
      }
      if(md->bd.ed[i][1]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
    }
    // exchange open region domain to main domain
    if(md->bd.sd[i]==0){
      md->bd.sd[i]=md->bd.md[i];
      md->bd.md[i]=0;
      if(etype==ELT3){
        tmpi=md->bd.ed[i][1];
        md->bd.ed[i][1]=md->bd.ed[i][2];
        md->bd.ed[i][2]=tmpi;
      }
      else {
        tmpi=md->bd.ed[i][3];
        md->bd.ed[i][3]=md->bd.ed[i][1];
        md->bd.ed[i][1]=tmpi;
      }
    }
  }
  fclose(fp);
  md->bd.NN=nc;
 
}

void print_mesh_data(DOMD *md)
{
  printf("-- mesh data --\n");
  printf("mesh data file name    : %s\n",md->msh_fn);
  printf("node number            : %8d\n",md->bd.Nn);
  printf("defined element number : %8d\n",md->bd.Ne*2);
}

void malloc_node(BOUD *bd)
{
  int i,N=bd->Nn;

  bd->rn=(double **)m_alloc2(N+1,sizeof(double*),"setup2.c, malloc_node(), bd->rn");
  for(i=0;i<=N;i++){
    bd->rn[i]=(double *)m_alloc2(3,sizeof(double),"setup2.c, malloc_node(), bd->rn[i]");
  }
}

void mfree_node(BOUD *bd)
{
  int i,N=bd->Nn;

  for(i=0;i<=N;i++) free(bd->rn[i]);
  free(bd->rn);
  bd->Nn=0;
}

void malloc_elem(BOUD *bd)
{
  int i,j,Ne=bd->Ne;

  bd->ed =(int **)m_alloc2(Ne+1,sizeof(int *),"setup2.c, malloc_elem(), bd->ed");
  bd->eni=(int **)m_alloc2(Ne+1,sizeof(int *),"setup2.c, malloc_elem(), bd->eni");
  for(i=0;i<=Ne;i++){
    bd->ed [i]=(int *)m_alloc2(4,sizeof(int ),"setup2.c, malloc_elem(), bd->ed[i]");
    bd->eni[i]=(int *)m_alloc2(4,sizeof(int ),"setup2.c, malloc_elem(), bd->eni[i]");
  }

  bd->md=(int *)m_alloc2(Ne+1,sizeof(int),"setup2.c, malloc_elem(), bd->md");
  bd->sd=(int *)m_alloc2(Ne+1,sizeof(int),"setup2.c, malloc_elem(), bd->sd");
  bd->gd=(int *)m_alloc2(Ne+1,sizeof(int),"setup2.c, malloc_elem(), bd->gd");

  // element constant
  bd->cr =(double ***)m_alloc2(Ne+1,sizeof(double **),"setup2.c, malloc_elem(), bd->cr");
  bd->cw =(double ***)m_alloc2(Ne+1,sizeof(double **),"setup2.c, malloc_elem(), bd->cw");
  bd->ren=(double ***)m_alloc2(Ne+1,sizeof(double **),"setup2.c, malloc_elem(), bd->ren");
  bd->wen=(double ***)m_alloc2(Ne+1,sizeof(double **),"setup2.c, malloc_elem(), bd->wen");
  bd-> Ui=(double complex ***)m_alloc2(Ne+1,sizeof(double complex **),"setup2.c, malloc_elem(), bd->Ui");
  bd->dUi=(double complex ***)m_alloc2(Ne+1,sizeof(double complex **),"setup2.c, malloc_elem(), bd->dUi");
  bd-> Ui0=(double complex ***)m_alloc2(Ne+1,sizeof(double complex **),"setup2.c, malloc_elem(), bd->Ui0");
  bd->dUi0=(double complex ***)m_alloc2(Ne+1,sizeof(double complex **),"setup2.c, malloc_elem(), bd->dUi0");
  bd-> Uip=(double complex ***)m_alloc2(Ne+1,sizeof(double complex **),"setup2.c, malloc_elem(), bd->Uip");
  for(i=0;i<=Ne;i++){
    bd->cr[i]=(double **)m_alloc2(3,sizeof(double *),"setup2.c, malloc_elem(), bd->cr[i]");
    bd->cw[i]=(double **)m_alloc2(3,sizeof(double *),"setup2.c, malloc_elem(), bd->cw[i]");
    for(j=0;j<3;j++){
      bd->cr[i][j]=(double *)m_alloc2(4,sizeof(double),"setup2.c, malloc_elem(), bd->cr[i][j]");
      bd->cw[i][j]=(double *)m_alloc2(3,sizeof(double),"setup2.c, malloc_elem(), bd->cw[i][j]");
    }

    bd->ren[i]=(double **)m_alloc2(4,sizeof(double *),"setup2.c, malloc_elem(), bd->ren[i]");
    bd->wen[i]=(double **)m_alloc2(4,sizeof(double *),"setup2.c, malloc_elem(), bd->wen[i]");
    bd-> Ui[i]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_elem(), bd->Ui[i]");
    bd->dUi[i]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_elem(), bd->dUi[i]");
    bd-> Ui0[i]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_elem(), bd->Ui0[i]");
    bd->dUi0[i]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_elem(), bd->dUi0[i]");
    bd-> Uip[i]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_elem(), bd->Uip[i]");
    for(j=0;j<4;j++){
      bd->ren[i][j]=(double *)m_alloc2(3,sizeof(double),"setup2.c, malloc_elem(), bd->ren[i][j]");
      bd->wen[i][j]=(double *)m_alloc2(3,sizeof(double),"setup2.c, malloc_elem(), bd->wen[i][j]");
      bd-> Ui[i][j]=(double complex *)m_alloc2(4,sizeof(double complex),"setup2.c, malloc_elem(), bd->Ui[i][j]");
      bd->dUi[i][j]=(double complex *)m_alloc2(4,sizeof(double complex),"setup2.c, malloc_elem(), bd->dUi[i][j]");
      bd-> Ui0[i][j]=(double complex *)m_alloc2(4,sizeof(double complex),"setup2.c, malloc_elem(), bd->Ui0[i][j]");
      bd->dUi0[i][j]=(double complex *)m_alloc2(4,sizeof(double complex),"setup2.c, malloc_elem(), bd->dUi0[i][j]");
      bd-> Uip[i][j]=(double complex *)m_alloc2(4,sizeof(double complex),"setup2.c, malloc_elem(), bd->Uip[i][j]");
    }
  }
}

void mfree_elem(BOUD *bd)
{
  int i,j,Ne=bd->Ne;

  for(i=0;i<=Ne;i++){
    free(bd->ed[i]);
    free(bd->eni[i]);
  }
  free(bd->ed);
  free(bd->eni);

  free(bd->md);
  free(bd->sd);
  free(bd->gd);

  for(i=0;i<=Ne;i++){
    for(j=0;j<3;j++){
      free(bd->cr[i][j]);
      free(bd->cw[i][j]);
    }
    free(bd->cr[i]);
    free(bd->cw[i]);

    for(j=0;j<4;j++){
      free(bd->ren[i][j]);
      free(bd->wen[i][j]);
      free(bd->Ui[i][j]);
      free(bd->dUi[i][j]);
      free(bd->Ui0[i][j]);
      free(bd->dUi0[i][j]);
      free(bd->Uip[i][j]);
    }
    free(bd->ren[i]);
    free(bd->wen[i]);
    free(bd->Ui[i]);
    free(bd->dUi[i]);
    free(bd->Ui0[i]);
    free(bd->dUi0[i]);
    free(bd->Uip[i]);
  }
  free(bd->cr);
  free(bd->cw);
  free(bd->ren);
  free(bd->wen);
  free(bd->Ui);
  free(bd->dUi);
  free(bd->Ui0);
  free(bd->dUi0);
  free(bd->Uip);

  bd->Ne=0;
}

void rotation_translation_obj(double *rv,double th,double *tv,DOMD *md)
{
  double ct,st,r[3],M[9],nv[3];
  size_t s,i;

  nv[0]=rv[0];
  nv[1]=rv[1];
  nv[2]=rv[2];
  vuni_d(nv);

  // rotation matrix
  st=sin(th);
  ct=cos(th);
  M[0]=ct+nv[0]*nv[0]*(1.0-ct);
  M[1]=nv[0]*nv[1]*(1.0-ct)-nv[2]*st;
  M[2]=nv[2]*nv[0]*(1.0-ct)+nv[1]*st;
  M[3]=nv[0]*nv[1]*(1.0-ct)+nv[2]*st;
  M[4]=ct+nv[1]*nv[1]*(1.0-ct);
  M[5]=nv[1]*nv[2]*(1.0-ct)-nv[0]*st;
  M[6]=nv[2]*nv[0]*(1.0-ct)-nv[1]*st;
  M[7]=nv[1]*nv[2]*(1.0-ct)+nv[0]*st;
  M[8]=ct+nv[2]*nv[2]*(1.0-ct);

  for(s=1;s<=md->bd.Nn;s++){
    for(i=0;i<3;i++) r[i]=M[3*i+0]*md->bd.rn[s][0]+M[3*i+1]*md->bd.rn[s][1]+M[3*i+2]*md->bd.rn[s][2]+tv[i];
    for(i=0;i<3;i++) md->bd.rn[s][i]=r[i];
  }
}

void init_elem_const(BOUD *bd)
{
  int i,j,d,Ne,a,b;
  double rc[3][4];

  Ne=bd->Ne;

  // geometric constant
  for(i=1;i<=Ne;i++){
    for(d=0;d<3;d++)
      for(j=0;j<4;j++) rc[d][j]=bd->rn[bd->ed[i][j]][d];

    if(bd->ed[i][3]!=0){ // bi-linear element
      for(d=0;d<3;d++){
        bd->cr[i][d][0]=0.25*( rc[d][0]+rc[d][1]+rc[d][2]+rc[d][3]);
        bd->cr[i][d][1]=0.25*(-rc[d][0]+rc[d][1]+rc[d][2]-rc[d][3]);
        bd->cr[i][d][2]=0.25*(-rc[d][0]-rc[d][1]+rc[d][2]+rc[d][3]);
        bd->cr[i][d][3]=0.25*( rc[d][0]-rc[d][1]+rc[d][2]-rc[d][3]);
      }

      a=1;      b=2;
      bd->cw[i][0][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][0][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][0][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
      a=2;      b=0;
      bd->cw[i][1][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][1][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][1][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
      a=0;      b=1;
      bd->cw[i][2][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][2][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][2][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
    }
    else { // linear triangular element
      for(d=0;d<3;d++){
        bd->cr[i][d][0]=1.0/3.0*( rc[d][0]+rc[d][1]+rc[d][2]);
        bd->cr[i][d][1]=1.0/3.0*(-rc[d][0]+2.0*rc[d][1]-rc[d][2]);
        bd->cr[i][d][2]=1.0/sqrt(3.0)*( rc[d][2]-rc[d][0]);
        bd->cr[i][d][3]=0.0;
      }

      bd->cw[i][0][0]=( (rc[1][1]-rc[1][0])*(rc[2][2]-rc[2][0]) - (rc[2][1]-rc[2][0])*(rc[1][2]-rc[1][0]) );
      bd->cw[i][0][1]=0.0;
      bd->cw[i][0][2]=0.0;
      bd->cw[i][1][0]=( (rc[2][1]-rc[2][0])*(rc[0][2]-rc[0][0]) - (rc[0][1]-rc[0][0])*(rc[2][2]-rc[2][0]) );
      bd->cw[i][1][1]=0.0;
      bd->cw[i][1][2]=0.0;
      bd->cw[i][2][0]=( (rc[0][1]-rc[0][0])*(rc[1][2]-rc[1][0]) - (rc[1][1]-rc[1][0])*(rc[0][2]-rc[0][0]) );
      bd->cw[i][2][1]=0.0;
      bd->cw[i][2][2]=0.0;
    }
  }

  // element constant
  // gaussian quadrature node and weight
  bd->zt_44[0]=-P44_N;  bd->zt_44[1]= P44_N;  bd->zt_44[2]= P44_N;  bd->zt_44[3]=-P44_N;
  bd->et_44[0]=-P44_N;  bd->et_44[1]=-P44_N;  bd->et_44[2]= P44_N;  bd->et_44[3]= P44_N;
  bd->wt_44[0]= P44_W;  bd->wt_44[1]= P44_W;  bd->wt_44[2]= P44_W;  bd->wt_44[3]= P44_W;

  bd->zt_49[0]=-P49_N;   bd->zt_49[1]= P49_N;   bd->zt_49[2]= P49_N;
  bd->zt_49[3]=-P49_N;   bd->zt_49[4]= 0.0;     bd->zt_49[5]= P49_N;
  bd->zt_49[6]= 0.0;     bd->zt_49[7]=-P49_N;   bd->zt_49[8]= 0.0;
  bd->et_49[0]=-P49_N;   bd->et_49[1]=-P49_N;   bd->et_49[2]= P49_N;
  bd->et_49[3]= P49_N;   bd->et_49[4]=-P49_N;   bd->et_49[5]= 0.0;
  bd->et_49[6]= P49_N;   bd->et_49[7]= 0.0;     bd->et_49[8]= 0.0;
  bd->wt_49[0]= P49_W0;  bd->wt_49[1]= P49_W0;  bd->wt_49[2]= P49_W0;
  bd->wt_49[3]= P49_W0;  bd->wt_49[4]= P49_W1;  bd->wt_49[5]= P49_W1;
  bd->wt_49[6]= P49_W1;  bd->wt_49[7]= P49_W1;  bd->wt_49[8]= P49_W2;

  bd->zt_34[0]=-P34_N0;  bd->zt_34[1]= 2.0*P34_N0;  bd->zt_34[2]=-P34_N0;  bd->zt_34[3]= 0.0;
  bd->et_34[0]=-P34_N1;  bd->et_34[1]= 0.0;         bd->et_34[2]= P34_N1;  bd->et_34[3]= 0.0;
  bd->wt_34[0]= P34_W0;  bd->wt_34[1]= P34_W0;      bd->wt_34[2]= P34_W0;  bd->wt_34[3]=-P34_W1;

  bd->zt_37[0]=-P37_N0;  bd->zt_37[1]= 2.0*P37_N0;  bd->zt_37[2]=-P37_N0;
  bd->zt_37[3]= P37_N1;  bd->zt_37[4]=-2.0*P37_N1;  bd->zt_37[5]= P37_N1;  bd->zt_37[6]= 0.0;
  bd->et_37[0]=-P37_N2;  bd->et_37[1]= 0.0;         bd->et_37[2]= P37_N2;
  bd->et_37[3]= P37_N3;  bd->et_37[4]= 0.0;         bd->et_37[5]=-P37_N3;  bd->et_37[6]= 0.0;
  bd->wt_37[0]= P37_W0;  bd->wt_37[1]= P37_W0;      bd->wt_37[2]= P37_W0;
  bd->wt_37[3]= P37_W1;  bd->wt_37[4]= P37_W1;      bd->wt_37[5]= P37_W1;  bd->wt_37[6]= P37_W2;

  // gauss-legendre GLN point rule
  gauleg(-1.0,1.0,bd->xli,bd->wli,GLN);
  gauleg(-1.0,1.0,bd->xhi,bd->whi,GHN);
}

void malloc_sub_domain(DOMD *md)
{
  int *Nc,i,j,k;
  Nc=(int *)m_alloc2(md->MN+1,sizeof(int),"setup2.c, malloc_sub_domain(), Nc");
  for(i=0;i<=md->MN;i++) Nc[i]=0;

  for(i=1;i<=md->bd.Ne;i++){
    Nc[md->bd.md[i]]++;
    Nc[md->bd.sd[i]]++;
  }

  md->bd.sb=(SUBD *)m_alloc2(md->MN+1,sizeof(SUBD),"setup2.c, malloc_sub_domain(), md->bd.sb");
  for(i=0;i<=md->MN;i++){
    md->bd.sb[i].Ne=Nc[i];
    md->bd.sb[i].sid=(int *)m_alloc2(Nc[i]+1,sizeof(int),"setup2.c, malloc_sub_domain(), md->bd.sb[i],sid");
    md->bd.sb[i].U =(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup2.c, malloc_sub_domain(), md->bd.sb[i].U");
    md->bd.sb[i].dU=(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dU");
    md->bd.sb[i].E =(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup2.c, malloc_sub_domain(), md->bd.sb[i].E");
    md->bd.sb[i].dE=(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dE");
    md->bd.sb[i].H =(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup2.c, malloc_sub_domain(), md->bd.sb[i].H");
    md->bd.sb[i].dH=(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dH");
    for(j=0;j<=Nc[i];j++){
      md->bd.sb[i].U [j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_sub_domain(), md->bd.sb[i].U[j]");
      md->bd.sb[i].dU[j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dU[j]");
      md->bd.sb[i].E [j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_sub_domain(), md->bd.sb[i].E[j]");
      md->bd.sb[i].dE[j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dE[j]");
      md->bd.sb[i].H [j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_sub_domain(), md->bd.sb[i].H[j]");
      md->bd.sb[i].dH[j]=(double complex **)m_alloc2(4,sizeof(double complex *),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dH[j]");
      for(k=0;k<4;k++){
        md->bd.sb[i].U [j][k]=(double complex *)m_alloc2(4,sizeof(double complex),"setup2.c, malloc_sub_domain(), md->bd.sb[i].U[j][k]");
        md->bd.sb[i].dU[j][k]=(double complex *)m_alloc2(4,sizeof(double complex),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dU[j][k]");
        md->bd.sb[i].E [j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup2.c, malloc_sub_domain(), md->bd.sb[i].E[j][k]");
        md->bd.sb[i].dE[j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dE[j][k]");
        md->bd.sb[i].H [j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup2.c, malloc_sub_domain(), md->bd.sb[i].H[j][k]");
        md->bd.sb[i].dH[j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"setup2.c, malloc_sub_domain(), md->bd.sb[i].dH[j][k]");
      }
    }
  }

  free(Nc);
}

void mfree_sub_domain(DOMD *md)
{
  int i,j,k;

  for(i=0;i<=md->MN;i++){
    for(j=0;j<=md->bd.sb[i].Ne;j++){
      for(k=0;k<4;k++){
        free(md->bd.sb[i].U[j][k]); free(md->bd.sb[i].dU[j][k]);
        free(md->bd.sb[i].E[j][k]); free(md->bd.sb[i].dE[j][k]);
        free(md->bd.sb[i].H[j][k]); free(md->bd.sb[i].dH[j][k]);
      }
      free(md->bd.sb[i].U[j]);      free(md->bd.sb[i].dU[j]);
      free(md->bd.sb[i].E[j]);      free(md->bd.sb[i].dE[j]);
      free(md->bd.sb[i].H[j]);      free(md->bd.sb[i].dH[j]);
    }
    free(md->bd.sb[i].sid);
    free(md->bd.sb[i].U);    free(md->bd.sb[i].dU);
    free(md->bd.sb[i].E);    free(md->bd.sb[i].dE);
    free(md->bd.sb[i].H);    free(md->bd.sb[i].dH);
  }

  free(md->bd.sb);
}

void init_sub_domain(DOMD *md)
{
  int d,i,c;

  for(d=0;d<=md->MN;d++){
    c=1;
    for(i=1;i<=md->bd.Ne;i++){
      if(md->bd.md[i]==d){
        md->bd.sb[d].sid[c]=i;
        c++;
      }
      else if(md->bd.sd[i]==d){
        md->bd.sb[d].sid[c]=-i;
        c++;
      }
    }

  }
}


void init_boundary_data2(DOMD *md)
{
  double cr[3][4],cw[3][3],r[3],w[3];
  int i,j,l,m;

  // element node data
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
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        bil_rw_zeta_eta(r,w,md->bd.zt_44[j],md->bd.et_44[j],cr,cw);
        for(l=0;l<3;l++){
          md->bd.ren[i][j][l]=r[l];
          md->bd.wen[i][j][l]=w[l];
        }
      }
    }
  }
}

double fid_calc_solid_angle(int type,double r[4][3],int *flg)
{
  double n01[3],n12[3],n20[3],nt[3],sa,ca;
  double Omega,a0,a1,a2,D;
  int st,i;

  Omega=0.0;
  *flg=0;

  st=0;
  if(type==ELT3) for(i=0;i<3;i++)  st+=vuni_d(r[i]);
  else for(i=0;i<4;i++) st+=vuni_d(r[i]);
  if(st<0){
    *flg=-1;    return 0; // on boundary
  }

  vcrs_d(n01,r[0],r[1]);
  vcrs_d(n12,r[1],r[2]);
  vcrs_d(n20,r[2],r[0]);

  vcrs_d(nt,n01,n20);
  sa=vabs_d(nt);
  ca=-vdot_d(n01,n20);
  a0=atan2(sa,ca);

  vcrs_d(nt,n12,n01);
  sa=vabs_d(nt);
  ca=-vdot_d(n12,n01);
  a1=atan2(sa,ca);

  vcrs_d(nt,n20,n12);
  sa=vabs_d(nt);
  ca=-vdot_d(n20,n12);
  a2=atan2(sa,ca);

  D=vdot_d(r[0],n12);
  if(D>0.0) Omega+=a0+a1+a2-M_PI;
  else if(D<0.0) Omega-=a0+a1+a2-M_PI;
  else { // on boundary
    *flg=-1;    return 0;
  }

  if(ELT4==type){
    n01[0]=-n20[0];    n01[1]=-n20[1];    n01[2]=-n20[2];

    vcrs_d(n12,r[2],r[3]);
    vcrs_d(n20,r[3],r[0]);

    vcrs_d(nt,n01,n20);
    sa=vabs_d(nt);
    ca=-vdot_d(n01,n20);
    a0=atan2(sa,ca);

    vcrs_d(nt,n12,n01);
    sa=vabs_d(nt);
    ca=-vdot_d(n12,n01);
    a1=atan2(sa,ca);

    vcrs_d(nt,n20,n12);
    sa=vabs_d(nt);
    ca=-vdot_d(n20,n12);
    a2=atan2(sa,ca);

    D=vdot_d(r[0],n12);
    if(D>0.0) Omega+=a0+a1+a2-M_PI;
    else if(D<0.0) Omega-=a0+a1+a2-M_PI;
    else { // on boundary
      *flg=-1;    return 0;
    }
  }
  return Omega;
}

int domain_id(double *rt,DOMD *md)
{
  double rv[4][3],omega;
  double *Og=(double *)m_alloc2(md->MN+1,sizeof(double),"setup2.c, domain_id()");
  int i,j,k,d,flg;

  for(d=0;d<md->MN+1;d++){
    for(i=1;i<=md->bd.sb[d].Ne;i++){
      if(md->bd.sb[d].sid[i]>0){
        // read node data
        for(j=0;j<4;j++)
          for(k=0;k<3;k++) rv[j][k]=md->bd.rn[md->bd.ed[md->bd.sb[d].sid[i]][j]][k]-rt[k];
        omega=fid_calc_solid_angle(check_element_type(md->bd.sb[d].sid[i],&(md->bd)),rv,&flg);
        if(flg<0){ // on boundary
          free(Og);
          return d;
        }
        Og[d]+=omega;
        Og[md->bd.sd[md->bd.sb[d].sid[i]]]-=omega;

      } // end if
    }
    if(d==0 && fabs(Og[d])<2.0*M_PI){ // opened region
      free(Og);
      return d;
    }
    else if(Og[d]>2.0*M_PI){ // closed region
      free(Og);
      return d;
    }
  }

  free(Og);
  return -1; // error
}



