#include "bem3_emf_b2.h"

int main(int argc,char *argv[]) 
{
  MOBJ mo;
  FILE *fp1,*fp2;
  double complex e[3],h[3];
  double rang,dr,r[3],*ie,*ih;
  int max,i,j,type;

  mo_dat_read(argv[1],&mo); // read data file outputed by d3b2_bv_solver
  mo_print_data(&mo);       // print data
  
  max=200;
  rang=1.0*mo.md[0].mw.lambda_0;
  dr=rang*2.0/(double)(max-1);
  type=1;
  
  ie=(double *)m_alloc2(max,sizeof(double),"exampl2.c,ie");
  ih=(double *)m_alloc2(max,sizeof(double),"exampl2.c,ih");
  
  // x=0 plane
  if((fp1=fopen("Ie_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# y z electric_field_intensity");
  if((fp2=fopen("Ih_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp2,"%s\n","# y z magnetic_field_intensity");
  r[0]=0.0;
  for(i=0;i<max;i++){
    r[1]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(e,h) // omp parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      mo_EH_t(e,h,r,type,&mo); // total field 
      ie[j]=creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2]));
      ih[j]=creal(h[0]*conj(h[0]))+creal(h[1]*conj(h[1]))+creal(h[2]*conj(h[2]));
    } // end parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[1],r[2],ie[j]);
      fprintf(fp2,"%g %g %15.14e\n",r[1],r[2],ih[j]);
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
  }
  fclose(fp1);
  fclose(fp2);

  // y=0 plane
  if((fp1=fopen("Ie_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# x z electric_field_intensity");
  if((fp2=fopen("Ih_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp2,"%s\n","# x z magnetic_field_intensity");
  r[1]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(e,h) // omp parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      mo_EH_t(e,h,r,type,&mo);
      ie[j]=creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2]));
      ih[j]=creal(h[0]*conj(h[0]))+creal(h[1]*conj(h[1]))+creal(h[2]*conj(h[2]));
    }// end parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[2],ie[j]);
      fprintf(fp2,"%g %g %15.14e\n",r[0],r[2],ih[j]);
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
  }
  fclose(fp1);
  fclose(fp2);

  // z=0 plane  
  if((fp1=fopen("Ie_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  } 
  fprintf(fp1,"%s\n","# x y electric_field_intensity");
  if((fp2=fopen("Ih_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  } 
  fprintf(fp2,"%s\n","# x y magnetic_field_intensity");
  r[2]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(e,h) // omp parallel
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      mo_EH_t(e,h,r,type,&mo);
      ie[j]=creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2]));
      ih[j]=creal(h[0]*conj(h[0]))+creal(h[1]*conj(h[1]))+creal(h[2]*conj(h[2]));
    }
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[1],ie[j]);
      fprintf(fp2,"%g %g %15.14e\n",r[0],r[1],ih[j]);
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
  }
  fclose(fp1);
  fclose(fp2);
  
  printf("Intensity plot is finished\n");
  
  free(ie);
  free(ih);
  mo_finalize(&mo);
  return 0;
}
