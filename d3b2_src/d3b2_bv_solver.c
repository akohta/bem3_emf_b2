#include "bem3_emf_b2.h"

int main(int argc,char **argv)
{
  MOBJ mo;

  mo_initialize(argc,argv,&mo);
  mo_print_data(&mo);
  //mo_print_data_mksa(&mo);

  mo_solve_bv(&mo);
  mo_dat_write(argv[2],&mo);
  mo_output_node_particles(argv[2],&mo);

  mo_finalize(&mo);
  return 0;
}
