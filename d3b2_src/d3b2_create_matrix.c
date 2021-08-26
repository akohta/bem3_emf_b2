/*
 * d3b2_create_matrix.c
 *
 *  Created on: Feb 16, 2019
 *      Author: ohta
 */

#include "bem3_emf_b2.h"


int main(int argc,char **argv)
{
  DOMD md;

  read_domd2(argc,argv,&md);
  print_domd2(&md);
  //print_domd2_mksa(&md);
  initialize_domd2(&md);

  create_matrix(&md,argv[3]);

  finalize_domd2(&md);
  return 0;
}

