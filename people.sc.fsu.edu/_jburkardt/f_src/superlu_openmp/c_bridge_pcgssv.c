# include <stdlib.h>
# include <stdio.h>

# include "pcsp_defs.h"
# include "slu_mt_util.h"
# include "slu_mt_Cnames.h"

/******************************************************************************/

void c_bridge_pcgssv_ ( int *nprocs, int *n, int *nnz, int *nrhs, complex *values,
  int *rowind, int *colptr, complex *b, int *ldb, int *info )

/******************************************************************************/
{
  SuperMatrix A, B, L, U;
  SCformat *Lstore;
  NCformat *Ustore;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      panel_size, permc_spec, i;
  superlu_memusage_t superlu_memusage;
/* 
  Adjust to 0-based indexing 
*/
  for (i = 0; i < *nnz; ++i) --rowind[i];

  for (i = 0; i <= *n; ++i) --colptr[i];

  cCreate_CompCol_Matrix ( &A, *n, *n, *nnz, values, rowind, colptr, 
    SLU_NC, SLU_C, SLU_GE );

  cCreate_Dense_Matrix ( &B, *n, *nrhs, b, *ldb, SLU_DN, SLU_C, SLU_GE );

  if ( !(perm_r = intMalloc(*n)) ) 
  {
    SUPERLU_ABORT("Malloc fails for perm_r[].");
  }

  if ( !(perm_c = intMalloc(*n)) ) 
  {
    SUPERLU_ABORT("Malloc fails for perm_c[].");
  }
/*
  Get column permutation vector perm_c[], according to permc_spec:
  0: natural ordering 
  1: minimum degree ordering on structure of A'*A
  2: minimum degree ordering on structure of A'+A
  3: approximate minimum degree for unsymmetric matrices
*/    	
  permc_spec = 1;
  get_perm_c ( permc_spec, &A, perm_c );

  panel_size = sp_ienv(1);
    
  pcgssv ( *nprocs, &A, perm_c, perm_r, &L, &U, &B, info );
    
  if ( *info == 0 ) 
  {

	Lstore = ( SCformat * ) L.Store;
	Ustore = ( NCformat * ) U.Store;
    printf("#NZ in factor L = %d\n", Lstore->nnz);
    printf("#NZ in factor U = %d\n", Ustore->nnz);
    printf("#NZ in L+U = %d\n", Lstore->nnz + Ustore->nnz - L.ncol);
	
	superlu_dQuerySpace ( *nprocs, &L, &U, panel_size, &superlu_memusage );

	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
	       superlu_memusage.expansions);
  } 
  else
  {
	printf ( "cgssv() error returns INFO= %d\n", *info );
/* 
  Factorization completes 
*/
	if ( info <= n ) 
    { 
	  superlu_dQuerySpace(*nprocs, &L, &U, panel_size, &superlu_memusage);

	  printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		superlu_memusage.for_lu/1e6, 
        superlu_memusage.total_needed/1e6,
	    superlu_memusage.expansions);
	}
  }

  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_SuperMatrix_Store(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
/* 
  Restore to 1-based indexing 
*/
  for ( i = 0; i < *nnz; ++i ) ++rowind[i];

  for ( i = 0; i <= *n; ++i ) ++colptr[i];

  return;
}
