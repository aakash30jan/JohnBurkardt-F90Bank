# include <iostream>
# include "TasmanianSparseGrid.hpp"

TasGrid::TasmanianSparseGrid tsgFCgrid;

extern "C" 
{
  extern void tsgfbridgegetpoints_( int*, int*, double* );

  extern void tsgfbridgegetneededpoints_( int*, int*, double* );

  void tsgfcmakelocalpolynomialgrid_( int* dimensions, int* outputs, 
    int* depth, int* order, int* boundary )
  {
    if ( (*boundary) == 1 )
    {
      tsgFCgrid.makeLocalPolynomialGrid( *dimensions, *outputs, *depth, 
        *order, TasGrid::rule_pwpolynomial );
    }
    else if ( (*boundary) == 0 )
    {
      tsgFCgrid.makeLocalPolynomialGrid( *dimensions, *outputs, *depth, 
        *order, TasGrid::rule_pwpolynomial0 );
    }
    else
    {
      std::cout << "ERROR: Wrong boundary type" << endl;
    }
  }

  void tsgfcgetnumpoints_( int *n )
  {
  //std::cout << " c nump = " << tsgFCgrid.getNumPoints() << endl;
    *n = tsgFCgrid.getNumPoints();
  }

  int tsgfcgetnumdimensions_( int *d )  
  {
  //std::cout << " d numd = " << tsgFCgrid.getNumDimensions() << endl;
    *d = tsgFCgrid.getNumDimensions();
  }

  int tsgfcgetnumoutputs_( int *o )
  {
  //std::cout << " d numd = " << tsgFCgrid.getNumDimensions() << endl;
    *o = tsgFCgrid.getNumOutputs();
  }

  void tsgfcgetpoints_()
  {
    double *c_points = 0;
    tsgFCgrid.getPoints( c_points );
        
  //std::cout << " c nump = " << nump << endl;
        
    int nump = tsgFCgrid.getNumPoints();
    int numd = tsgFCgrid.getNumDimensions();
        
    tsgfbridgegetpoints_( &nump, &numd, c_points );
    c_points = 0;
  }

  void tsgfcgetneededpoints_()
  {
    double *c_points = 0;
    tsgFCgrid.getPoints( c_points );
        
  //std::cout << " c nump = " << nump << endl;
   
    int nump = tsgFCgrid.getNumPoints();
    int numd = tsgFCgrid.getNumDimensions();
        
    tsgfbridgegetneededpoints_( &nump, &numd, c_points );
    c_points = 0;
  }

  void tsgfcloadneededpoints_( double *needed_points )
  {
    tsgFCgrid.loadNeededPoints( needed_points );
  }

//void c_tsgdeletegrid_()
//{
//
//}

}
