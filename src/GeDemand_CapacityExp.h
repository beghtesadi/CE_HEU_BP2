//#ifndef __SCIP_VARDATA_VACCINEBP__
//#define __SCIP_VARDATA_VACCINEBP__
#include <vector>
#include <iostream>
#include "scip/scip.h"

using namespace std;


/*GEnerates demand and probabilities*/
void GenerateDemand (vector<int> &x , vector<double> &y, double mean, double K);


double CalProb (int n, double mean);

/*Find minimum of two integers*/
int CalMin (int  a, int b);
//#endif
