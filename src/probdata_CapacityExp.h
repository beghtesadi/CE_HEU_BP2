/**@file   probdata_vaccinebp.h
 * @brief  Problem data 
 * @author Bahareh Eghtesadi
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_VACCINEBP__
#define __SCIP_PROBDATA_VACCINEBP__

#include <vector>

#include "scip/scip.h"
using namespace std;


struct SCIP_ProbData
{
   vector< SCIP_VAR* >                          z_var;             /**< pointer to problem data */
   vector< SCIP_VAR* >                          vars;               /**< all exist variables */
   //vector< vector < vector <SCIP_CONS* > > >    sh_con;            /**< set partitioning constraints for each job exactly one */
   //vector<vector<SCIP_CONS* > >                  z_con;
   //SCIP_CONS*                                   alpha_con;
   //vector< vector<SCIP_CONS* > >                w_con;
   int                                           n_vars;           /**< number of variables */
   int                                           NScen;
   int                                           T;
   double                                        K;
   int                                           I;
   vector < vector < int > >                    jmat;
   //int                                           n_shcon;          /**< number of items */ 
   //int                                           n_zcon;
   //int                                           n_alphacon;
   //int                                           n_wcon;     
};


/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                                        scip,               /**< SCIP data structure */
   SCIP_PROBDATA*                              probdata,  
   vector< SCIP_VAR* > &                         z_var,              /**< pointer to problem data */
   vector< SCIP_VAR* > &                         vars,               /**< all exist variables */
   //vector< vector < vector <SCIP_CONS* > > > &    sh_con,            /**< set partitioning constraints for each job exactly one */
   //vector<vector<SCIP_CONS* > > &                z_con,
   //SCIP_CONS*                                   alpha_con,
   //vector< vector<SCIP_CONS* > > &               w_con,
   int                                          n_vars,              /**< number of variables */
   int                                          NScen,
   int                                          T,
   double                                       K,
   int                                          I,
   vector < vector < int > >                    jmat
);


/** returns sh_con  */
/*extern
SCIP_CONS* SCIPprobdataGetshcon(
   SCIP_PROBDATA*        probdata            
   );

/** returns z_con  */
/*extern
SCIP_CONS* SCIPprobdataGetZcon(
   SCIP_PROBDATA*        probdata            
   );

/** returns zw_con */
/*extern
SCIP_CONS* SCIPprobdataGetwcon(
   SCIP_PROBDATA*        probdata            
   );

/** returns alpha_con*/
/*extern
SCIP_CONS* SCIPprobdataGetalpha(
   SCIP_PROBDATA*        probdata            
   );*/

/** returns number of variables */
extern
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );
extern
int SCIPprobdataGetNScen(SCIP_PROBDATA * probdata );

/** returns w variables */
extern
SCIP_RETCODE SCIPprobdataGetVars(
   SCIP_PROBDATA*                  probdata,
   vector< SCIP_VAR*>              & wvars             
   );

/** returns z variables */
extern
SCIP_RETCODE SCIPprobdataGetZVars(
    SCIP_PROBDATA*                  probdata,
    vector< SCIP_VAR*>              & zvar                
   );

extern
SCIP_RETCODE SCIPprobdatagetjmat (
   SCIP_PROBDATA*                probdata,
   vector < vector < int > >     & jjmat
);
extern
int SCIPprobdataGetT(SCIP_PROBDATA * probdata );

extern
double SCIPprobdataGetK(SCIP_PROBDATA * probdata );

extern
int SCIPprobdataGetI(SCIP_PROBDATA * probdata );



/** adds given variable to the problem data */
extern
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var,                /**< variables to add */
   SCIP_VAR*             z_var
   );

#endif
