/**@file   vardata_vaccinebp.h
 * @brief  Variable data containing the ids of constraints in which the variable appears
 * @author Bahareh Eghtesadi
 *
 * This file implements the handling of the variable data which is attached to each file. See SCIP_VarData and \ref PRICER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VARDATA_VACCINEBP__
#define __SCIP_VARDATA_VACCINEBP__
#include <vector>

#include "scip/scip.h"

using namespace std;


struct SCIP_VarData
{
	vector < SCIP_Real >                           x_column;
   SCIP_Real                                       z_column;
   SCIP_Bool                                       vartype;
   int                                             stage;
   int                                             state;
}; 

/** create variable data */
extern
SCIP_RETCODE SCIPvardataCreatevaccinebp(
   SCIP*                                           scip,               /**< SCIP data structure */
   SCIP_VARDATA*                                   vardata,            /**< pointer to vardata */
   vector < SCIP_Real >                            x_column,
   SCIP_Real                                       z_column,
   SCIP_Bool                                       vartype,
   int                                             MinNN,
   int                                             I,
   double                                          K,
   int                                             stage,
   int                                             state
   );


/** get x_column */
extern
SCIP_RETCODE SCIPvardataGetxcolumn(
   SCIP_VARDATA*         vardata,
   vector< SCIP_Real >  & x            
   );

/** get y_column */
/*extern
SCIP_RETCODE SCIPvardataGetycolumn(
   SCIP_VARDATA*                             vardata,
   vector< vector< vector< SCIP_Real > > >   &y
   );*/


/** get z_column */
extern
SCIP_Real SCIPvardataGetzcolumn(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** get vartype */
extern
SCIP_Bool SCIPvardataGetvartype(
   SCIP_VARDATA*                             vardata            /**< variable data */
   );
extern
int SCIPvardatagetstage(
    SCIP_VARDATA*                            vardata            /**< variable data */
   );

/**get state*/
extern
int SCIPvardatagetstate(
    SCIP_VARDATA*                            vardata            /**< variable data */
   );



/** creates variable */
extern
SCIP_RETCODE SCIPcreateVarvaccinebp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );


#endif
