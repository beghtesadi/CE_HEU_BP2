/**@file   vardata_vaccinebp.cpp
 * @brief  Variable data containing the ids of constraints in which the variable appears
 * @author Bahareh Eghtesadi
 *
 * This file implements the handling of the variable data which is attached to each file. See SCIP_VarData and \ref PRICER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/scip.h"
#include <vector>
#include <iostream>
#include "probdata_CapacityExp.h"
#include "vardata_CapacityExp.h"

using namespace std;

/** @brief Variable data which is attached to all variables.
 *
 *  This variables data is used to know in which constraints this variables appears. Therefore, the variable data
 *  contains the ids of constraints in which the variable is part of. Hence, that data give us a column view.
 */


/**@name Local methods
 *
 * @{
 */

/** create a vardata */
static
SCIP_RETCODE vardataCreate(
		SCIP*                                          scip,               /**< SCIP data structure */
		SCIP_VARDATA*                                  vardata,            /**< pointer to vardata */
		vector < SCIP_Real >                           x_column,
		SCIP_Real                                      z_column,
		SCIP_Bool                                      vartype,
		int                                            MinNN,
		int                                            I,
		double                                         K,
		int                                            stage,
		int                                            state
)
{
	for (int i=0; i< I ; i++)
	{
		vardata->x_column.push_back(x_column[i]);
	}

   vardata->z_column = z_column;
   vardata->vartype = vartype;
   vardata->stage = stage;
   vardata->state = state;


   return SCIP_OKAY;
}

/** frees user data of variable */
static
SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*        vardata             /**< vardata to delete */
   )
{
   //SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** frees user data of transformed variable (called when the transformed variable is freed) */
static
SCIP_DECL_VARDELTRANS(vardataDelTrans)
{
   SCIP_CALL( vardataDelete(scip, *vardata) );

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** create variable data */
SCIP_RETCODE SCIPvardataCreatevaccinebp(
   SCIP*                                            scip,               /**< SCIP data structure */
   SCIP_VARDATA*                                    vardata,            /**< pointer to vardata */
   vector < SCIP_Real >                             x_column,
   SCIP_Real                                        z_column,
   SCIP_Bool                                        vartype,
   int                                              MinNN,
   int                                              I,
   double                                           K,
   int                                            stage,
   int                                            state 
   )
{
   SCIP_CALL( vardataCreate(scip, vardata, x_column, z_column, vartype, MinNN, I, K, stage, state) );
   return SCIP_OKAY;
}

/** get x_column */
SCIP_RETCODE SCIPvardataGetxcolumn(
   SCIP_VARDATA*         vardata,
   vector< SCIP_Real >  & x            
   )
{
   for (int i=0; i< (int) vardata->x_column.size() ; i++)
   {
    x.push_back( vardata->x_column[i] );
	   //x[i] = vardata->x_column[i];
   }
   //cout << vardata->x_column[0] << endl;
   //cout << x[0] << endl;

   //cout << vardata->x_column[1] << endl;
   //cout << x[1] << endl;

   return SCIP_OKAY;
}
/** get y_column */
/*SCIP_RETCODE SCIPvardataGetycolumn(
   SCIP_VARDATA*                             vardata,
   vector< vector< vector< SCIP_Real > > >   &y
   )
{
   for (int i=0; i< (int) vardata->y_column.size() ; i++)
   { 
     vector<vector< SCIP_Real > > rowrow;
     for (int j=0; j< (int) vardata->y_column[i].size(); j++)
      {
        vector< SCIP_Real > row;
        for (int jj=0 ; jj< (int) vardata->y_column[i][j].size(); j++)
         {
            row.push_back(vardata->y_column[i][j][jj]);
          }
        rowrow.push_back(row);
      }
     y.push_back(rowrow);
   }//for loop
   return SCIP_OKAY;
}*/
/** get z_column */
SCIP_Real SCIPvardataGetzcolumn(
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   return vardata->z_column;
}

SCIP_Bool SCIPvardataGetvartype(
   SCIP_VARDATA*                             vardata            /**< variable data */
   )
{ 
  return vardata->vartype;
}

/**get stage*/
int SCIPvardatagetstage(
    SCIP_VARDATA*                            vardata            /**< variable data */
   )
{ 
  return vardata->stage;
}
/**get state*/
int SCIPvardatagetstate(
    SCIP_VARDATA*                            vardata            /**< variable data */
   )
{ 
  return vardata->state;
}


/** creates variable */
SCIP_RETCODE SCIPcreateVarvaccinebp(
   SCIP*                                      scip,               /**< SCIP data structure */
   SCIP_VAR**                                 var,                /**< pointer to variable object */
   const char*                                name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real                                  obj,                /**< objective function value */
   SCIP_Bool                                  initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool                                  removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*                              vardata             /**< user data for this specific variable */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   /* create a basic variable object */
   SCIP_CALL( SCIPcreateVarBasic(scip, var, name, 0.0, SCIPinfinity(scip), obj, SCIP_VARTYPE_CONTINUOUS) );
   assert(*var != NULL);

   /* set callback functions */
   SCIPvarSetData(*var, vardata);
   SCIPvarSetDeltransData(*var, vardataDelTrans);

   /* set initial and removable flag */
   SCIP_CALL( SCIPvarSetInitial(*var, initial) );
   SCIP_CALL( SCIPvarSetRemovable(*var, removable) );

   SCIPvarMarkDeletable(*var);

   SCIPdebug( SCIPprintVar(scip, *var, NULL) );

   return SCIP_OKAY;
}


/** prints vardata to file stream */
//Deleted
/**@} */
