/**@file   probdata_vaccinebp.cpp
 * @brief  Problem data 
 * @author Bahareh Eghtesadi
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 *
 * @page PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the binpacking problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <string.h>
#include <vector>

#include "probdata_CapacityExp.h"
#include "vardata_CapacityExp.h"
#include "pricer_CapacityExp.h"
#include "GeDemand_CapacityExp.h"

#include "scip/scip.h"

using namespace std;
using namespace scip;


/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the binpacking, all variables which are created, and all
 * constrsaints.
 */
/*struct SCIP_ProbData
{
   vector< SCIP_VAR* >                          z_var;             /
   vector< SCIP_VAR* >                          vars;               
   vector< vector < vector <SCIP_CONS* > > >    sh_con;            
   vector<vector<SCIP_CONS* >>                  z_con;
   SCIP_CONS*                                   alpha_con;
   vector< vector<SCIP_CONS* > >                w_con;
   int                                           n_vars;          
   int                                           NScen;
   int                                           T;
   double                                        K;     
//};


/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "addedvar"
#define EVENTHDLR_DESC         "event handler for catching added variables"

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecAddedVar)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   SCIPdebugMessage("exec method of event handler for added variable to probdata\n");

   /* add new variable to probdata */
   SCIP_CALL( SCIPprobdataAddVar(scip, SCIPgetProbData(scip), SCIPeventGetVar(event), NULL) );

   return SCIP_OKAY;
}

/**@} */


/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                                         scip,               /**< SCIP data structure */
   SCIP_PROBDATA*                                probdata,  
   vector< SCIP_VAR* > &                         z_var,              /**< pointer to problem data */
   vector< SCIP_VAR* > &                         vars,               /**< all exist variables */
   //vector< vector < vector <SCIP_CONS* > > >  &  sh_con,            /**< set partitioning constraints for each job exactly one */
   //vector<vector<SCIP_CONS* > > &                z_con,
   //SCIP_CONS*                                    alpha_con,
   //vector< vector<SCIP_CONS* > >  &              w_con,
   int                                           n_vars,              /**< number of variables */
   int                                           NScen,
   int                                           T,
   double                                        K ,
   int                                           I,
   vector < vector < int > >                    jmat  
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   //SCIP_CALL( SCIPallocMemory(scip, probdata) );

   for (int i=0; i<NScen ; i++)
    {
       assert(z_var[i] != NULL);
       probdata->z_var.push_back(z_var[i]);
    }
   
   if( n_vars > 0 )
   {
      /* copy variable array */
      for (int j=0 ; j<n_vars ; j++)
       {
           assert(vars[j] != NULL);
           probdata->vars.push_back(vars[j]);
           //cout << "probdata->vars[j]" << probdata->vars[j] << endl;
           assert(probdata->vars[j] != NULL);

       }  
   }

   /* duplicate sh_con and w_con*/
  /* for(int t=2; t<T+1; t++)
  {
    vector<vector< SCIP_CONS* > > shrowrow;
    vector< SCIP_CONS*> wrow;
    for (int j=0 ; j<pow(K,t-1) ; j++)
      {
          wrow.push_back(w_con[t-2][j]);
          vector< SCIP_CONS*> shrow;
          for (int q=0 ; q<pow(K,T-t) ; q++)
            { 
              shrow.push_back(sh_con[t-2][j][q]);
            }
          shrowrow.push_back(shrow);
      }
    probdata->w_con.push_back(wrow);
    probdata->sh_con.push_back(shrowrow);
   }
   //Duplicate z_con
   for (int i=0; i<NScen; i++)
   {  
      vector< SCIP_CONS*> zrow;
      for (int j=2; j<T+1 ; j++)
       {  
         zrow.push_back(z_con[i][j-2]);
       }
      probdata->z_con.push_back(zrow); 
   } */ 
   //Duplicate jmat     
   for (int t=1; t<T+1 ; t++)
   {  
      vector< int > jmatrow;
      for (int j=0; j<pow(K,T-1) ; j++)
       {  
         jmatrow.push_back(jmat[t-1][j]);
       }
      probdata->jmat.push_back(jmatrow); 
   }       
     
   //probdata->alpha_con = alpha_con;
   /*(*probdata)->n_shcon = n_shcon;
   (*probdata)->n_zcon = n_zcon;
   (*probdata)->n_wcon = n_wcon;*/
   probdata->n_vars = n_vars;
   probdata->T = T;
   probdata->K = K;
   probdata->I = I;
   probdata->NScen = NScen;
   

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< pointer to problem data */
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* release all variables */
   for( i = 0; i < probdata->n_vars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, & probdata->vars[i]) );
   }
   
   for (int j=0 ; j< probdata->NScen ; j++)
   {
      SCIP_CALL( SCIPreleaseVar(scip, & probdata->z_var[j]) );
   }

   /* release all constraints */
   /* for(int t=2; t<(probdata->T)+1; t++)
  {
    for (int j=0 ; j<pow(probdata->K,t-1) ; j++)
      {
          for (int q=0 ; q<pow(probdata->K,probdata->T-t) ; q++)
            { 
              SCIP_CALL( SCIPreleaseCons(scip, & probdata->sh_con[t-2][j][q]) );
            }
          SCIP_CALL( SCIPreleaseCons(scip, & probdata->w_con[t-2][j]) );
      }
   }
   //Duplicate z_con
   for (int i=0; i< probdata->NScen; i++)
   {  
      for (int j=2; j< probdata->T+1 ; j++)
       {  
         SCIP_CALL( SCIPreleaseCons(scip, & probdata->z_con[i][j-2]) );
       } 
   } */      
   /* free memory of arrays */
  // SCIPfreeMemoryArray(scip, &(probdata)->vars);

   /* free probdata */
   //SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigvaccinebp)
{
   SCIPdebugMessage("free original problem data\n");

   SCIP_CALL( probdataFree(scip, *probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransvaccine)
{
   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, *targetdata, sourcedata->z_var ,  sourcedata->vars, sourcedata->n_vars, sourcedata->NScen,
                          sourcedata-> T, sourcedata-> K, sourcedata-> I, sourcedata-> jmat) );

   /* transform all constraints */
   //SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nitems, (*targetdata)->conss, (*targetdata)->conss) );

   /* transform all variables */
   //SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nvars, (*targetdata)->vars, (*targetdata)->vars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransvaccinebp)
{
   SCIPdebugMessage("free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, *probdata) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolvaccinebp)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* catch variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolvaccinebp)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* drop variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );


   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
/*SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                         scip,               
   const char*                   probname,           
   int                           T,
   double                        K,
   int                           I,
   vector<double>                purchcost,
   vector<int>                   vialsize,
   double                        Shortbeta,
   double                        Wastbeta,
   int                           Nper_to_spoil,
   int                           NScen,
   int                           Nnodes,
   vector<int>                   Demand,
   vector<double>                Dprob,
   vector<double>                Scenprob,
   vector<double>                Nodeprobe,
   vector<double>                Nodedemand,
   vector<vector<int> >          nodemat,
   vector<vector<int> >          jmat
   )
{
   

   return SCIP_OKAY;
}*/
 

/*Creats probdata*/
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                                        scip,               /**< SCIP data structure */
   SCIP_PROBDATA*                               probdata,  
   vector< SCIP_VAR* > &                         z_var,              /**< pointer to problem data */
   vector< SCIP_VAR* > &                         vars,               /**< all exist variables */
   //vector< vector < vector <SCIP_CONS* > > > &   sh_con,            /**< set partitioning constraints for each job exactly one */
   //vector<vector<SCIP_CONS* > >  &               z_con,
   //SCIP_CONS*                                   alpha_con,
   //vector< vector<SCIP_CONS* > >  &              w_con,
   int                                          n_vars,              /**< number of variables */
   int                                          NScen,
   int                                          T,
   double                                       K,
   int                                          I,
   vector < vector < int > >                    jmat
)
{
   /* create event handler if it does not exist yet */
   if( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL )
   {
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL) );
   }

   SCIP_CALL( probdataCreate(scip, probdata, z_var, vars, n_vars, NScen,T, K, I, jmat) );
   
   return SCIP_OKAY;
}
/** returns number of variables */
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata           
   )
{
   return probdata->n_vars;
}

/*returns NScen*/
int SCIPprobdataGetNScen(SCIP_PROBDATA * probdata )
{
   return probdata->NScen;
}


/** returns w variables */
SCIP_RETCODE SCIPprobdataGetVars(
   SCIP_PROBDATA*                  probdata,
   vector< SCIP_VAR*>              & wvars             
   )
{
	//cout << "wvars" << probdata->vars[0] << endl;
   for (int i =0 ; i <SCIPprobdataGetNVars(probdata) ; i++)
   {
     wvars.push_back(probdata->vars[i]);
    } 
   //cout << "wvars" << wvars[0] << endl;
   
   return SCIP_OKAY;
}

/** returns z variables */
SCIP_RETCODE SCIPprobdataGetZVars(
    SCIP_PROBDATA*                  probdata,
    vector< SCIP_VAR*>              & zvar  
)
{
  for (int i =0 ; i <SCIPprobdataGetNScen(probdata) ; i++)
   {
     zvar.push_back(probdata->z_var[i]);
    } 
   
   return SCIP_OKAY;
}

SCIP_RETCODE SCIPprobdatagetjmat (
   SCIP_PROBDATA*                probdata,
   vector < vector < int > >     & jjmat
)
{
  for (int t=1; t<SCIPprobdataGetT(probdata)+1 ; t++)
   {  
      vector< int > jmatrow;
      for (int j=0; j<pow(SCIPprobdataGetK(probdata),SCIPprobdataGetT(probdata)-1) ; j++)
       {  
         jmatrow.push_back(probdata->jmat[t-1][j]);
       }
      jjmat.push_back(jmatrow); 
      
   } 
  return SCIP_OKAY;      
}

int SCIPprobdataGetT(SCIP_PROBDATA * probdata )
{
  return probdata->T;
}

double SCIPprobdataGetK(SCIP_PROBDATA * probdata )
{
  return probdata->K;
}

int SCIPprobdataGetI(SCIP_PROBDATA * probdata )
{
  return probdata->I;
}


/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var,                /**< variables to add */
   SCIP_VAR*             z_var
   )
{
   if (z_var != NULL)
   {
     SCIP_CALL( SCIPcaptureVar(scip, z_var) );
     probdata->z_var.push_back( z_var);
   }
   if (var != NULL)
   {
     /* caputure variables */
     SCIP_CALL( SCIPcaptureVar(scip, var) );

     //cout << "wvar" << var << endl;
     (probdata->vars).push_back(var);
     //cout << "wvar" << probdata->vars[0] << endl;
     probdata->n_vars++;
   }
   SCIPdebugMessage("added variable to probdata; nvars = %d\n", probdata->n_vars);

   return SCIP_OKAY;
}
