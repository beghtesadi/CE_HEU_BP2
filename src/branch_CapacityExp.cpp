#include <assert.h>
#include <string.h>
#include <vector>

#include "branch_CapacityExp.h"
#include "conshdlr_CapacityExp.h"
#include "probdata_CapacityExp.h"
#include "vardata_CapacityExp.h"

#include "objscip/objscip.h"

using namespace std;
using namespace scip;

/**@name Branching rule properties
 *
 * @{
 */
ObjBranchrulevaccinebp::ObjBranchrulevaccinebp(
             SCIP *             scip,
             const char * 	    name
):
ObjBranchrule( scip, name, "branchs on Original vars", 500000 , -1 , 1.0)
{}
   

ObjBranchrulevaccinebp::~ObjBranchrulevaccinebp()
{}
/**@} */

/**@name Callback methods
 *
 * @{
 */

/** branching execution method for fractional LP solutions */
SCIP_DECL_BRANCHEXECLP(ObjBranchrulevaccinebp::scip_execlp)		
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   int nlpcands;
   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);


   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   vector< SCIP_VAR* > zcand;
   vector< SCIP_VAR* > wcand;
   SCIP_Real zcandv;
   SCIP_NODE* childdown;  // Node with variable rounded down
   SCIP_NODE* childup;
   SCIP_CONS* consdown;
   SCIP_CONS* consup;

   SCIP_VARDATA* vardata;

   SCIP_CALL(SCIPprobdataGetZVars(probdata, zcand));

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), "vaccinebp_Branchrule") == 0);
   assert(result != NULL);

   //SCIPdebugMessage("start branching at node %"SCIP_LONGINT_FORMAT", depth %d\n", SCIPgetNNodes(scip), SCIPgetDepth(scip));

   Objconshdlrvaccinebp* objconsvaccine;
   const char * objconshdlr_name = "vaccinebp_Conshdlr"; 
   objconsvaccine = (Objconshdlrvaccinebp*)(SCIPfindObjConshdlr(scip, objconshdlr_name));
   
   SCIP_SOL* sol;
   sol = SCIPgetBestSol (scip);

   *result = SCIP_DIDNOTRUN;
   int nzcand = 0 ;
   for (int i=0 ; i< SCIPprobdataGetNScen(probdata) ; i++)
   {
     zcandv = SCIPvarGetLPSol(zcand[i]); //SCIPgetSolVal(scip, sol, zcand[i]); //SCIPvarGetLPSol(zcand[i]);
     cout << "i" << i << "zcand" << zcandv << endl;
     if (zcandv <= 0.000000001 || 1- zcandv <= 0.000000001)
        continue;
     else
         {
        /* create the branch-and-bound tree child nodes of the current node */
         SCIP_CALL( SCIPcreateChild(scip, &childdown, 0.0, SCIPgetLocalTransEstimate(scip)) );
         SCIP_CALL( SCIPcreateChild(scip, &childup, 0.0, SCIPgetLocalTransEstimate(scip)) );
        /* create corresponding constraints */
         SCIP_CALL(objconsvaccine->SCIPcreateconshdlr(scip, &consdown, "Down", childdown, true, zcandv , 1, SCIPprobdataGetT(probdata), i, 0, 0));
        /* add constraints to nodes */
         SCIP_CALL( SCIPaddConsNode(scip, childdown, consdown, NULL) );
         SCIP_CALL(SCIPchgVarUbNode(scip, childdown, zcand[i],0));
         //SCIP_CALL( SCIPaddConsNode(scip, childup, consup, NULL) );
         SCIP_CALL(SCIPchgVarLbNode(scip, childup,zcand[i],1));
        /* release constraints */
        SCIP_CALL( SCIPreleaseCons(scip, &consdown) );
        nzcand ++;
        *result = SCIP_BRANCHED;
        return SCIP_OKAY;
        }
    }

   cout << nlpcands << endl;
   assert(*result == SCIP_BRANCHED);
   return SCIP_OKAY;
}
