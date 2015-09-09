#include <assert.h>
#include <string.h>
#include <iostream>
#include "objscip/objscip.h"

#include "conshdlr_CapacityExp.h"
#include "conshdlr2_CapacityExp.h"
#include "probdata_CapacityExp.h"
#include "vardata_CapacityExp.h"
#include "pricer_CapacityExp.h"

using namespace scip;
using namespace std;

Objconshdlrvaccinebpsec::Objconshdlrvaccinebpsec(
		SCIP*                    scip,
		const char *             name

):
		  ObjConshdlr(scip, name, "enforces the integrality conditions",-10000 ,10, 1000, -1, -1, 1, 0, FALSE,FALSE,FALSE,
				  FALSE, SCIP_PROPTIMING_BEFORELP)
{}

Objconshdlrvaccinebpsec::~Objconshdlrvaccinebpsec()
{}

struct SCIP_ConsData
{
	int    a;
};


SCIP_RETCODE Objconshdlrvaccinebpsec::consdataFree(
		SCIP*                 scip,               /**< SCIP data structure */
		SCIP_CONSDATA**       consdata            /**< pointer to the constraint data */
)
{

	return SCIP_OKAY;
}

//Callback Methods 
SCIP_DECL_CONSDELETE(Objconshdlrvaccinebpsec::scip_delete)
{

	return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(Objconshdlrvaccinebpsec::scip_trans)
{

	return SCIP_OKAY;
}

// separation method of constraint handler for LP solution
SCIP_DECL_CONSSEPALP(Objconshdlrvaccinebpsec::scip_sepalp)
{
	*result = SCIP_FEASIBLE;
	return SCIP_OKAY;
}

// separation method of constraint handler for arbitrary primal solution
SCIP_DECL_CONSSEPASOL(Objconshdlrvaccinebpsec::scip_sepasol)
{
	*result = SCIP_FEASIBLE;
	return SCIP_OKAY;
}

// constraint enforcing method of constraint handler for LP solutions
SCIP_DECL_CONSENFOLP(Objconshdlrvaccinebpsec::scip_enfolp)
{
	SCIP_PROBDATA* probdata;
	probdata = SCIPgetProbData(scip);
	assert(probdata != NULL);
	SCIP_SOL* sol;
	sol = SCIPgetBestSol (scip);

	vector< SCIP_VAR* > wcand;
	SCIP_Real zcandv;
	SCIP_NODE* childdown;  // Node with variable rounded down
	SCIP_NODE* childup;
	SCIP_CONS* consdown;
	SCIP_CONS* consup;

	SCIP_VARDATA* vardata;

	Objconshdlrvaccinebp* objconsvaccine;
	const char * objconshdlr_name = "vaccinebp_Conshdlr";
	objconsvaccine = (Objconshdlrvaccinebp*)(SCIPfindObjConshdlr(scip, objconshdlr_name));


	vector < vector < vector < SCIP_Real > > > wcandv(SCIPprobdataGetT(probdata)-1 );
	vector < vector < int > > nwcand(SCIPprobdataGetT(probdata)-1);
	vector < vector < vector < SCIP_VARDATA* > > > xvardata(SCIPprobdataGetT(probdata)-1);

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		wcandv[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		nwcand[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		xvardata[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
	}
	SCIP_CALL( SCIPprobdataGetVars(probdata,wcand));
	for (int i=0 ; i< SCIPprobdataGetNVars(probdata) ; i++)
	{
		vardata = SCIPvarGetData(wcand[i]);
		int stage = SCIPvardatagetstage(vardata);
		int state = SCIPvardatagetstate(vardata);
		wcandv[stage-2][state].push_back(SCIPvarGetLPSol(wcand[i])); //SCIPvarGetLPSol(wcand[i])   SCIPgetSolVal(scip, sol, wcand[i])
		nwcand[stage-2][state] ++ ;
	}
	vector < vector < vector < SCIP_Real > > > zcolumn_w(SCIPprobdataGetT(probdata)-1);
	vector < vector < vector < vector < SCIP_Real > > > > xcolumn_w(SCIPprobdataGetT(probdata)-1);

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		xcolumn_w[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		zcolumn_w[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
		{
			xcolumn_w[t-2][j].resize(nwcand[t-2][j]);
			zcolumn_w[t-2][j].resize(nwcand[t-2][j]);
			cout <<  "nwcand" <<  nwcand[t-2][j] << endl;
		}

	}
	for (int i=0 ; i< SCIPprobdataGetNVars(probdata) ; i++)
	{
		vardata = SCIPvarGetData(wcand[i]);
		int stage = SCIPvardatagetstage(vardata);
		int state = SCIPvardatagetstate(vardata);
		xvardata[stage-2][state].push_back(vardata);
	}

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
		{
			for (int v=0 ; v < nwcand[t-2][j] ; v++ )
			{
				zcolumn_w[t-2][j][v] = SCIPvardataGetzcolumn(xvardata[t-2][j][v]);
				//cout << wcandv[t-2][j][v] << endl;
				cout << "zval" << zcolumn_w[t-2][j][v] << endl;
				cout << "wval" << wcandv[t-2][j][v] << endl;
			}
			SCIP_Real zcandfrac = 0;
			for (int n=0; n<nwcand[t-2][j] ; n++)
			{
				//cout << "xx" << wcandv[t-2][j][n] << endl;
				zcandfrac +=  zcolumn_w[t-2][j][n] * wcandv[t-2][j][n];
			}
			if (zcandfrac-floor(zcandfrac) > 0.000000001 && ceil(zcandfrac)-zcandfrac > 0.000000001)
			{
				/* create the branch-and-bound tree child nodes of the current node */
				SCIP_CALL( SCIPcreateChild(scip, &childdown, 0.0, SCIPgetLocalTransEstimate(scip)) );
				SCIP_CALL( SCIPcreateChild(scip, &childup, 100.0, SCIPgetLocalTransEstimate(scip)) );
				/* create corresponding constraints */
				SCIP_CALL(objconsvaccine->SCIPcreateconshdlr(scip, &consdown, "Down", childdown, true, zcandfrac, 0 ,  t, j, 0, 0));
				SCIP_CALL(objconsvaccine->SCIPcreateconshdlr(scip, &consup, "Up", childup, true, zcandfrac, 0 , t, j, 1, 0));
				/* add constraints to nodes */
				SCIP_CALL( SCIPaddConsNode(scip, childdown, consdown, NULL) );
				SCIP_CALL( SCIPaddConsNode(scip, childup, consup, NULL) );
				cout  << consdown << "---" << consup << endl;
				/* release constraints */
				SCIP_CALL( SCIPreleaseCons(scip, &consdown) );
				SCIP_CALL( SCIPreleaseCons(scip, &consup) );
				*result = SCIP_BRANCHED;
				return SCIP_OKAY;
			}
		}
	}

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
		{
			for (int v=0 ; v < nwcand[t-2][j] ; v++ )
			{
				SCIP_CALL(SCIPvardataGetxcolumn(xvardata[t-2][j][v], xcolumn_w[t-2][j][v]));
				//cout << wcandv[t-2][j][v] << endl;
				cout << "xx1-" << xcolumn_w[t-2][j][v][0] << endl;
				cout << "xx5-" << xcolumn_w[t-2][j][v][1] << endl;
				cout << "xx10-" << xcolumn_w[t-2][j][v][2] << endl;
				cout << "wval-" << wcandv[t-2][j][v] << endl;
			}
		}
	}

	int nnwcand =0;
	for (int i=0 ; i < SCIPprobdataGetI(probdata) ; i++)
	{
		for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
		{
			for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
			{
				if (nwcand[t-2][j] > 0 )
				{
					//vector< SCIP_Real> xx;
					SCIP_Real xcandfrac = 0;
					for (int n=0; n<nwcand[t-2][j] ; n++)
					{
						//xx.push_back(xcolumn_w[t-2][j][n][i]);
						//cout << "xx" << wcandv[t-2][j][n] << endl;
						xcandfrac +=  xcolumn_w[t-2][j][n][i] * wcandv[t-2][j][n];
					}
					if (xcandfrac-floor(xcandfrac) > 0.000000001 && ceil(xcandfrac)-xcandfrac > 0.000000001)
					{
						// create the branch-and-bound tree child nodes of the current node
						SCIP_CALL( SCIPcreateChild(scip, &childdown, 100.0, SCIPgetLocalTransEstimate(scip)) );
						SCIP_CALL( SCIPcreateChild(scip, &childup, 0.0, SCIPgetLocalTransEstimate(scip)) );
						// create corresponding constraints
						SCIP_CALL(objconsvaccine->SCIPcreateconshdlr(scip, &consdown, "Down", childdown, false, xcandfrac, false,  t, j, 0,i));
						SCIP_CALL(objconsvaccine->SCIPcreateconshdlr(scip, &consup, "Up", childup, false, xcandfrac, false, t, j, 1,i));
						// add constraints to nodes
						SCIP_CALL( SCIPaddConsNode(scip, childdown, consdown, NULL) );
						SCIP_CALL( SCIPaddConsNode(scip, childup, consup, NULL) );
						// release constraints
						SCIP_CALL( SCIPreleaseCons(scip, &consdown) );
						SCIP_CALL( SCIPreleaseCons(scip, &consup) );
						nnwcand ++;
						*result = SCIP_BRANCHED;
						return SCIP_OKAY;
					}

				}
			}
		}
	}

	*result = SCIP_FEASIBLE;
	return SCIP_OKAY;
}

//constraint enforcing method of constraint handler for pseudo solutions
SCIP_DECL_CONSENFOPS(Objconshdlrvaccinebpsec::scip_enfops)
{
	*result = SCIP_FEASIBLE;
	return SCIP_OKAY;
}

//feasibility check method of constraint handler for primal solutions
SCIP_DECL_CONSCHECK(Objconshdlrvaccinebpsec::scip_check)
{
	SCIP_PROBDATA* probdata;
	probdata = SCIPgetProbData(scip);
	assert(probdata != NULL);

	vector< SCIP_VAR* > wcand;
	SCIP_Real zcandv;
	SCIP_NODE* childdown;  // Node with variable rounded down
	SCIP_NODE* childup;
	SCIP_CONS* consdown;
	SCIP_CONS* consup;

	SCIP_VARDATA* vardata;

	Objconshdlrvaccinebp* objconsvaccine;
	const char * objconshdlr_name = "vaccinebp_Conshdlr";
	objconsvaccine = (Objconshdlrvaccinebp*)(SCIPfindObjConshdlr(scip, objconshdlr_name));


	vector < vector < vector < SCIP_Real > > > wcandv(SCIPprobdataGetT(probdata)-1 );
	vector < vector < int > > nwcand(SCIPprobdataGetT(probdata)-1);
	vector < vector < vector < SCIP_VARDATA* > > > xvardata(SCIPprobdataGetT(probdata)-1);

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		wcandv[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		nwcand[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		xvardata[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
	}
	SCIP_CALL( SCIPprobdataGetVars(probdata,wcand));
	for (int i=0 ; i< SCIPprobdataGetNVars(probdata) ; i++)
	{
		vardata = SCIPvarGetData(wcand[i]);
		int stage = SCIPvardatagetstage(vardata);
		int state = SCIPvardatagetstate(vardata);
		wcandv[stage-2][state].push_back( SCIPgetSolVal(scip,sol ,wcand[i])); //SCIPvarGetLPSol(wcand[i])   SCIPgetSolVal(scip, sol, wcand[i])
		nwcand[stage-2][state] ++ ;
	}
	vector < vector < vector < SCIP_Real > > > zcolumn_w(SCIPprobdataGetT(probdata)-1);
	vector < vector < vector < vector < SCIP_Real > > > > xcolumn_w(SCIPprobdataGetT(probdata)-1);

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		xcolumn_w[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		zcolumn_w[t-2].resize(pow(SCIPprobdataGetK(probdata),t-1));
		for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
		{
			xcolumn_w[t-2][j].resize(nwcand[t-2][j]);
			zcolumn_w[t-2][j].resize(nwcand[t-2][j]);
			cout <<  "nwcand" <<  nwcand[t-2][j] << endl;
		}

	}
	for (int i=0 ; i< SCIPprobdataGetNVars(probdata) ; i++)
	{
		vardata = SCIPvarGetData(wcand[i]);
		int stage = SCIPvardatagetstage(vardata);
		int state = SCIPvardatagetstate(vardata);
		xvardata[stage-2][state].push_back(vardata);
	}

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
		{
			for (int v=0 ; v < nwcand[t-2][j] ; v++ )
			{
				zcolumn_w[t-2][j][v] = SCIPvardataGetzcolumn(xvardata[t-2][j][v]);
				//cout << wcandv[t-2][j][v] << endl;
				cout << "zval" << zcolumn_w[t-2][j][v] << endl;
				cout << "wval" << wcandv[t-2][j][v] << endl;
			}
			SCIP_Real zcandfrac = 0;
			for (int n=0; n<nwcand[t-2][j] ; n++)
			{
				//cout << "xx" << wcandv[t-2][j][n] << endl;
				zcandfrac +=  zcolumn_w[t-2][j][n] * wcandv[t-2][j][n];
			}
			if (zcandfrac-floor(zcandfrac) > 0.000000001 && ceil(zcandfrac)-zcandfrac > 0.000000001)
			{

				*result = SCIP_INFEASIBLE;
				return SCIP_OKAY;
			}
		}
	}

	for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
	{
		for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
		{
			for (int v=0 ; v < nwcand[t-2][j] ; v++ )
			{
				SCIP_CALL(SCIPvardataGetxcolumn(xvardata[t-2][j][v], xcolumn_w[t-2][j][v]));
				//cout << wcandv[t-2][j][v] << endl;
				cout << "xx1-" << xcolumn_w[t-2][j][v][0] << endl;
				cout << "xx5-" << xcolumn_w[t-2][j][v][1] << endl;
				cout << "xx10-" << xcolumn_w[t-2][j][v][2] << endl;
				cout << "wval-" << wcandv[t-2][j][v] << endl;
			}
		}
	}

	int nnwcand =0;
	for (int i=0 ; i < SCIPprobdataGetI(probdata) ; i++)
	{
		for (int t=2 ; t < SCIPprobdataGetT(probdata)+1 ; t++)
		{
			for (int j=0 ; j < pow(SCIPprobdataGetK(probdata),t-1) ; j++)
			{
				if (nwcand[t-2][j] > 0 )
				{
					//vector< SCIP_Real> xx;
					SCIP_Real xcandfrac = 0;
					for (int n=0; n<nwcand[t-2][j] ; n++)
					{
						//xx.push_back(xcolumn_w[t-2][j][n][i]);
						//cout << "xx" << wcandv[t-2][j][n] << endl;
						xcandfrac +=  xcolumn_w[t-2][j][n][i] * wcandv[t-2][j][n];
					}
					if (xcandfrac-floor(xcandfrac) > 0.000000001 && ceil(xcandfrac)-xcandfrac > 0.000000001)
					{

						//nnwcand ++;
						*result = SCIP_INFEASIBLE;
						return SCIP_OKAY;
					}

				}
			}
		}
	}
	*result = SCIP_FEASIBLE;

	return SCIP_OKAY;
}

//variable rounding lock method of constraint handler
SCIP_DECL_CONSLOCK(Objconshdlrvaccinebpsec::scip_lock)
{
	return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
SCIP_DECL_CONSPROP(Objconshdlrvaccinebpsec::scip_prop)
{
	return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
//#define consLockSamediff NULL

SCIP_DECL_CONSACTIVE(Objconshdlrvaccinebpsec::scip_active)
{

	return SCIP_OKAY;

}

SCIP_DECL_CONSDEACTIVE(Objconshdlrvaccinebpsec::scip_deactive)
{
	return SCIP_OKAY;
}
