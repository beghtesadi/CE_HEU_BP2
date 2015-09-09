#include "pricer_CapacityExp.h"
#include "probdata_CapacityExp.h"
#include "vardata_CapacityExp.h"
#include "conshdlr_CapacityExp.h"
#include "branch_CapacityExp.h"
#include "GeDemand_CapacityExp.h"
#include "scip/scipdefplugins.h"


#include <math.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "scip/cons_linear.h"


using namespace std;
using namespace scip;

#define PRICER_NAME            "vaccinebp"
#define PRICER_DESC            "variable pricer template"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */


/** Constructs the pricer object with the data needed
 *
 *  An alternative is to have a problem data class which allows to access the data.
 */
ObjPricervaccinebp::ObjPricervaccinebp(
   SCIP*                                                                scip, 
   const char*                                                         p_name,
   int                                                                  p_Nper_to_spoil,
   double                                                               p_Wastbeta,
   int                                                                  p_T,
   int                                                                  p_Nnodes,
   double                                                               p_K,
   int                                                                  p_I,
   vector<double>                                                     & p_Nodeprobe,
   vector<double>                                                     & p_purchcost,
   vector < int >                                                     & p_vialsize,
   vector < vector < int > >                                          & p_jmat,
   vector< vector<int > >                                             & p_nodemat,
   vector< SCIP_VAR* >                                                & p_z_var,
   vector< SCIP_VAR* >                                                & p_wj_var,
   vector< vector < vector <SCIP_CONS* > > >                          & p_sh_con,
   vector<vector<SCIP_CONS* > >                                       & p_z_con,
   SCIP_CONS*                                                         & p_alpha_con,
   vector< vector<SCIP_CONS* > >                                      & p_w_con,
   int                                                                  p_n_vars,
   bool                                                                 p_result
 ):     
   ObjPricer(scip, p_name, "Finds columns with negative reduced cost.", 0, TRUE),
   _Nper_to_spoil(p_Nper_to_spoil),
   _Wastbeta(p_Wastbeta),
   _T(p_T),
   _Nnodes(p_Nnodes),
   _K(p_K),
   _I(p_I),
   _Nodeprobe(p_Nodeprobe),
   _purchcost(p_purchcost),
   _vialsize(p_vialsize),
   _jmat(p_jmat),
   _nodemat(p_nodemat),
   _z_var(p_z_var),
   _wj_var(p_wj_var),
   _sh_con(p_sh_con),
   _z_con(p_z_con),
   _alpha_con(p_alpha_con),
   _w_con(p_w_con),
   _n_vars(p_n_vars),
   _result(p_result)
{}


/** Destructs the pricer object. */
ObjPricervaccinebp::~ObjPricervaccinebp()
{}

SCIP_DECL_PRICERINIT(ObjPricervaccinebp::scip_init)
{
   for (int i=0; i<pow(_K,_T-1) ; i++)
     {  //gets the transformed var and z_con
      SCIP_CALL( SCIPgetTransformedVar(scip, _z_var[i], &_z_var[i]) );
      for (int j=2; j<_T+1 ; j++)
        {
          SCIP_CALL( SCIPgetTransformedCons(scip, _z_con[i][j-2], &_z_con[i][j-2]) );
        }
     }
   cout << "n_vars"<< _n_vars << endl;
   for(int i=0; i<_n_vars ; i++)
   {
    SCIP_CALL( SCIPgetTransformedVar(scip, _wj_var[i], &_wj_var[i]) ); 
   }
  
 
   for(int t=2; t<_T+1; t++)
  {
    for (int j=0 ; j<pow(_K,t-1) ; j++)
      {
        SCIP_CALL( SCIPgetTransformedCons(scip, _w_con[t-2][j], &_w_con[t-2][j]));
        for (int q=0 ; q<pow(_K,_T-t) ; q++)
            {  
              SCIP_CALL( SCIPgetTransformedCons(scip, _sh_con[t-2][j][q], &_sh_con[t-2][j][q]) );
            }
      }
   }
   //gets transformed for alpha_con
   SCIP_CALL( SCIPgetTransformedCons(scip, _alpha_con, &_alpha_con));

   return SCIP_OKAY;
}


/** perform pricing*/
SCIP_RETCODE ObjPricervaccinebp::pricing(SCIP* scip)              
{ //Get the duals of s_sh_con
  vector< vector < vector <SCIP_Real > > > pi(_T-1);
  for(int t=2; t<_T+1; t++)
  {
    pi[t-2].resize(pow(_K,t-1)); 
    for (int j=0 ; j<pow(_K,t-1) ; j++)
      {
         pi[t-2][j].resize(pow(_K,_T-t));
         for (int q=0 ; q<pow(_K,_T-t) ; q++)
            {
              pi[t-2][j][q] = SCIPgetDualsolLinear(scip, _sh_con[t-2][j][q]);
              cout << "sh duals" << pi[t-2][j][q] << endl;
              /*if (pi[t-2][j][q] < 0)
            	  pi[t-2][j][q] = 0;*/
             }
       }
   }

   //Get the duals of 5c
   vector<vector<SCIP_Real> > gamma(pow(_K,_T-1)); 
   
   for (int i=0; i<pow(_K,_T-1); i++)
   {  
      gamma[i].resize(_T-1);
      for (int j=2; j<_T+1 ; j++)
       {
        gamma[i][j-2] = SCIPgetDualsolLinear(scip, _z_con[i][j-2]);
         cout  << "gamma" << gamma[i][j-2] << endl;
        }
   }
   //Gets the dual of w_con
   vector< vector< SCIP_Real> > mu(_T-1);
   vector< vector< SCIP_Real> > redcost(_T-1);
   for (int t=2; t< _T+1 ; t++)
   {
    mu[t-2].resize(pow(_K,t-1));
    redcost[t-2].resize(pow(_K,t-1));
    for (int j=0; j < pow (_K, t-1); j++)
     {
       mu[t-2][j]= SCIPgetDualsolLinear(scip, _w_con[t-2][j]);
       cout << "mu" << mu[t-2][j] << endl;
     }
    }


   vector< vector < vector < SCIP_Real > > > x_newcolumn (_I);
   vector< vector < vector < vector< vector< SCIP_Real > > > > > y_newcolumn (_I); 
   vector< vector < SCIP_Real > > z_newcolumn (_T-1);
   vector< vector < SCIP_Real > > objval (_T-1);
   vector< vector < SCIP_VARDATA* > > vardata (_T-1);


   SCIP_VARDATA* newcolumn_vardata_ptr;
    
  for (int i=0; i<_I ; i++)
  {
    x_newcolumn[i].resize(_T-1);
    y_newcolumn[i].resize(_T-1);
    for (int t=2; t<_T+1; t++)
    {
     x_newcolumn[i][t-2].resize(pow(_K, t-1));
     y_newcolumn[i][t-2].resize(pow(_K, t-1));
     for (int j=0; j<pow(_K, t-1); j++)
     {
       x_newcolumn[i][t-2][j] = 0;
       int YNN;
       YNN = CalMin(_T-t,_Nper_to_spoil) ;
       y_newcolumn[i][t-2][j].resize(YNN+1);
       for (int kk=0; kk < YNN+1 ; kk++)
        {
          y_newcolumn[i][t-2][j][kk].resize(pow(_K,kk));
          for (int q=0 ; q<pow(_K,kk) ; q++)
          { 
            y_newcolumn[i][t-2][j][kk][q] = 0;
          }
        }
      } 
     }//t 
    }// y[i][t][j][kk][jj]
   for (int t=2; t<_T+1; t++)
  {
    z_newcolumn[t-2].resize(pow(_K, t-1));
    objval[t-2].resize(pow(_K, t-1));
    vardata[t-2].resize(pow(_K, t-1));
    for (int j =0; j<pow(_K, t-1) ; j++)
     {
       newcolumn_vardata_ptr = new SCIP_VARDATA;
       z_newcolumn[t-2][j] = 0;
       vardata[t-2][j] = newcolumn_vardata_ptr;
      }
   }

   find_new_column (scip, pi , gamma,  x_newcolumn, y_newcolumn, z_newcolumn, objval); //for compiling
   vector< SCIP_Real > xcolumn (_I);
   vector< vector< vector < SCIP_Real > > > ycolumn (_I);
   int n_added_var =0;
   for (int t=2; t<_T+1; t++)
   {
    for (int j=0 ; j<pow(_K,t-1); j++)
     {
      int NN= CalMin(_T-t, _Nper_to_spoil);
      for (int i=0 ; i<_I; i++){
       ycolumn[i].resize(NN+1);
       xcolumn[i] = x_newcolumn[i][t-2][j];
       for(int kk=0 ; kk<NN+1; kk++){
        ycolumn[i][kk].resize(pow(_K, kk));
        for(int jj=0 ;jj<pow(_K, kk); jj++){
          ycolumn[i][kk][jj] = y_newcolumn[i][t-2][j][kk][jj];
          }
         }
        }
      redcost[t-2][j] = objval[t-2][j] -mu[t-2][j];
      cout << "objval" << objval[t-2][j] << endl;
      //see if reduced cost is negative
      
      if ( SCIPisNegative(scip, redcost[t-2][j]) )
      {
    	  //cout << "xcolumn" << xcolumn[1] << endl;
    	  //cout<< "x_newcolumn[i][t-2][j]" << x_newcolumn[1][t-2][j] << endl;
          add_newcolumn_variable(scip, xcolumn, ycolumn, z_newcolumn[t-2][j], vardata[t-2][j], t, j);
          SCIP_CALL( SCIPvardataCreatevaccinebp(scip, vardata[t-2][j], xcolumn, z_newcolumn[t-2][j], false ,NN, _I, _K, t, j ));
          n_added_var ++;
       }
      }
    }
   if (n_added_var > 0 )
   {
	  return SCIP_OKAY;
   }


   SCIP_CALL( SCIPwriteTransProblem(scip, "CapacityExp.lp", "lp", FALSE) );


   return SCIP_OKAY;
}

SCIP_DECL_PRICERREDCOST(ObjPricervaccinebp::scip_redcost)
{
   SCIPdebugMessage("call scip_redcost ...\n");

   // set result pointer, see above 
   *result = SCIP_SUCCESS;
   _result=0;

   // call pricing routine 
   SCIP_CALL( pricing(scip) );

   if (_result)
	   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
} 


/** add new variable to problem */
SCIP_RETCODE ObjPricervaccinebp::add_newcolumn_variable(
   SCIP*                                                           scip,                                  
   vector< SCIP_Real >                                           & x_newcolumn,
   vector< vector< vector < SCIP_Real > > >                      & y_newcolumn, 
   SCIP_Real                                                      z_newcolumn,
   SCIP_VARDATA*                                                  vardata,
   int                                                             t,
   int                                                             j 
    )       
{
  SCIP_Real wobjcoef =0 ;
  char var_name[255];
  SCIP_VAR * w_var;
  int lenght = pow(_K,_T-1)/pow(_K, t-1);
  for (int i=0 ; i<_I ; i++)
      {
	    //cout << x_newcolumn[i] << endl;
        wobjcoef = wobjcoef +  x_newcolumn[i] * _purchcost[i] * _Nodeprobe[_nodemat[t-1][j*lenght]-2];
      }
  //cout << "wobjcoef" << wobjcoef << endl;
  SCIPsnprintf(var_name, 255, "Alpha%d_%d", t, j);
  SCIP_CALL( SCIPcreateVarvaccinebp(scip, &w_var, var_name, wobjcoef, FALSE, TRUE, vardata) );
  SCIPdebugMessage("new variable <%s>\n", var_name);

  // add new variable to the list of variables to price into LP (score: leave 1 here) 
   SCIP_CALL( SCIPaddPricedVar(scip, w_var, 1.0) );

  SCIP_PROBDATA* probdata;
  probdata = SCIPgetProbData(scip);
  SCIP_CALL(SCIPprobdataAddVar( scip, probdata, w_var, NULL ));
  
  
  // Add the new variable to the constraint 5b
  
  for (int tt=t; tt< t+CalMin(_T-t,_Nper_to_spoil)+1 ; tt++)
  {
    int lenn= pow(_K,tt-1)/pow(_K,t-1);
    for (int jj=0 ; jj <lenn ; jj++)
    {
      for (int q=0 ; q< pow(_K, _T-tt) ; q++)
      {
        SCIP_Real coefsh = 0;
        for (int i=0 ; i<_I; i++)
        {
          coefsh = coefsh + y_newcolumn [i][tt-t][jj]; 
        }
        //cout << "coefsh" << coefsh << endl;
        //cout << _sh_con[tt-2][j*lenn+jj][q] << endl;
        SCIP_CALL( SCIPaddCoefLinear(scip, _sh_con[tt-2][j*lenn+jj][q], w_var, coefsh ) );
      }
    }
  }

   // Add to constraint 5c
   for (int jj=0 ; jj< pow(_K,_T-t); jj++)
   {
     SCIP_CALL( SCIPaddCoefLinear(scip, _z_con[j*pow(_K,_T-t)+jj][t-2], w_var, -1*z_newcolumn ) );
   }

   //Add to constraint 5e
   SCIP_CALL( SCIPaddCoefLinear(scip, _w_con[t-2][j], w_var, 1 ) );

  //cleanup 

   SCIP_CALL( SCIPreleaseVar(scip, &w_var) );

   return SCIP_OKAY;

}

//Finds new columns
SCIP_RETCODE ObjPricervaccinebp::find_new_column(
   SCIP*                                      scip,
   vector< vector < vector <SCIP_Real > > >  & pi,      
   vector< vector < SCIP_Real > > & gamma,
   vector< vector <vector< SCIP_Real > > >         & x_newcolumn,
   vector< vector <vector< vector< vector< SCIP_Real > > > > > & y_newcolumn, 
   vector< vector < SCIP_Real > > & z_newcolumn,
   vector< vector < SCIP_Real > > & objval
 )
{
	vector < vector < SCIP_Real > > _zvar_lb (_T-1);
	vector < vector < SCIP_Real > > _zvar_ub (_T-1);

	vector < vector < vector < SCIP_Real > > > _xvar_lb(_I);
	vector < vector < vector < SCIP_Real > > > _xvar_ub(_I);

	vector < SCIP_Real > x_ub ={50, 10, 5};
	//vector < SCIP_Real > x_ub ={ 50, 10, 5};
	//vector < SCIP_Real > x_ub ={10};

	int p_NN;
	for (int i =0 ; i < _I ; i++)
	   {
	     _xvar_lb[i].resize(_T-1);
	     _xvar_ub[i].resize(_T-1);
	     //_yvar_lb[i].resize(_T-1);
	     //_yvar_ub[i].resize(_T-1);
	     for (int t=2 ; t<_T+1 ; t++)
	     {
	       p_NN = CalMin(_T-t, _Nper_to_spoil);
	       _xvar_lb[i][t-2].resize(pow(_K,t-1));
	       _xvar_ub[i][t-2].resize(pow(_K,t-1));
	       for (int j=0 ; j< pow(_K,t-1) ; j++)
	       {
	         _xvar_lb[i][t-2][j] = 0;
	         _xvar_ub[i][t-2][j] = x_ub[i];
	        }
	      }
	    }

	//Initialize _zvar,
	 for (int t=2 ; t<_T+1 ; t++)
	   {
	     p_NN = CalMin(_T-t, _Nper_to_spoil);
	     _zvar_lb[t-2].resize(pow(_K,t-1));
	     _zvar_ub[t-2].resize(pow(_K,t-1));
	       for (int j=0 ; j< pow(_K,t-1) ; j++)
	       {
	         _zvar_lb[t-2][j] = 0;
	         _zvar_ub[t-2][j] = 1;
	        }
	      }
  SCIP_CONSHDLR* myconshdlr;
  myconshdlr = SCIPfindConshdlr(scip, "vaccinebp_Conshdlr");
  SCIP_CALL(addBranchingDecisionConss(scip, myconshdlr, _xvar_lb, _xvar_ub, _zvar_lb, _zvar_ub));
  if (_result)
	  return SCIP_OKAY;
  vector < vector < vector <SCIP_Real > > > xvarlb(_T-1);
  vector < vector < vector <SCIP_Real > > > xvarub(_T-1);
  for (int t=2; t<_T+1 ; t++)
  {
	  xvarlb[t-2].resize(pow(_K,t-1));
	  xvarub[t-2].resize(pow(_K,t-1));
	  for (int j=0 ; j< pow(_K,t-1) ; j++)
	  {
		  for (int i=0; i<_I ; i++)
		  {
			  xvarlb[t-2][j].push_back(_xvar_lb[i][t-2][j]);
			  cout << "zlb" << _xvar_lb[i][t-2][j] << endl;
			  xvarub[t-2][j].push_back(_xvar_ub[i][t-2][j]);
			  cout << "xub" << _xvar_ub[i][t-2][j] << endl;
		  }
	  }
  }
  for (int t=2 ; t<_T+1; t++)
  {   
    int lenght = pow(_K,_T-1)/pow(_K,t-1); 
    for (int j=0; j<pow(_K,t-1) ; j++)
    {
    	int NNN = CalMin(_T-t, _Nper_to_spoil);
    	SCIP_Real y_ub =50;
    	//vector < SCIP_Real > x_ub = { 50, 10, 5};
    	SCIP_Real zcoeff = 0;
    	for (int q=0; q< pow(_K,_T-t); q++)
    	{
    		zcoeff += gamma[j*lenght+q][t-2];
    	}
    	vector <vector < int > > KK(_I);
    	vector <vector < int > > JJ(_I);
    	//vector <vector < int > > KK_beta(_I);
    	//vector <vector < int > > JJ_beta(_I);
    	//vector <bool > allzero(_I , false);
    	vector < vector < vector < SCIP_Real > > > ycoeff (_I);
    	vector < SCIP_Real > zval (_I);
    	SCIP_CALL(findyandz(pi, gamma,  ycoeff, KK,  JJ, zval, xvarlb[t-2][j], xvarub[t-2][j], _zvar_lb[t-2][j], _zvar_ub[t-2][j], NNN, j, t, y_ub));
    	//cout << "findzyz isdone" << endl;

    	//Set x and y to their lb:
    	for (int i =0; i<_I ; i++)
    	{
    		for (int kk=0 ; kk< (int)KK[i].size() ; kk++)
    		{
    			y_newcolumn[i][t-2][j][KK[i][kk]][JJ[i][kk]] = _xvar_lb[i][t-2][j]*_vialsize[i];
    		}
    		x_newcolumn[i][t-2][j] = _xvar_lb[i][t-2][j];
    	}

    	cout << "zlb" << _zvar_lb[t-2][j] << endl;
    	cout << "zub" << _zvar_ub[t-2][j] << endl;

    	for (int i =0; i<_I ; i++)
    	{
    		cout << "zval" << zval[i] << endl;
    		SCIP_Real a = 0;
    		//SCIP_Real b = 0;
    		for (int kk=0 ; kk< (int)KK[i].size() ; kk++)
    		{
    			a += ycoeff[i][KK[i][kk]][JJ[i][kk]];
    		}
    		if (_xvar_lb[i][t-2][j] == _xvar_ub[i][t-2][j])   //if x is fixed
    		{
    			x_newcolumn[i][t-2][j] = _xvar_lb[i][t-2][j];
    			for (int kk=0 ; kk< (int)KK[i].size() ; kk++)
    			{
    				y_newcolumn[i][t-2][j][KK[i][kk]][JJ[i][kk]] = _xvar_lb[i][t-2][j]* _vialsize[i] ;
    			}
    			if (_xvar_lb[i][t-2][j] ==0 )
    				zval[i] = _zvar_lb[t-2][j];
    		}
    		else if ( _Nodeprobe[_nodemat[t-1][j*lenght]-2] * _purchcost[i] * _xvar_ub[i][t-2][j] <= (a* _xvar_ub[i][t-2][j]*_vialsize[i])) // if x is not fixed, try the upper bound
    		{
    			for (int kk=0 ; kk< (int)KK[i].size() ; kk++)
    			{
    				y_newcolumn[i][t-2][j][KK[i][kk]][JJ[i][kk]] = _xvar_ub[i][t-2][j]*_vialsize[i];
    			}
    			/*for (int kk=0 ; kk< (int)KK_beta[i].size() ; kk++)
    			{
    				y_newcolumn[i][t-2][j][KK_beta[i][kk]][JJ_beta[i][kk]] = y_ub-_Wastbeta;
    			}*/
    			x_newcolumn[i][t-2][j] = _xvar_ub[i][t-2][j];
    		}
    		else if (_xvar_lb[i][t-2][j] == 0)
    		{
    			zval[i] = _zvar_lb[t-2][j]; //Since we are not going to open anything
    		}
    	}
    	cout << z_newcolumn[t-2][j] << endl;

    	//Dteremin znewcolumn
    	for (int ii=0 ; ii<_I; ii++)
    	{
    		if (z_newcolumn[t-2][j] < zval[ii])
    			z_newcolumn[t-2][j] = zval[ii];
    	}

    	cout << z_newcolumn[t-2][j] << endl;
    	for (int i =0 ; i< _I ; i++)
    	{
    		objval[t-2][j] += _Nodeprobe[_nodemat[t-1][j*lenght]-2] * _purchcost[i] * x_newcolumn[i][t-2][j];
    		cout << x_newcolumn[i][t-2][j] << endl;
    		//cout << "xcoeff" << _Nodeprobe[_nodemat[t-1][j*lenght]-2] * _purchcost[i] << endl;
    		cout << "objchg" << _Nodeprobe[_nodemat[t-1][j*lenght]-2] * _purchcost[i] * x_newcolumn[i][t-2][j] << endl;
    		for (int kk=0; kk <NNN+1; kk++)
    		{
    			for (int jj=0; jj < pow(_K,kk) ; jj++)
    			{
    				objval[t-2][j] -= ycoeff[i][kk][jj] * y_newcolumn[i][t-2][j][kk][jj];
    			}
    		}
    		cout << objval[t-2][j] << endl;
    	}
    	objval[t-2][j] += zcoeff * z_newcolumn[t-2][j];
  }//j

 }//t

return SCIP_OKAY;
} 


/** add branching decisions constraints to the sub SCIP */
SCIP_RETCODE ObjPricervaccinebp::addBranchingDecisionConss(
   SCIP*                                                                  scip,
   SCIP_CONSHDLR*                                                         conshdlr,
   vector < vector < vector < SCIP_Real > > >                             & _xvar_lb,
   vector < vector < vector < SCIP_Real > > >                             & _xvar_ub,
   vector < vector < SCIP_Real > >                                        & _zvar_lb,
   vector < vector < SCIP_Real > >                                        & _zvar_ub
   )
{
   SCIP_CONS** conss;
   SCIP_CONS* cons;
   int nconss;
   int stage;
   int state;
   int nodetype;
   SCIP_Bool candtype;
   SCIP_Real candfrac;
   bool ztype;
   int  vialtype;

   int c;

   assert( scip != NULL );
   //assert( subscip != NULL );
   cout << "conshdlr" << conshdlr << endl;
   assert( conshdlr != NULL );

   // collect all branching decision constraints
   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);
   Objconshdlrvaccinebp* objconsvaccine;
   const char * objconshdlr_name = "vaccinebp_Conshdlr"; 
   objconsvaccine = (Objconshdlrvaccinebp*)(SCIPfindObjConshdlr(scip, objconshdlr_name));
    //loop over all branching decision constraints and apply the branching decision if the corresponding constraint is
    // active
   
   cout << "nconss" << nconss << endl;
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];

       // ignore constraints which are not active since these are not laying on the current active path of the search
       // tree
   
      if( !SCIPconsIsActive(cons) )
         continue;

      // collect the two item ids and the branching type (SAME or DIFFER) on which the constraint branched 
      stage = objconsvaccine->SCIPgetstage(scip, cons);
      state = objconsvaccine->SCIPgetstate(scip, cons);
      ztype = objconsvaccine->SCIPgetztype(scip , cons);
      candtype = objconsvaccine->SCIPgetcandtype(scip, cons);
      nodetype = objconsvaccine->SCIPgetnodetype(scip, cons);
      candfrac = objconsvaccine->SCIPgetcandfrac(scip , cons);
      vialtype = objconsvaccine->SCIPgetvialtype(scip , cons);


      //SCIPdebugMessage("create varbound for %s(%d,%d)\n", type == SAME ? "same" : "diff",
         //SCIPprobdataGetIds(SCIPgetProbData(scip))[id1], SCIPprobdataGetIds(SCIPgetProbData(scip))[id2]);

       // depending on the branching type select the correct left and right hand side for the linear constraint which
       // enforces this branching decision in the pricing problem MIP
       
 
      if (candtype == true && nodetype == 0)
      {
    	  if (ztype)
    	  {
    		  for (int t=2; t < _T+1 ; t++)
    		  {
    			  _zvar_ub[t-2][_jmat[t-1][state]] = 0.0 ;
    			  //assert(_zvar_ub[t-2][_jmat[t-1][state]] > _zvar_lb[t-2][_jmat[t-1][state]]);
    			  if (_zvar_ub[t-2][_jmat[t-1][state]] < _zvar_lb[t-2][_jmat[t-1][state]])
    				  _result = 1;
    		  }
    	  }
    	  else
    	  {
    		  _zvar_ub[stage-2][state] = 0.0 ;
    		  //assert(_zvar_ub[stage-2][state] > _zvar_lb[stage-2][state]);
    		  if (_zvar_ub[stage-2][state] < _zvar_lb[stage-2][state])
    			  _result = 1;
    	  }
      }
      else if (candtype == true && nodetype == 1)
      {
          _zvar_lb[stage-2][state] = 1.0 ;
          //assert(_zvar_lb[stage-2][state] < _zvar_ub[stage-2][state]);
          if (_zvar_lb[stage-2][state] > _zvar_ub[stage-2][state])
        	  _result = 1;

      }
      else if (candtype == false && nodetype == 0)
      {
          _xvar_ub[vialtype][stage-2][state] = floor(candfrac) ;
      }
      else if (candtype == false && nodetype == 1)
      {
          _xvar_lb[vialtype][stage-2][state] = ceil(candfrac) ;
      }
      else
      {
         SCIPerrorMessage("unknow constraint type <%d>\n, type");
         return SCIP_INVALIDDATA;
      }

      cout << _xvar_lb[0][0][0] << endl;

     
   }

   return SCIP_OKAY;
}
SCIP_RETCODE ObjPricervaccinebp::findyandz(
		vector< vector < vector < SCIP_Real > > >    pi,      /**< dual variable value*/
		vector< vector < SCIP_Real > >               gamma,
		vector < vector < vector < SCIP_Real > > >  & ycoeff,
		vector < vector < int > >                   & KK,
		vector < vector < int > >                   & JJ,
		vector < SCIP_Real >                        & zval,
		vector < SCIP_Real >                          x_lb,
		vector < SCIP_Real >                          x_ub,
		SCIP_Real                                    z_lb,
		SCIP_Real                                    z_ub,
		int                                          NNN,
		int                                          j,
		int                                          t,
		SCIP_Real                                    y_ub
		)
{
	int lenght = pow(_K,_T-1)/pow(_K,t-1);
	SCIP_Real zcoeff = 0;
	//SCIP_Real zval = 0;
	vector < vector < SCIP_Real > > ycoeff_max (_I);
	vector < vector < int > > kk_max (_I);
	vector < vector < int > > jj_max (_I);
    vector < SCIP_Real > b (_I);
    vector < SCIP_Real > a (_I);
	for (int q=0; q< pow(_K,_T-t); q++)
	{
		zcoeff += gamma[j*lenght+q][t-2];
	}
	assert (z_lb <= z_ub);
	for (int i =0 ; i<_I ; i++)
	{
		cout << x_lb[i] << endl;
		cout << x_ub[i] << endl;
	}
	//Find zval
	//if (z_lb == z_ub )
		//zval = z_lb;

	//Calculate ycoeff
	for (int i=0 ; i<_I ; i++)
	{
		ycoeff[i].resize(NNN+1);
		//ycoeff_max[i].resize(pow(_K,NNN), -100000);
		//kk_max[i].resize(pow(_K,NNN));
		//jj_max[i].resize(pow(_K,NNN));
		for (int kk=0; kk <NNN+1; kk++)
		{
			ycoeff[i][kk].resize(pow(_K,kk));
			int len = pow(_K,t+kk-1)/pow(_K,t-1);
			for(int jj=0; jj < pow(_K,kk) ; jj++)
			{
				for (int q=0; q<pow(_K,_T-t-kk) ; q++)
				{
					ycoeff[i][kk][jj] += pi[t+kk-2][j*len + jj][q];
				}
				cout << "ycoe" << ycoeff[i][kk][jj] << endl;
			}
		}
	}

	//Calculate ycoeff_max, and set kk,jj_max
	vector < double > y_max_1(_I);
	for (int i=0; i<_I ; i++)
	{
		for (int q=0 ; q<pow(_K,NNN) ; q++)
		{
			y_max_1[i] += ycoeff[i][NNN][q];
		}
		cout << "ymax1" << y_max_1[i] << endl;
		cout << "ycoeff" << ycoeff[i][0][0] << endl;
		if (ycoeff[i][0][0] >= y_max_1[i])
		{
			ycoeff_max[i].push_back(ycoeff[i][0][0]) ;
			//cout << "max" << ycoeff_max[i][q] << endl;
			kk_max[i].push_back(0) ;
			jj_max[i].push_back(0) ;
		}
		else
		{
			for (int qq=0 ; qq<pow(_K,NNN) ; qq++)
			{
				ycoeff_max[i].push_back(ycoeff[i][1][qq]);
				kk_max[i].push_back(1) ;
				jj_max[i].push_back(qq) ;
			}
		}
	}
	/*cout << "ymax" << ycoeff_max[i][q] << endl;
	for (int kk=0 ;kk<NNN+1 ; kk++)
	{
		if (t==3&&j==1)
			cout<<"test";
		cout << "ycoeff" << ycoeff[i][kk][kk*q] << endl;
		double temp_10 = ycoeff[i][kk][kk*q];
		double temp_20 = ycoeff_max[i][q];
		double temp_1 = floor(temp_10 * 1000000000.0) / 1000000000.0;
		double temp_2 = floor(temp_20 * 1000000000.0) / 1000000000.0;
		if ( temp_1 >= temp_2)
		{
			ycoeff_max[i][q] = ycoeff[i][kk][kk*q];
			//cout << "max" << ycoeff_max[i][q] << endl;
			kk_max[i][q] = kk ;
			jj_max[i][q] = kk*q ;
			//cout << jj_max[i][q] << endl;
			//cout << kk_max[i][q] << endl;
		}
	}
}
}//////////*/ //up to here!

	cout << "max" << ycoeff_max[0][0] << "kk" << kk_max[0][0] << "jj" << jj_max[0][0] << endl;
	cout << "maz" << ycoeff_max[0][1] << "kk" << kk_max[0][1] << "jj" << jj_max[0][1] << endl;

	//Check if any of ycoeff_max is zero
	vector < int > nzero(_I);
	vector < vector < int > > KK_zero(_I);
	vector < vector < int > > JJ_zero(_I);
	for (int ii=0; ii < _I ; ii++)
	{
		for (int jj=0 ; jj< (int) kk_max[ii].size() ; jj++)
		{
			if (ycoeff[ii][kk_max[ii][jj]][jj_max[ii][jj]] < 0.000000000000001)
			{
				KK_zero[ii].push_back(kk_max[ii][jj]);
				JJ_zero[ii].push_back(jj_max[ii][jj]);
				nzero[ii] ++;
			}
		}
	}//here

	for (int i=0; i<_I ; i++)
	{
		if (x_lb[i] == x_ub[i]) //See if x is fixed
		{
			if (z_lb == z_ub && z_lb == 0) //See if z is fixed
			{
				SCIP_CALL(Useall(KK[i], JJ[i], kk_max[i], jj_max[i]));
				zval[i] = 0;
			}
			else if (z_lb == z_ub && z_lb == 1)
			{
				if (nzero[i] > 0)
					SCIP_CALL(Usesome(i , KK[i], JJ[i] , kk_max[i], jj_max[i], ycoeff));
				else
					SCIP_CALL(Useall(KK[i], JJ[i], kk_max[i], jj_max[i]));
				zval[i] = 1;
			}
			else // z is not fixed
			{
				zval[i] = Choosetouse(nzero[i] , i , zcoeff, KK[i], JJ[i], kk_max[i], jj_max[i], ycoeff);
			}
		}
		else //x is not fixed
		{
			if (z_lb == z_ub && z_lb == 0)
			{
				SCIP_CALL(Useall(KK[i], JJ[i], kk_max[i], jj_max[i]));
				zval[i] = 0;
			}
			else if (z_lb == z_ub && z_lb == 1)
			{
				if (nzero[i] > 0)
					SCIP_CALL(Usesome(i , KK[i], JJ[i] , kk_max[i], jj_max[i], ycoeff));
				else
					SCIP_CALL(Useall(KK[i], JJ[i], kk_max[i], jj_max[i]));
				zval[i] = 1;
			}
			else
			{
				zval[i] = Choosetouse(nzero[i] , i , zcoeff, KK[i], JJ[i], kk_max[i], jj_max[i], ycoeff);
			}
		}
	}



	/*for (int i =0; i<_I ; i++)
	{
		for (int q=0 ; q< pow(_K,NNN) ; q++)
		{
			if (q==1)
				if (kk_max[i][0] == kk_max[i][1] && jj_max[i][0] == jj_max[i][1])
					break;
			//cout << "ymax" << ycoeff_max[i] << endl;
			if (ycoeff_max[i][q] >= 0)
			{
				a[i] += ycoeff_max[i][q];
				kk_a[i].push_back(kk_max[i][q]);
				jj_a[i].push_back(jj_max[i][q]);
			}
			else
			{
				b[i] += ycoeff_max[i][q];
				kk_b[i].push_back(kk_max[i][q]);
				jj_b[i].push_back(jj_max[i][q]);
			}
		}
		if (a[i] == 0)
		{
			cout << "SET EVERYTHING TO ZERO" << endl;
			allzero[i] = true;
			//return zval;
		}
		else if (b[i] < 0)
		{
			//See if z is set to zero, if so at all scenarios have to be opened
			if (z_lb == z_ub && z_lb == 0)
			{
				for (int q=0 ; q< (int) kk_a[i].size() ; q++)
				{
					KK[i].push_back(kk_a[i][q]);
					JJ[i].push_back(jj_a[i][q]);
				}
				for (int q=0 ; q< (int) kk_b[i].size() ; q++)
				{
					KK_beta[i].push_back(kk_b[i][q]);
					JJ_beta[i].push_back(jj_b[i][q]);
				}
				//zval = 0;
			}
			//see if  z is set 1
			else if (z_lb == z_ub && z_ub == 1)
			{
				for (int q = 0; q < (int) kk_a.size(); q++) {
					KK[i].push_back(kk_a[i][q]);
					JJ[i].push_back(jj_a[i][q]);
				}
				//zval =1;
			}
			//if z is not fixed yet
			else if (z_lb != z_ub)
			{
				if (zcoeff == 0)
				{
					for (int q=0 ; q< (int) kk_a[i].size() ; q++)
					{
						KK[i].push_back(kk_a[i][q]);
						JJ[i].push_back(jj_a[i][q]);
					}
					//zval = 1;
				}
				else if (zcoeff > 0)
				{
					if (b[i]*(y_ub-_Wastbeta) < zcoeff)
					{
						for (int q=0 ; q< (int) kk_a[i].size() ; q++)
						{
							KK[i].push_back(kk_a[i][q]);
							JJ[i].push_back(jj_a[i][q]);
						}
						for (int q=0 ; q< (int) kk_b[i].size() ; q++)
						{
							KK_beta[i].push_back(kk_b[i][q]);
							JJ_beta[i].push_back(jj_b[i][q]);
						}
					}
					else
					{
						for (int q = 0; q < (int) kk_a.size(); q++)
						{
							KK[i].push_back(kk_a[i][q]);
							JJ[i].push_back(jj_a[i][q]);
						}
						//zval = 1;
					}
				}
				else
					cout << "There is a negative zcoeff!!"<< endl;
			}
			else
			{
				cout << "There is sth wrong!" << endl;
			}
		}
	}*/


	return SCIP_OKAY;
}

SCIP_RETCODE ObjPricervaccinebp::Useall (
		vector < int >                & KK1,
		vector < int >                & JJ1,
		vector < int >                  kk_max1,
		vector < int >                  jj_max1)
//vector < SCIP_Real >          & zval)
{

	for (int jj=0 ; jj< (int) kk_max1.size() ; jj++)
	{
		KK1.push_back(kk_max1[jj]);
		JJ1.push_back(jj_max1[jj]);
	}

	return SCIP_OKAY;
}
SCIP_RETCODE ObjPricervaccinebp::Usesome(
		int                                          i,
		vector < int >                             & KK1,
		vector < int >                             & JJ1,
		vector < int >                               kk_max1,
		vector < int >                               jj_max1,
		vector < vector < vector < SCIP_Real > > >   ycoeff	)
{
	for (int jj=0 ; jj< (int) kk_max1.size() ; jj++)
	{
		if (ycoeff[i][kk_max1[jj]][jj_max1[jj]] >= 0.000000000000001)
		{
			KK1.push_back(kk_max1[jj]);
			JJ1.push_back(jj_max1[jj]);
		}
	}
	return SCIP_OKAY;
}

SCIP_Real ObjPricervaccinebp::Selectrandom (
		int                                            i,
		vector < int >                               & KK1,
		vector < int >                               & JJ1,
		vector < int >                               kk_max1,
		vector < int >                               jj_max1,
		vector < vector < vector < SCIP_Real > > >   ycoeff	)
{
	int randval =  rand() % 100 +1 ;
	//If it is less than 50, set z=0
	if (randval <= 50)
	{

			for (int q=0 ; q< (int)kk_max1.size()  ; q++)
			{
				KK1.push_back(kk_max1[q]);
				JJ1.push_back(jj_max1[q]);
			}

		return 0;
	}
	else
	{

			for (int jj=0 ; jj< (int) kk_max1.size() ; jj++)
			{
				if (ycoeff[i][kk_max1[jj]][jj_max1[jj]] >= 0.000000000000001)
				{
					KK1.push_back(kk_max1[jj]);
					JJ1.push_back(jj_max1[jj]);
				}
			}

		return 1;
	}
}

SCIP_Real ObjPricervaccinebp::Choosetouse(
		int                                          nzero,
		int                                          i,
		SCIP_Real                                    zcoeff,
		vector < int >                             & KK1,
		vector < int >                             & JJ1,
		vector < int >                               kk_max1,
		vector < int >                               jj_max1,
		vector < vector < vector < SCIP_Real > > >   ycoeff
)
{
	if ( nzero > 0)
	{
		if ((int) kk_max1.size() == nzero)
		{
			return 0;
		}
		else
		{
			if (zcoeff == 0)
			{
			   return Selectrandom ( i, KK1, JJ1, kk_max1 , jj_max1, ycoeff);
			}
		}
	}
	else
	{
		for (int q=0 ; q < (int)kk_max1.size() ; q++)
		{
			KK1.push_back(kk_max1[q]);
			JJ1.push_back(jj_max1[q]);
		}
		//cout << "kk" << KK1.size() <<endl;
		//cout << "kk" << (int) KK1.size() <<endl;
		//cout << "kksize" << KK[i].size() << endl;
		return 0;
	}
}

