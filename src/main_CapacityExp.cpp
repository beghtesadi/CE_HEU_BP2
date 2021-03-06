/**@file   main_vaccinebp.cpp
 * @brief  Main file for vaccinebp pricing example
 * @author Bahareh Eghtesadi
 *
 *  This file contains the \ref main() main function of the projects. This includes all inputs of the problem and all the default plugins of
 *  \SCIP and the once which belong to that projects. After that is starts the interactive shell of \SCIP or processes
 *  the shell arguments if given.
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "scip/scipshell.h"

/* user defined includes */
#include "pricer_CapacityExp.h"
#include "probdata_CapacityExp.h"
#include "vardata_CapacityExp.h"
#include "GeDemand_CapacityExp.h"
#include "branch_CapacityExp.h"
#include "conshdlr_CapacityExp.h"
#include "conshdlr2_CapacityExp.h"



/* namespace usage */
using namespace std;
using namespace scip;


SCIP_RETCODE Printresults( SCIP* scip, int T , int I , int NScen ,
		vector< SCIP_VAR* >        z_var,
		vector< double>             Nodeprobe,
		vector< vector< int > >     nodemat,
		vector < double >           purchcost
		);
double mean=5;
double K=2;


int main ()
{

  vector < double > purchcost = { 1,4,7}; //{4000}; //{1000,4000,7000};
  vector < int > vialsize = {1, 5,10}; //{5}; //{1,5,10};
  int I=3;
  double ConfRate=0.75;
  double Shortbeta=1;
  double Wastbeta=1;
  double M=100;
  int Nper_to_spoil=1;
  int T=4; //T is the number of stages
  int NScen=pow(K, T-1);
  int Nnodes=0;
  for (int i=0; i<T ; i++){
	Nnodes=Nnodes+pow(K,i);}

  vector<double>Dprob;
  vector<int>Demand;
  vector<double>Nodedemand;
  vector<int>ScenDemand;
  vector<double>Scenprob;
  vector<double>Scenprob2;
  vector<double>Nodeprobe;

  GenerateDemand(Demand, Dprob, mean, K);

  //generating the nodes matrix
  vector< vector<int> > nodemat;
  vector< vector<int> > jmat;

	int n=1;
	for (int t=1 ; t<T+1 ; t++){
		vector<int>row; //creat an empty vector
                vector<int>jrow;
		int len=NScen/pow(K,t-1);
		for (int i=0 ; i<pow(K,t-1) ; i++){
			for (int j=0 ; j<len ; j++){
				row.push_back(n);
                                jrow.push_back(i);}
		n++;}
		nodemat.push_back(row);
                jmat.push_back(jrow);}


	for (int i=0 ; i<pow(K,T-1) ; i++){
		Scenprob2.push_back(0);
	    Scenprob.push_back(0);}


	for (int j=0 ; j<K ; j++){
		Nodedemand.push_back(Demand[j]);
		Nodeprobe.push_back(Dprob[j]);
		Scenprob[j]=Dprob[j];}

	for (int t=3; t<T+1 ; t++){
		int nn=0;
		for (int i=0 ; i<pow(K,t-2) ; i++){
			for (int j=0 ; j < K ; j++){
				Nodedemand.push_back(Demand[j]);
				Nodeprobe.push_back(Dprob[j]);
				Scenprob2[nn]=Scenprob[i]* Dprob[j];
				nn++;}}
		int ssize= (int)Scenprob2.size();
		for(int ii=0 ; ii<ssize ; ii++){
			Scenprob[ii]=Scenprob2[ii];}
	}

	//Generating ScenDemand
	for (int i=0 ; i<pow(K,T-2) ; i++){
		for (int j=0; j<K ; j++){
			ScenDemand.push_back(Demand[j]);}}

  SCIP* scip = NULL;
  const char* probname = "vaccinebp";


  SCIP_CALL( SCIPcreate(&scip) );

  // include default plugins 
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   
   // set verbosity parameter 
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", FALSE) );
   // for column generation instances, disable restarts 
   SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );
   SCIP_CALL(SCIPsetHeuristics(scip, SCIP_PARAMSETTING_DEFAULT , TRUE));

   
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   assert(scip != NULL);

   // create problem in SCIP and add non-NULL callbacks via setter functions 
   SCIP_CALL( SCIPcreateProb(scip, probname, 0, 0, 0, 0, 0, 0, 0) );

   // set objective sense 
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );



   SCIP_ProbData* probdata = new SCIP_ProbData;
   cout<< probdata<<endl;
   vector< SCIP_VARDATA* > zvardata;
   //SCIP_VAR* wvar =NULL
   int n_vars=0;
   char var_name[255];
   char con_name[255];
   char name[255];
 
   
   // tell SCIP that the objective will be always integral 
   //SCIP_CALL( SCIPsetObjIntegral(scip) );

   //SCIP_CALL( SCIPallocBufferArray(scip, &conss, nitems) );

  vector< SCIP_VAR* > z_var(NScen);
  //vector< <vector < vector< SCIP_Real > > > zycol (I, vector < vector< SCIP_Real> > (1, vector< SCIP_Real> (1)) );
  vector < SCIP_Real > z_xcol (I);
  
  SCIP_Real zzcol = 0;
  SCIP_Bool zvartype = true;
  SCIP_VARDATA* zvardata_ptr;
  for (int i=0; i<NScen; i++)
  {
     zvardata_ptr = new SCIP_VARDATA;
     zvardata.push_back(zvardata_ptr);
     SCIP_VAR* var;
     SCIPsnprintf(var_name, 255, "Z%d", i );
     SCIP_CALL( SCIPcreateVar(scip,
                     &var,                      // returns new index
                     var_name,                  // name
                     0.0,                       // lower bound
                     1.0,                       // upper bound
                     0,                         // objective
                     SCIP_VARTYPE_BINARY,       // variable type
                     true,                     // initial
                     false,                    // forget the rest ...
                     0, 0, 0, 0, 0) );
    SCIP_CALL( SCIPaddVar(scip, var) );
    z_var[i] = var;
    //Add the var to the probdata and vardata
    SCIP_CALL( SCIPvardataCreatevaccinebp (scip, zvardata[i], z_xcol, zzcol, zvartype, 0, I, K, T, i) );
    SCIPvarSetData(var, zvardata[i]);
   }
  
  //ADD constraint 5b

  vector< vector < vector <SCIP_CONS* > > > sh_con(T-1);
 
  //SCIP_Real mwj =0 ;
  for(int t=2; t<T+1; t++)
  {
    sh_con[t-2].resize(pow(K,t-1));
    int lenght=NScen/pow(K,t-1); 
    for (int j=0 ; j<pow(K,t-1) ; j++)
      {
         sh_con[t-2][j].resize(pow(K,T-t));
         for (int q=0 ; q<pow(K,T-t) ; q++)
            {
                SCIP_CONS* con;
                SCIPsnprintf(con_name, 255, "sh%d_%d_%d", t, j, q);
                SCIP_VAR* index = z_var[j * pow(K,T-t) + q];
                SCIP_Real coeff = Nodedemand[nodemat[t-1][j*lenght]-2];
                SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 1, &index, &coeff,
                     Nodedemand[nodemat[t-1][j*lenght]-2]-Shortbeta,   
                     SCIPinfinity(scip),     
                     true,                   
                     false,                  
                     true,                   
                     true,                  
                     true,                   
                     false,                  
                     true,                   
                     false,                 
                     false,                 
                     false) );               
                SCIP_CALL( SCIPaddCons(scip, con) );
                sh_con[t-2][j][q] = con;
             }
       }
   }
   //Add constraint 5c
   vector<vector<SCIP_CONS* > >z_con (NScen); 
 
   //char con_name[255];  
   for (int i=0; i<NScen; i++)
   {  
      z_con[i].resize(T-1);
      for (int j=2; j<T+1 ; j++)
       {
          SCIP_CONS* con;
          SCIPsnprintf(con_name, 255, "z%d", i);
          SCIP_VAR* index = z_var[i];
          SCIP_Real coeff = 1;
          SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 1, &index, &coeff,
                     0,                     
                     SCIPinfinity(scip),    
                     true,                   
                     false,                  
                     true,                   
                     true,                  
                     true,                   
                     false,                  
                     true,                   
                     false,                  
                     false,                  
                     false) );               
                SCIP_CALL( SCIPaddCons(scip, con) );
                z_con[i][j-2] = con;
        }
    }

   //Add constraint 5d 
   SCIP_CONS*  alpha_con = NULL;
   SCIP_CALL( SCIPcreateConsLinear(scip, &alpha_con,"AlphaCon" , 0, NULL, NULL,
                     -SCIPinfinity(scip),    
                     1-ConfRate,                    
                     true,                  
                     false,                  
                     true,                  
                     true,                   
                     true,                   
                     false,                  
                     true,                   
                     false,                  
                     false,                 
                     false) );              

    for (int i = 0; i < NScen; ++i)
	 SCIP_CALL( SCIPaddCoefLinear(scip, alpha_con, z_var[i], Scenprob[i]) );
    SCIP_CALL( SCIPaddCons(scip, alpha_con) );

    //Add constraint 5e 
   //SCIP_Real wje =1;
   vector< vector<SCIP_CONS* > > w_con(T-1);
   for (int t=2; t< T+1 ; t++)
   {
    w_con[t-2].resize(pow (K, t-1));
    for (int j=0; j < pow (K, t-1); j++)
     {
      SCIP_CONS* con = NULL;
      SCIPsnprintf(con_name, 255, "w%d_%d", t, j);
      SCIP_CALL( SCIPcreateConsLinear( scip, &con, con_name, 0, NULL, NULL,
    		                           1.0,
                                       1.0,   
                                       true,  
                                       false, 
                                       true,  
                                       true,  
                                       true,  
                                       false, 
                                       true,  
                                       false, 
                                       false, 
                                       false) );
      SCIP_CALL( SCIPaddCons(scip, con) );
      w_con[t-2][j]=con;
     }
   }
   
  vector< SCIP_VAR*> null_var(n_vars, NULL);
   //Createproblemdata
   SCIP_CALL( SCIPprobdataCreate (scip, probdata, z_var, null_var,  n_vars, NScen ,T, K, I, jmat) );

    // create start solution 
   vector< SCIP_VAR*> wj_var;
   //Creating vardata
   vector< vector< SCIP_VARDATA*> > w_vardata;
   vector< vector< SCIP_VARDATA*> > w_vardatazero;
   SCIP_VARDATA* wvardata_ptr;
   SCIP_VARDATA* wvardatazero_ptr;
   for (int t=2; t<T+1; t++)
   { 
     vector < SCIP_VARDATA* > row_wvardata;
     vector < SCIP_VARDATA* > row_wvardatazero;
     for(int j=0; j<pow(K,t-1) ; j++)
     {
      wvardata_ptr = new SCIP_VARDATA;
      wvardatazero_ptr = new SCIP_VARDATA;
      row_wvardata.push_back(wvardata_ptr);
      row_wvardatazero.push_back(wvardatazero_ptr);
     }
     w_vardata.push_back(row_wvardata);
     w_vardatazero.push_back(row_wvardatazero);
    }
   SCIP_Bool wvartype = false;


   for (int t=2; t<T+1; t++)
   { 
     int lenght=NScen/pow(K,t-1); 
     int NNN = CalMin (T-t,Nper_to_spoil);
     for(int j=0; j<pow(K,t-1) ; j++)
     {
       SCIP_VAR * var;
       vector< SCIP_Real> xcolumnz (I);
       vector< vector< vector< SCIP_Real > > > ycolumnz (I);  
       SCIP_Real zerocolumnz = 0; 
     
       for (int i=0; i<I; i++){
         xcolumnz[i] = 0;
         ycolumnz[i].resize(NNN+1);
         for (int kk=0; kk<NNN+1; kk++){
          ycolumnz[i][kk].resize(pow(K,kk)); 
          for (int q = 0 ; q< pow(K,kk) ; q++)
           {
            ycolumnz[i][kk][q]= 0;
            }
           }
          }
	   SCIP_Real zcolumnz;
	   if (t==T && j==0)
           zcolumnz = 1 ;
	   else
            {
              zcolumnz = 0 ;
            }
	   if (t==T && j==0)
		   ycolumnz[1][0][0] = 0;  //ycolumnz[0][0][0] = 0; //ycolumnz[1][0][0] = 0;
	   else
		   ycolumnz[1][0][0] = Nodedemand[nodemat[t-1][j*lenght]-2];  //ycolumnz[0][0][0] = Nodedemand[nodemat[t-1][j*lenght]-2]; //ycolumnz[1][0][0] = Nodedemand[nodemat[t-1][j*lenght]-2];
	   xcolumnz[1] = 1;    //xcolumnz[0] = 1; // xcolumnz[1] = 1;
       SCIPsnprintf(name, 255, "ww_%d", nodemat[t-1][j*lenght]);
       SCIPdebugMessage("create variable for node %d \n", nodemat[t-1][j*lenght]);
       SCIP_Real objj = Nodeprobe[nodemat[t-1][j*lenght]-2]*purchcost[1]; ////SCIP_Real objj = Nodeprobe[nodemat[t-1][j*lenght]-2]*purchcost[1];

       SCIP_CALL( SCIPcreateVarvaccinebp(scip, &var, name, objj , TRUE, TRUE, w_vardatazero[t-2][j]) );  ///?????????????var pointer and vardata


      // create the variable data for the variable
      SCIP_CALL( SCIPvardataCreatevaccinebp(scip, w_vardatazero[t-2][j], xcolumnz, zcolumnz, wvartype, NNN, I, K, t, j ));
       
      // add variable to the problem 
      SCIP_CALL( SCIPaddVar(scip, var) );
      if (var == NULL)
    	  cout<< "Var is NULL" << endl;

      // store variable in the problme data 
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, var, NULL) );

      // add variable to corresponding constraints 
      SCIP_Real wcoef = Nodedemand[nodemat[t-1][j*lenght]-2];
	  SCIP_Real wcoef2=0;
      SCIP_Real wje = 1;
      if (t==T && j==0)
		  SCIP_CALL( SCIPaddCoefLinear(scip, sh_con[t-2][j][0], var, wcoef2) );
	  else
	  {
         for (int q=0 ; q<pow(K,T-t) ; q++)
         {
           SCIP_CALL( SCIPaddCoefLinear(scip, sh_con[t-2][j][q], var, wcoef) );
         }
	  }
      SCIP_CALL( SCIPaddCoefLinear(scip, w_con[t-2][j], var, wje ) );
      SCIP_Real mwj =0 ;
      for (int jj=0 ; jj< pow(K,T-t); jj++)
      {
		  if (j*pow(K,T-t)+jj==0 && t==T)
			  SCIP_CALL( SCIPaddCoefLinear(scip, z_con[j*pow(K,T-t)+jj][t-2], var, -1 ) );
		  else
              SCIP_CALL( SCIPaddCoefLinear(scip, z_con[j*pow(K,T-t)+jj][t-2], var, 0 ) );
      }
     

      // add the variable data to the variable 
      //SCIPvarSetData(var, w_vardata[t-2][j]);
      wj_var.push_back(var);
      n_vars++;
    }
   }  //// For t, j
  
   // set user problem data 
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   //include conshdlr
    const char* vaccinebp_conshdlr_NAME = "vaccinebp_Conshdlr";

    Objconshdlrvaccinebp* vaccinebp_conshdlr_ptr = new Objconshdlrvaccinebp(scip,vaccinebp_conshdlr_NAME );

    assert(vaccinebp_conshdlr_ptr != NULL);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, vaccinebp_conshdlr_ptr, true ));

    //include second conshdlr
    const char* vaccinebp_conshdlrsec_NAME = "vaccinebp_Conshdlrsec";

    Objconshdlrvaccinebpsec* vaccinebp_conshdlrsec_ptr = new Objconshdlrvaccinebpsec(scip,vaccinebp_conshdlrsec_NAME );

    assert(vaccinebp_conshdlrsec_ptr != NULL);

    SCIP_CALL(SCIPincludeObjConshdlr(scip, vaccinebp_conshdlrsec_ptr, true ));


   //include branch rule
   const char* vaccinebp_branchrule_name = "vaccinebp_Branchrule";

   ObjBranchrulevaccinebp* vaccinebp_branchrule_ptr = new ObjBranchrulevaccinebp(scip, vaccinebp_branchrule_name);
   //cout << "vaccinebp_branchrule" << vaccinebp_branchrule_ptr << endl;

   assert(vaccinebp_branchrule_ptr != NULL);
   SCIP_CALL(SCIPincludeObjBranchrule(scip, vaccinebp_branchrule_ptr, true));


   // include pricer
   static const char* vaccinebp_PRICER_NAME = "vaccinebp_Pricer";
   int p_n_pricing = 0;
      
   ObjPricervaccinebp* vaccinebp_pricer_ptr = new ObjPricervaccinebp (scip, vaccinebp_PRICER_NAME, Nper_to_spoil, Wastbeta,T, Nnodes, K, I, Nodeprobe, purchcost, vialsize, jmat, nodemat,
		   z_var, wj_var, sh_con, z_con, alpha_con, w_con, n_vars, 0);

    SCIP_CALL( SCIPincludeObjPricer(scip, vaccinebp_pricer_ptr, true) );

   // activate pricer 
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip,vaccinebp_PRICER_NAME)) );


   /***********************
    * Version information *
    ***********************/


 // Up to here  

   SCIP_CALL( SCIPwriteOrigProblem(scip, "CapacityExp_init.lp", "lp", FALSE) );

   /*************
    *  Solve    *
    *************/


   SCIP_CALL( SCIPsolve(scip) );


   /**************
    * Statistics *
    *************/
   
   FILE* fp = fopen("/home/bahar/SCIPproject/vaccinebp-cons-heu/results","w");
  SCIP_CALL( SCIPprintStatistics(scip, fp) );

   SCIP_CALL( SCIPprintBestSol(scip, fp, FALSE) );

   SCIP_CALL (Printresults(scip,  T ,  I , NScen , z_var, Nodeprobe , nodemat , purchcost));

   /********************
    * Deinitialization *
    ********************/

   /* free local buffer arrays 
   SCIPfreeBufferArray(scip, &w_con);
   SCIPfreeBufferArray(scip, &sh_con);
   SCIPfreeBufferArray(scip, &z_con);
   SCIPfreeBufferArray(scip, &alpha_con);*/
   
   delete probdata;
   delete wvardata_ptr;
   delete zvardata_ptr;
   /*for (int i=0; i<NScen ; i++)
   {
    delete zvardata[i];
    }*/
   

  SCIP_CALL( SCIPfree(&scip) );


   BMScheckEmptyMemory();

   return 0;

}//for main 
SCIP_RETCODE Printresults( SCIP* scip, int T , int I , int NScen ,
		vector< SCIP_VAR* >        z_var,
		vector< double>             Nodeprobe,
		vector< vector< int > >     nodemat,
		vector < double >           purchcost
		)
{
	SCIP_SOL* sol;
	sol = SCIPgetBestSol (scip);

	SCIP_PROBDATA* probdata;
	probdata = SCIPgetProbData(scip);
	assert(probdata != NULL);
	vector< SCIP_VAR* > wcand;
	SCIP_VARDATA* vardata;
	vector < vector < vector < SCIP_Real > > > wcandv(T-1 );
	vector < vector < int > > nwcand(T-1);
	vector < vector < vector < SCIP_VARDATA* > > > xvardata(T-1);
	for (int t=2 ; t < T+1 ; t++)
	{
		wcandv[t-2].resize(pow(K,t-1));
		nwcand[t-2].resize(pow(K,t-1));
		xvardata[t-2].resize(pow(K,t-1));
	}

	SCIP_CALL( SCIPprobdataGetVars(probdata,wcand));
	for (int i=0 ; i< SCIPprobdataGetNVars(probdata) ; i++)
	{
		vardata = SCIPvarGetData(wcand[i]);
		int stage = SCIPvardatagetstage(vardata);
		int state = SCIPvardatagetstate(vardata);
		wcandv[stage-2][state].push_back(SCIPgetSolVal(scip, sol, wcand[i]));
		if (SCIPgetSolVal(scip, sol, wcand[i]) > 0)
			cout << "alpha" << SCIPgetSolVal(scip, sol, wcand[i]) <<  "   " << SCIPvarGetLPSol(wcand[i]) << endl;
		nwcand[stage-2][state] ++ ;
	}
	//vector < vector < vector < SCIP_Real > > > zcolumn_w(SCIPprobdataGetT(probdata)-1);
	vector < vector < vector < vector < SCIP_Real > > > > xcolumn_w(T-1);

	for (int t=2 ; t < T+1 ; t++)
	{
		xcolumn_w[t-2].resize(pow(K,t-1));
		//zcolumn_w[t-2].resize(pow(K,t-1));
		for (int j=0 ; j < pow(K,t-1) ; j++)
		{
			xcolumn_w[t-2][j].resize(nwcand[t-2][j]);
			//zcolumn_w[t-2][j].resize(nwcand[t-2][j]);
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
					SCIP_CALL(SCIPvardataGetxcolumn(xvardata[t-2][j][v], xcolumn_w[t-2][j][v]));
					//cout << wcandv[t-2][j][v] << endl;
					//cout << "xx1-" << xcolumn_w[t-2][j][v][0] << endl;
					//cout << "xx5-" << xcolumn_w[t-2][j][v][1] << endl;
					//cout << "xx10-" << xcolumn_w[t-2][j][v][2] << endl;
					//cout << "wval-" << wcandv[t-2][j][v] << endl;
				}
			}
		}
	SCIP_Real objval =0;
		for (int i=0 ; i < I ; i++)
		{
			for (int t=2 ; t < T+1 ; t++)
			{
				int lenght = pow(K,T-1)/pow(K, t-1);
				for (int j=0 ; j < pow(K,t-1) ; j++)
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
							objval += xcolumn_w[t-2][j][n][i] * wcandv[t-2][j][n]* purchcost[i] * Nodeprobe[nodemat[t-1][j*lenght]-2];
						}
						if (xcandfrac > 0)
						{
							cout << "x_" <<  i << "_" << t << "_" << j << " " << xcandfrac << endl;
						}
					}
				}
			}
		}
		cout << "objval " << objval << endl;

		for (int i=0 ; i < NScen ; i++)
		{
			if (SCIPgetSolVal(scip, sol, z_var[i]) > 0 )
				cout << "z_" << i << SCIPgetSolVal(scip, sol, z_var[i]) << endl;
		}
		return SCIP_OKAY;
}

