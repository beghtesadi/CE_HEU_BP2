/**@file pricer_vaccinebp.h
 * @brief pricer plugin
 * @author Bahareh Eghtesadi
 */
#ifndef __SCIP_PRICER_vaccinebp_H__
#define __SCIP_PRICER_vaccinebp_H__

#include "objscip/objscip.h"
#include "scip/pub_var.h"

#include <vector>
#include <list>

using namespace std;
using namespace scip;


/** pricer class */
class ObjPricervaccinebp : public ObjPricer
{
public:

   /** Constructs the pricer object with the data needed */
   ObjPricervaccinebp(
   SCIP*                                                                scip, 
   const char*                                                         p_name,         /**< SCIP pointer */
   //vector < vector < SCIP* > >                                        & subscip,
   //vector < vector < SCIP_VAR* > >                                    & p_zvar,
   //vector < vector < SCIP_Real > >                                    & p_zvar_lb,
   //vector < vector < SCIP_Real > >                                    & p_zvar_ub,
   //vector < vector < vector < SCIP_VAR* > > >                         & p_xvar, 
   //vector < vector < vector < SCIP_Real > > >                         & p_xvar_lb,
   //vector < vector < vector < SCIP_Real > > >                         & p_xvar_ub,
   //vector < vector < vector < vector < vector< SCIP_VAR* > > > > >   & p_yvar, 
   //vector < vector < vector < vector < vector< SCIP_Real > > > > >   & p_yvar_lb,
   //vector < vector < vector < vector < vector< SCIP_Real > > > > >   & p_yvar_ub,
   //vector < vector < vector < SCIP_CONS* > > >                        & p_wcon,             // Wastage constraint
   //vector < vector < vector < vector < SCIP_CONS* > > > >             & p_secon,     //Second constrain       
   int                                                                  p_Nper_to_spoil,
   double                                                               p_Wastbeta,
   int                                                                  p_T,
   int                                                                  p_Nnodes,
   double                                                               p_K,
   int                                                                  p_I,
   vector<double>                                                     & p_Nodeprobe,
   vector< double >                                                   & p_purchcost,
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
   );

   /** Destructs the pricer object. */
   virtual ~ObjPricervaccinebp();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

   /** farkas pricing method of variable pricer for infeasible LPs */
   //virtual SCIP_DECL_PRICERFARKAS(scip_farkas);

   /** perform pricing */
   SCIP_RETCODE pricing(SCIP* scip);            /**< whether we perform Farkas pricing */


   /** add tour variable to problem */
   SCIP_RETCODE add_newcolumn_variable(
   SCIP*                                                           scip,                                  
   vector< SCIP_Real >                                           & x_newcolumn,
   vector< vector< vector < SCIP_Real > > >                      & y_newcolumn, 
   SCIP_Real                                                      z_newcolumn,
   SCIP_VARDATA*                                                  vardata,
   int                                                             t,
   int                                                             j 
    );

   /** return negative reduced cost tour (uses restricted shortest path dynamic programming algorithm) */
   SCIP_RETCODE find_new_column(
   SCIP*                                                           scip, 
   vector< vector < vector <SCIP_Real > > >         & pi,      
   vector< vector < SCIP_Real > >                   & gamma,
   vector< vector <vector< SCIP_Real > > >         & x_newcolumn,
   vector< vector <vector< vector< vector< SCIP_Real > > > > > & y_newcolumn, 
   vector< vector < SCIP_Real > >                 & z_newcolumn,
   vector< vector < SCIP_Real > >                  & objval );

   /*SCIP_RETCODE SCIPchangepricerlb(
         SCIP*                scip, 
         SCIP_Bool            candtype,
         int                  vialtype,
         int                  stage,
         int                  state,
         SCIP_Real            newbound
   );

   SCIP_RETCODE SCIPchangepricerub(
         SCIP*                scip, 
         SCIP_Bool            candtype,
         int                  vialtype,
         int                  stage,
         int                  state,
         SCIP_Real            newbound
   );*/
   SCIP_RETCODE addBranchingDecisionConss(
      SCIP*                                                                  scip,
      SCIP_CONSHDLR*                                                         conshdlr,
      vector < vector < vector < SCIP_Real > > >                             & _xvar_lb,
      vector < vector < vector < SCIP_Real > > >                             & _xvar_ub,
      vector < vector < SCIP_Real > >                                        & _zvar_lb,
      vector < vector < SCIP_Real > >                                        & _zvar_ub
      );

   SCIP_RETCODE findyandz(
   		vector< vector < vector < SCIP_Real > > >    pi,      /**< dual variable value*/
   		vector< vector < SCIP_Real > >               gamma,
   		vector < vector < vector < SCIP_Real > > >     &ycoeff,
   		vector < vector < int > >                   & KK,
   		vector < vector < int > >                   & JJ,
   		//vector < vector < int > >                   & KK_beta,
   		//vector < vector < int > >                   & JJ_beta,
   		vector < SCIP_Real >                        & zval,
   		vector < SCIP_Real >                          x_lb,
   		vector < SCIP_Real >                          x_ub,
   		SCIP_Real                                    z_lb,
   		SCIP_Real                                    z_ub,
   		int                                          NNN,
   		int                                          j,
   		int                                          t,
   		SCIP_Real                                    y_ub
   		);

   SCIP_RETCODE Useall (
   		vector < int >                & KK1,
   		vector < int >                & JJ1,
   		vector < int >                  kk_max1,
   		vector < int >                  jj_max1);

   SCIP_RETCODE Usesome(
   		int                                          i,
   		vector < int >                             & KK1,
   		vector < int >                             & JJ1,
   		vector < int >                               kk_max1,
   		vector < int >                               jj_max1,
   		vector < vector < vector < SCIP_Real > > >   ycoeff	);

   SCIP_Real Selectrandom (
   		int                                            i,
   		vector < int >                               & KK1,
   		vector < int >                               & JJ1,
   		vector < int >                               kk_max1,
   		vector < int >                               jj_max1,
   		vector < vector < vector < SCIP_Real > > >   ycoeff	);

   SCIP_Real Choosetouse(
   		int                                          nzero,
   		int                                          i,
   		SCIP_Real                                    zcoeff,
   		vector < int >                             & KK1,
   		vector < int >                             & JJ1,
   		vector < int >                               kk_max1,
   		vector < int >                               jj_max1,
   		vector < vector < vector < SCIP_Real > > >   ycoeff);


private:

   //vector < vector < SCIP* > >                                        _subscip;
   //vector < vector < SCIP_VAR* > >                                    _p_zvar;
   //vector < vector < SCIP_Real > >                                    _p_zvar_lb;
   //vector < vector < SCIP_Real > >                                    _p_zvar_ub;
   //vector < vector < vector < SCIP_VAR* > > >                         _p_xvar; 
   //vector < vector < vector < SCIP_Real > > >                         _p_xvar_lb;
   //vector < vector < vector < SCIP_Real > > >                         _p_xvar_ub;
   //vector < vector < vector < vector < vector< SCIP_VAR* > > > > >   _p_yvar; 
   //vector < vector < vector < vector < vector< SCIP_Real > > > > >   _p_yvar_lb;
   //vector < vector < vector < vector < vector< SCIP_Real > > > > >   _p_yvar_ub;
   //vector < vector < vector < SCIP_CONS* > > >                        _p_wcon;             // Wastage constraint
   //vector < vector < vector < vector < SCIP_CONS* > > > >             _p_secon;     //Second constrain       
   int                                                                _Nper_to_spoil;
   double                                                             _Wastbeta;
   int                                                                _T;
   int                                                                _Nnodes;
   double                                                             _K;
   int                                                                _I;
   vector< double >                                                    _Nodeprobe;
   vector< double >                                                   _purchcost;
   vector < int >                                                     _vialsize;
   vector < vector < int > >                                          _jmat;  
   vector< vector<int > >                                             _nodemat;
   vector< SCIP_VAR* >                                                _z_var;
   vector< SCIP_VAR* >                                                _wj_var;
   vector< vector < vector <SCIP_CONS* > > >                          _sh_con;
   vector<vector<SCIP_CONS* > >                                       _z_con;
   SCIP_CONS*                                                         _alpha_con;
   vector< vector<SCIP_CONS* > >                                      _w_con;
   int                                                                _n_vars;
   bool                                                               _result;
};   
#endif

