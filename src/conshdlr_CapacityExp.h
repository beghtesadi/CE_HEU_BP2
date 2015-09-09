/*@file   Conshd_vaccinebp.h
 *@brief  C++ constraint handler 
 *@author Bahareh Eghtesadi*/

#ifndef __SCIP_CONSHDLR_VACCINEBP_H__
#define __SCIP_CONSHDLR_VACCINEBP_H__

#include "objscip/objscip.h"
#include "probdata_CapacityExp.h"

using namespace std;
using namespace scip;

class Objconshdlrvaccinebp : public scip::ObjConshdlr
{

public:
    
   Objconshdlrvaccinebp(
            SCIP*                 scip,
            const char *          name);

  virtual ~Objconshdlrvaccinebp();

  SCIP_RETCODE consdataCreate(
                SCIP*                scip,
                SCIP_CONSDATA**      consdata,  
                SCIP_NODE*           childnode, 
                SCIP_Bool            candtype, 
                SCIP_Real            fracval, 
                bool                 ztype,
                int                  stage, 
                int                  state, 
                int                  nodetype,
                int                  vialtype
   );

  SCIP_RETCODE checkVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR*             var,                /**< variables to check  */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   );

  SCIP_RETCODE consdataFixVariables(
   SCIP*                          scip,               /**< SCIP data structure */
   SCIP_CONSDATA*                 consdata,           /**< constraint data */
   vector< SCIP_VAR*>             vars,               /**< generated variables */
   int                            nvars,              /**< number of generated variables */
   SCIP_RESULT*                   result              /**< pointer to store the result of the fixing */
   );

  SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the constraint data */
   );

  virtual SCIP_DECL_CONSDELETE(scip_delete);

  virtual SCIP_DECL_CONSTRANS(scip_trans);

  virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

  virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);

  virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

  virtual SCIP_DECL_CONSENFOPS(scip_enfops);

  virtual SCIP_DECL_CONSCHECK(scip_check);

  virtual SCIP_DECL_CONSLOCK(scip_lock);

  virtual SCIP_DECL_CONSPROP(scip_prop);

  virtual SCIP_DECL_CONSACTIVE(scip_active);

  virtual SCIP_DECL_CONSDEACTIVE(scip_deactive);

  SCIP_RETCODE SCIPcreateconshdlr(
                SCIP*                scip, 
                SCIP_CONS**          cons,               /**< pointer to hold the created constraint */
                const char*          name,  
                SCIP_NODE*           childnode, 
                SCIP_Bool            candtype, 
                SCIP_Real            fracval, 
                bool                 ztype,
                int                  stage, 
                int                  state, 
                int                  nodetype,
                int                  vialtype);

int SCIPgetnodetype(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

int SCIPgetstage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

int SCIPgetstate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

bool SCIPgetztype(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

int SCIPgetcandtype(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );

SCIP_Real SCIPgetcandfrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );
int SCIPgetvialtype(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   );
};
#endif
