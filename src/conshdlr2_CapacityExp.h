/*@file   Conshd_vaccinebp.h
 *@brief  C++ constraint handler 
 *@author Bahareh Eghtesadi*/

#ifndef __SCIP_CONSHDLRSEC_VACCINEBP_H__
#define __SCIP_CONSHDLRSEC_VACCINEBP_H__

#include "objscip/objscip.h"
#include "probdata_CapacityExp.h"

using namespace std;
using namespace scip;

class Objconshdlrvaccinebpsec : public scip::ObjConshdlr
{

public:
    
   Objconshdlrvaccinebpsec(
            SCIP*                 scip,
            const char *          name);

  virtual ~Objconshdlrvaccinebpsec();


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

};
#endif
