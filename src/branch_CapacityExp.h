
#ifndef __SCIP_BRANCH_VACCINEBP_H__
#define __SCIP_BRANCH_VACCINEBP_H__

#include <iostream>
#include "objscip/objscip.h"

using namespace std;
using namespace scip;


/** pricer class */
class ObjBranchrulevaccinebp : public ObjBranchrule
{
public:
    ObjBranchrulevaccinebp(
             SCIP *             scip,
             const char * 	name );

    virtual ~ObjBranchrulevaccinebp();

    virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

};

#endif
