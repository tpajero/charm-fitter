/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_DY_RS_h
#define PDF_DY_RS_h

#include "TString.h"

#include "CharmUtils.h"
#include "PDF_Abs.h"

class PDF_DY_RS : public PDF_Abs
{
    public:
        PDF_DY_RS(TString measurement_id, const theory_config& th_cf);
        ~PDF_DY_RS();
        void buildPdf() override;
        void initObservables(const TString& setName);
        virtual void initParameters() override;
        virtual void initRelations() override;
        void setCorrelations(TString measurement_id) override;
        void setObservables(TString measurement_id) override;
        void setUncertainties(TString measurement_id) override;
    private:
        const theory_config th_cfg;
};

#endif

