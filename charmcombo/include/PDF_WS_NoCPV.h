/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_WSMixingNoCPV_h
#define PDF_WSMixingNoCPV_h

#include "TString.h"

#include "CharmUtils.h"
#include "PDF_Abs.h"

class PDF_WS_NoCPV : public PDF_Abs
{
    public:
        PDF_WS_NoCPV(TString measurement_id, const theory_config& th_cfg);
        ~PDF_WS_NoCPV();
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
