/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#ifndef PDF_K3pi_h
#define PDF_K3pi_h

#include "TString.h"

#include "CharmUtils.h"
#include "PDF_Abs.h"

class PDF_K3pi : public PDF_Abs
{
    public:
        PDF_K3pi(TString measurement_id, const theory_config& th_cfg);
        ~PDF_K3pi();
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
