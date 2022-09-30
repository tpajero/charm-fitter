/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#ifndef PDF_BES_Kpi_h
#define PDF_BES_Kpi_h

#include "TString.h"

#include "ParametersCharmCombo.h"
#include "PDF_Abs.h"

class PDF_BES_Kpi : public PDF_Abs {
    public:
        PDF_BES_Kpi(const theory_config& th_cfg);
        ~PDF_BES_Kpi();
        void buildPdf() override;
        void initObservables() override;
        void initParameters() override;
        void initRelations() override;
        void setCorrelations(TString c) override;
        void setObservables(TString c) override;
        void setUncertainties(TString c) override;
    private:
        const theory_config th_cfg;
};

#endif
