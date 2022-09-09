/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#ifndef PDF_AcpHH_LHCb_Run12_h
#define PDF_AcpHH_LHCb_Run12_h

#include "TString.h"

#include "ParametersCharmCombo.h"
#include "PDF_Abs.h"


class PDF_AcpHH_LHCb_Run12 : public PDF_Abs
{
    public:
        PDF_AcpHH_LHCb_Run12(const theory_config& th_cfg);
        ~PDF_AcpHH_LHCb_Run12();
        void buildPdf() override;
        void initObservables() override;
        void initParameters() override;
        void initRelations() override;
        void setCorrelations(TString c) override;
        void setObservables(TString c) override;
        void setUncertainties(TString c) override;
    private:
        // Helper functions to avoid boilerplate code
        void add_acpkk(RooArgList* theory, TString name, double avg_time);
        void add_dacp(RooArgList* theory, TString name, double avg_time_kk, double avg_time_pipi);
        const double d0_lifetime = 4.103e-13;
        const theory_config th_cfg;
};

#endif
