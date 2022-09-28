/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#ifndef PDF_DY_f_h
#define PDF_DY_f_h

#include "TString.h"

#include "ParametersCharmCombo.h"
#include "PDF_Abs.h"

class PDF_DY_f : public PDF_Abs
{
    public:
        PDF_DY_f(TString measurement_id, const theory_config& th_cf);
        ~PDF_DY_f();
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
