/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_BinFlip_h
#define PDF_BinFlip_h

#include "TString.h"

#include "ParametersCharmCombo.h"
#include "PDF_Abs.h"

class PDF_BinFlip : public PDF_Abs
{
    public:
        PDF_BinFlip(TString measurement_id, const theory_config& th_cf);
        ~PDF_BinFlip();
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
