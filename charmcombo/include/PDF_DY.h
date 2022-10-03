/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_DY_h
#define PDF_DY_h

#include "TString.h"

#include "CharmUtils.h"
#include "PDF_Abs.h"

class PDF_DY : public PDF_Abs {
    public:
        PDF_DY(TString measurement_id, const theory_config th_cf, const FSC fsc);
        ~PDF_DY();
        void buildPdf() override;
        void initObservables(const TString& setName);
        virtual void initParameters() override;
        virtual void initRelations() override;
        void setCorrelations(TString measurement_id) override;
        void setObservables(TString measurement_id) override;
        void setUncertainties(TString measurement_id) override;
    private:
        const FSC fsc;
        const theory_config th_cfg;
};

#endif
