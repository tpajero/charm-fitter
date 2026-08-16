#include <PDF_Charm.h>

#include <CharmParameters.h>

#include <RooArgList.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <TString.h>

void PDF_Charm::initParameters() {
  CharmParameters p;
  parameters = new RooArgList("parameters");
  for (const auto& name : getParameterNames()) parameters->add(*(p.get(name)));
}

void PDF_Charm::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *observables, *theory, covMatrix);
}

void PDF_Charm::initialise(const TString val_id, const TString unc_id, const TString cor_id, const bool buildCov) {
  initParameters();
  initRelations();
  initObservables();
  setObservables(val_id);
  setUncertainties(unc_id);
  setCorrelations(cor_id);
  if (buildCov)
    build();
  else
    buildPdf();
}
