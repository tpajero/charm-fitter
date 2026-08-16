#ifndef PDF_Charm_h
#define PDF_Charm_h

#include <PDF_Abs.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Utility class for charm PDFs, to avoid boilerplate code in each PDF class.
 */
class PDF_Charm : public PDF_Abs {
 public:
  using PDF_Abs::PDF_Abs;
  void initParameters() override;
  void buildPdf() override;

 protected:
  /**
   * Run the standard initialisation sequence needed by most PDFs.
   *
   * @param buildCov Whether to build the covariance matrix and PDF as part of this initialisation, or the PDF only.
   */
  void initialise(TString val_id, TString unc_id, TString cor_id, bool buildCov = true);

 private:
  /// Get the names of the parameters needed for the theory expressions.
  virtual std::set<std::string> getParameterNames() const = 0;
};

#endif
