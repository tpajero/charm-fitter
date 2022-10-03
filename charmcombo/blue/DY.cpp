/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2022
 *
 * This code is based on the BLUE package (https://blue.hepforge.org/).
 *
 * Run it from a root interactive session through:
 * 
 *    gSystem->Load("libBlue.so")
 *    .L DY.cpp
 *    DY(<n-combination>)
 *
 **/

#include "Blue.h"

void DY(Int_t flag = 0) {
  
    //----------------------------------------------------------------------------
    // Flag steers which of the results should be calculated
    //  0: The LHCb Run 1-2 combination
    //  1: The LHCb Run 1-2 combination for KK
    //  2: The LHCb Run 1-2 combination for PP
    //  3: The LHCb Run 2 D*+-tagged combination
    //  4: The world average
    //  5: The world average for KK (LHCb Run1-2 + CDF)
    //  6: The world average for PP (LHCb Run1-2 + CDF)
    //  7: The world average for KK as of 2020 (LHCb Run1-2mu + CDF)
    //  8: The world average for PP as of 2020 (LHCb Run1-2mu + CDF)
    //  9: The average of the Belle and BaBar measurements
    //----------------------------------------------------------------------------

    // number of estimates, uncertainties, observables
    static const Int_t num_est = 12;
    static const Int_t num_unc = 10;
    static const Int_t num_obs = 1;
    Int_t IWhichObs[num_est] = {num_est*0};

    static const TString names[num_est] = {
            "BaBar         ",  //  0
            "CDF_KK        ",  //  1
            "CDF_PP        ",  //  2
            "Belle         ",  //  3
            "Run1_mu_KK    ",  //  4
            "Run1_mu_PP    ",  //  5
            "Run1_prompt_KK",  //  6
            "Run1_prompt_PP",  //  7
            "Run2_mu_KK    ",  //  8
            "Run2_mu_PP    ",  //  9
            "Run2_prompt_KK",  // 10 
            "Run2_prompt_PP"}; // 11
    static const TString names_unc[num_unc] = {
            "  Stat",
            "m(hhh)",
            "   Sec",
            " m(KK)",
            " m(PP)",
            "Weight",
            "TimRes",
            "Mistag",
            "Run1Mu",
            " Other"};
    static const TString names_obs[num_obs] = {"     DY"};
  

    // array of estimates and uncertainties
    static const Int_t len_x_est = num_est * (num_unc+1);
    Double_t x_est[len_x_est] = {
            //             0      1      2      3      4       5      6      7     8       9
            //          Stat m(hhh)    Sec   m(KK)  m(PP) Weight TimeRes Mistag Run1Mu Other
             -8.8 ,    25.5 ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.  ,   0.  ,  0.  , 5.8 ,  // BaBar
             19.  ,    15.  ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.  ,   0.  ,  0.  , 0.  ,  // CDF KK
              1.  ,    18.  ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.  ,   0.  ,  0.  , 0.  ,  // CDF PP
              3.  ,    20.  ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.  ,   0.  ,  0.  , 7.  ,  // Belle
             13.40,     7.70,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.49,   0.  ,  2.55, 0.  ,  // Run 1 mu KK
              9.2 ,    14.5 ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.42,   0.  ,  2.48, 0.  ,  // Run 1 mu PP
              3.0 ,     3.2 ,  0.1 ,  0.8 ,  0.5 ,  0.  ,   0.2 ,  0.  ,   0.  ,  0.  , 0.  ,  // Run 1 prompt KK
             -4.6 ,     5.8 ,  0.1 ,  1.2 ,  0.  ,  0.  ,   0.2 ,  0.  ,   0.  ,  0.  , 0.  ,  // Run 1 prompt PP
          //  4.3 ,     3.6 ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.3 ,   0.3 ,  0.  , 0.3 ,  // Run 2 mu KK
          // -2.2 ,     7.0 ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.4 ,   0.6 ,  0.  , 0.3 ,  // Run 2 mu PP
              4.8 ,     4.0 ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.  ,   0.  ,  0.  , 0.3 ,  // Run 2 mu KK  (scale factor of 1.12 to account for dilution)
             -2.5 ,     7.8 ,  0.  ,  0.  ,  0.  ,  0.  ,   0.  ,  0.  ,   0.  ,  0.  , 0.3 ,  // Run 2 mu PP  (scale factor of 1.12 to account for dilution)
             -2.321,    1.524, 0.24,  0.13,  0.06,  0.  ,   0.05,  0.  ,   0.  ,  0.  , 0.  ,  // Run 2 prompt KK
             -4.014,    2.814, 0.34,  0.13,  0.  ,  0.03,   0.05,  0.  ,   0.  ,  0.  , 0.  ,  // Run 2 prompt PP
    };
    static const Double_t rho_val[num_unc] = {
                        0.,    1.,    1.,    1.,    1.,     1.,    1.,     1.,    1.,   0.};

    TMatrixD *inp_est = new TMatrixD(num_est, num_unc+1, &x_est[0]);

    // Statistical precision in systematic uncertainties
    static const Int_t len_s_unc = num_est * num_unc;
    static const Double_t s_unc[len_s_unc] = {
            // Stat
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // BaBar
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // CDF KK
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // CDF PP
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Belle
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 1 mu KK
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 1 mu PP
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 1 prompt KK
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 1 prompt PP
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 2 mu KK
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 2 mu PP
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 2 prompt KK
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  // Run 2 prompt PP
    };
    TMatrixD *inp_sta = new TMatrixD(num_est, num_unc, &s_unc[0]);


    static const TString for_val = "%2.2f";
    static const TString for_unc = "%2.2f";
    static const TString for_wei = "%2.2f";
    static const TString for_rho = "%2.2f";
    static const TString for_pul = "%2.2f";
    static const TString for_chi = "%2.2f";
    static const TString for_uni = "";

    Blue *my_blue = new Blue(num_est, num_unc);
    my_blue->SetFormat(for_val, for_unc, for_wei, for_rho, for_pul, for_chi, for_uni);
    my_blue->FillNamEst(&names[0]);
    my_blue->FillNamUnc(&names_unc[0]);
    my_blue->FillNamObs(&names_obs[0]);

    my_blue->FillEst(inp_est);
    my_blue->FillSta(inp_sta);
    for (Int_t i=0; i<num_unc; i++) {
        my_blue->FillCor(i, rho_val[i]);
    }

    switch(flag) {
        case 0:
            std::cout << "-------------------------------------------------\n"
                      << "LHCb Run 1+2 average\n"
                      << "-------------------------------------------------\n";
            for (auto i : {0, 1, 2, 3}) my_blue->SetInActiveEst(i);
            break;
        case 1:
            std::cout << "-------------------------------------------------\n"
                      << "LHCb Run 1+2:  D0 -> K+ K- average\n"
                      << "-------------------------------------------------\n";
            for (auto i : {0, 1, 2, 3, 5, 7, 9, 11}) my_blue->SetInActiveEst(i);
            break;
        case 2:
            std::cout << "-------------------------------------------------\n"
                      << "LHCb Run 1+2:  D0 -> pi+ pi- average\n"
                      << "-------------------------------------------------\n";
            for (auto i : {0, 1, 2, 3, 4, 6, 8, 10}) my_blue->SetInActiveEst(i);
            break;
        case 3:
            std::cout << "-------------------------------------------------\n"
                      << "LHCb Run 2 (D*+)-tagged (K+ K- & pi+ pi- average)\n"
                      << "-------------------------------------------------\n";
            for (int i=0; i<10; ++i) my_blue->SetInActiveEst(i);
            break;
        case 4:
            std::cout << "-------------------------------------------------\n"
                      << "World average\n"
                      << "-------------------------------------------------\n";
            break;
        case 5:
            std::cout << "-------------------------------------------------\n"
                      << "World average (K+ K-)\n"
                      << "-------------------------------------------------\n";
            for (auto i : {0, 2, 3, 5, 7, 9, 11}) my_blue->SetInActiveEst(i);
            break;
        case 6:
            std::cout << "-------------------------------------------------\n"
                      << "World average (pi+ pi-)\n"
                      << "-------------------------------------------------\n";
            for (auto i : {0, 1, 3, 4, 6, 8, 10}) my_blue->SetInActiveEst(i);
            break;
        case 7:
            std::cout << "-------------------------------------------------\n"
                      << "World average 2020 (K+ K-)\n"
                      << "-------------------------------------------------\n";
            for (auto i : {0, 2, 3, 5, 7, 9, 10, 11}) my_blue->SetInActiveEst(i);
            break;
        case 8:
            std::cout << "-------------------------------------------------\n"
                      << "World average 2020 (pi+ pi-)\n"
                      << "-------------------------------------------------\n";
            for (auto i : {0, 1, 3, 4, 6, 8, 10, 11}) my_blue->SetInActiveEst(i);
            break;
        case 9:
            std::cout << "-------------------------------------------------\n"
                      << "BaBar + Belle\n"
                      << "-------------------------------------------------\n";
            my_blue->SetInActiveEst(1);
            my_blue->SetInActiveEst(2);
            for (int i=4; i<12; ++i) my_blue->SetInActiveEst(i);
            break;
        default:
            std::cerr << "Unsupported flag value " << flag << ".\nExiting...\n";
            exit(EXIT_FAILURE);
    }

    my_blue->FixInp();
    my_blue->PrintEst();
    my_blue->Solve();
    my_blue->PrintResult();

    delete my_blue;
    inp_est->Delete();
    inp_sta->Delete();

    return;
}
