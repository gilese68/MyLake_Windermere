#include <iostream>
#include "mex.h"
#include <IPhreeqc.hpp>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

    size_t n_species;
    double *rows = (double *) mxGetData(prhs[0]);
    double *data = (double *) mxGetData(prhs[1]);
    int n_rows = (int) rows[0];
    float alkalinity;
    float C4;
    float N5;
    float S2;
    float S6;
    float Fe2;
    float Fe3;
    float Ca;
    float P;
    float N3;
    char sol_str[256];
    char alkalinity_str[256];
    char C4_str[256];
    char N5_str[256];
    char N3_str[256];
    char S2_str[256];
    char S6_str[256];
    char Fe2_str[256];
    char Fe3_str[256];
    char Ca_str[256];
    char P_str[256];
    char K_str[256];
    char Na_str[256];
    char Al_str[256];
    char Mg_str[256];
    char eq_str[256];
    char cell_str[256];
    double H0[n_rows];
    double HCO30[n_rows];
    double CO20[n_rows];
    double CO30[n_rows];
    double NH30[n_rows];
    double NH40[n_rows];
    double HS0[n_rows];
    double H2S0[n_rows];
    double OH0[n_rows];
    double CO2g0[n_rows];
    double Fe20[n_rows];
    double Ca20[n_rows];
    double NO30[n_rows];
    double SO40[n_rows];
    double PO40[n_rows];
    double FeS0[n_rows];
    double FeS20[n_rows];
    double FeOH30[n_rows];
    double FeOOH0[n_rows];
    double Ca3PO40[n_rows];
    double PO4adsa0[n_rows];
    double PO4adsb0[n_rows];

    double est_pH[n_rows];

    n_species = mxGetN(prhs[1]);

    for (int i = 0; i < n_rows; ++i)    {
        H0[i]       = data[ 0 * n_rows + i];
        HCO30[i]    = data[ 1 * n_rows + i];
        CO20[i]     = data[ 2 * n_rows + i];
        CO30[i]     = data[ 3 * n_rows + i];
        NH30[i]     = data[ 4 * n_rows + i];
        NH40[i]     = data[ 5 * n_rows + i];
        HS0[i]      = data[ 6 * n_rows + i];
        H2S0[i]     = data[ 7 * n_rows + i];
        OH0[i]      = data[ 8 * n_rows + i];
        CO2g0[i]   = data[ 9 * n_rows + i];
        Fe20[i]     = data[10 * n_rows + i];
        Ca20[i]     = data[11 * n_rows + i];
        NO30[i]     = data[12 * n_rows + i];
        SO40[i]     = data[13 * n_rows + i];
        PO40[i]     = data[14 * n_rows + i];
        FeS0[i]     = data[15 * n_rows + i];
        FeS20[i]    = data[16 * n_rows + i];
        FeOH30[i]   = data[17 * n_rows + i];
        FeOOH0[i]   = data[18 * n_rows + i];
        Ca3PO40[i]  = data[19 * n_rows + i];
        PO4adsa0[i] = data[20 * n_rows + i];
        PO4adsb0[i] = data[21 * n_rows + i];
    }


    IPhreeqc iphreeqc;

    if (iphreeqc.LoadDatabase("Sediment-v2.0/ph-module/phreeqc.dat") != 0 && iphreeqc.LoadDatabase("phreeqc.dat") != 0) {
        mexPrintf("%s", iphreeqc.GetErrorString());
    }

    /*
        SOLUTION 1
            -temp 8
            -units mmol/kgw
            Alkalinity 0.3 ## = [HCO3] + 2[CO3] + [OH] - [NH40] + 1.5[PO4] - [H+]
            C(+4)   0.1 ## = HCO3 + CO2aq +CO3
            N(+5) 0.1 ## = NO3
            S(-2)   0.1 ## = HS + H2S + 0.5 * FeS2
            Fe(+2)  0.1 ## = Fe2 + FeS2 + FeS
            Fe(+3)  0.1 ## = FeOH3 + FeOOH
            Ca  0.1 ## = Ca2 + Ca3PO4
            P 0.1 ## = PO4 + Ca3PO4 + PO4adsa + PO4adsb
        EQUILIBRIUM_PHASES 1
            Pyrite ## = FeS2
            Mackinawite ## =FeS
            Goethite ## =FeOOH
            Fe(OH)3(a) ## =FeOH3
            Vivianite ##


        SELECTED_OUTPUT
            -molalities Fe+2
        END
    */

    for (int i = 0; i < n_rows; ++i) {
        alkalinity = HCO30[i] + 2 * CO30[i] + OH0[i]  - H0[i] + HS0[i] - NH40[i];
        C4 =  CO30[i] + HCO30[i];// CO20[i];
        N5 = NO30[i];
        S2 = HS0[i] + H2S0[i]; // + 0.5 * FeS20[i];
        S6 = SO40[i]; // + 0.5 * FeS20[i];
        Fe2 = Fe20[i]; // + FeS20[i] +  FeS0[i]
        // Fe3 = 0; // FeOH30[i] + FeOOH0[i]
        Ca = Ca20[i]; // + Ca3PO40[i]
        P = PO40[i];  // + Ca3PO40[i]
        N3 = NH40[i] + NH30[i];

        sprintf(sol_str,        "SOLUTION %i", i + 1);
        sprintf(alkalinity_str, "     Alkalinity    %f", alkalinity);
        sprintf(C4_str,         "     C(+4)         %f", C4);
        sprintf(N5_str,         "     N(+5)         %f", N5);
        sprintf(N3_str,         "     NH4+         %f", N3);
        sprintf(S2_str,         "     S(-2)         %f", S2);
        sprintf(S6_str,         "     SO4         %f", S6);
        sprintf(Fe2_str,        "     Fe(+2)        %f", Fe2);
        // sprintf(Fe3_str,        "     Fe(+3)        %f", Fe3);
        sprintf(Ca_str,         "     Ca            %f", Ca);
        sprintf(P_str,          "     P             %f", P);
        // sprintf(K_str,          "     K             %f", 0.12);
        // sprintf(Na_str,          "     Na             %f", 1.6);
        // sprintf(Al_str,          "     Al             %f", 0.013);
        // sprintf(Mg_str,          "     Mg             %f", 0.14);

        iphreeqc.AccumulateLine(sol_str);
        iphreeqc.AccumulateLine("     temperature 8");
        iphreeqc.AccumulateLine("     pressure 5");
        iphreeqc.AccumulateLine("     -units mmol/kgw");
        iphreeqc.AccumulateLine(alkalinity_str);
        iphreeqc.AccumulateLine(C4_str);
        iphreeqc.AccumulateLine(N5_str);
        iphreeqc.AccumulateLine(N3_str);
        iphreeqc.AccumulateLine(S2_str);
        iphreeqc.AccumulateLine(Fe2_str);
        // iphreeqc.AccumulateLine(Fe3_str);
        iphreeqc.AccumulateLine(Ca_str);
        iphreeqc.AccumulateLine(P_str);
        // iphreeqc.AccumulateLine(K_str);
        // iphreeqc.AccumulateLine(Na_str);
        // iphreeqc.AccumulateLine(Al_str);
        // iphreeqc.AccumulateLine(Mg_str);
        // iphreeqc.AccumulateLine("END");
    }
    // for (int i = 0; i < n_rows; ++i) {
    // sprintf(eq_str,        "EQUILIBRIUM_PHASES %i", i+1);
    // iphreeqc.AccumulateLine(eq_str);
    // }

    // iphreeqc.AccumulateLine("END");
    iphreeqc.AccumulateLine("SELECTED_OUTPUT");
    iphreeqc.AccumulateLine("     -molalities Fe+2");
    // iphreeqc.AccumulateLine("END");
    // iphreeqc.AccumulateLine("RUN_CELLS");
    // sprintf(cell_str,       "     -cells 1-%i",n_rows);
    // iphreeqc.AccumulateLine(cell_str);
    iphreeqc.AccumulateLine("END");


    // show error if not converge ?
    if (iphreeqc.RunAccumulated() != 0) {
        mexPrintf("Phreeqc::Sediments::%s", iphreeqc.GetErrorString());
    }

    VAR v;
    ::VarInit(&v);
    // for (int i = 0; i < iphreeqc.GetSelectedOutputRowCount(); ++i) {
    //   for (int j = 0; j < iphreeqc.GetSelectedOutputColumnCount(); ++j) {
    //     if (iphreeqc.GetSelectedOutputValue(i, j, &v) == VR_OK) {
    //       switch (v.type) {
    //       case TT_LONG:
    //         // std::cout << v.lVal << " ";
    //         mexPrintf("%f\t", v.lVal);
    //         break;
    //       case TT_DOUBLE:
    //         // std::cout << v.dVal << " ";
    //         mexPrintf("%f\t", v.dVal);
    //         break;
    //       case TT_STRING:
    //         // std::cout << v.sVal << " ";
    //         mexPrintf("%s\t", v.sVal);
    //         break;
    //       }
    //     }
    //     ::VarClear(&v);
    //   }
    //   mexPrintf("\n");
    // }

    for (int i = 0; i < n_rows; ++i) {
        iphreeqc.GetSelectedOutputValue(i + 1, 6, &v) ;
        est_pH[i] = v.dVal;
        ::VarClear(&v);
        // mexPrintf("%f\n", est_pH[i]);
    }




    // Debugging
    // mexPrintf("=%f\n", data[63]);
    // mexPrintf("Number of elements = %i\n", n_species);
    // mexPrintf("Size of element = %i\n",n_rows);
    // mexPrintf("Element = %f\n",Ca3PO40(0));
    // mexPrintf("Element = %d\n", v.dVal);



    plhs[0] = mxCreateDoubleMatrix(1, n_rows, mxREAL);
    double * z = mxGetPr(plhs[0]);
    for (int i = 0; i < n_rows; ++i) {
        z[i] = est_pH[i];
    }
}
