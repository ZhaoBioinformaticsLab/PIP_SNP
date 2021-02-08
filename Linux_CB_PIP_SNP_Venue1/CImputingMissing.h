/*Author: Wenchao Zhang, Noble Research Institute
  Written on 12/13/2018                          */
#ifndef CIMPUTINGMISSING_H
#define CIMPUTINGMISSING_H

#include  "Common.h"
#include  "LDBin_DataObject.h"
using namespace std;
class CImputingMissing
{
    public:
        CImputingMissing(CLDBin_DataObject &ldbin_obj, int K);
        virtual ~CImputingMissing();
    int   LDBin_Start;
    int   LDBin_End;
    int   LDBin_Pos;
    int   LDBin_Marker_Num;
    int   LDBin_Individual_Num;
    int   Method_Corr;

    int   K_NN;  //An integer parameter for KNN method, which is determine the only limited individuals for imputing the missing genotype
    char  *Lpldbin_OriginalData_Transpose;  // (nxm format) To store the original LDBin Data, a data pointer by cascading all the individual line vector.
    char  *Lpldbin_ImputedData;             // (mxn format) To store the imputed LDBin Data, as a pointer cascading one marker line with another marker line
    float *Lpldbin_Synthesized;             // (1xn format) To store the synthesized LDBin Data, either from the integration(usually as the norm) of multiple SNP markers, or the optimal SNP marker

    vector <struct Genotype_Positon> Missing_Genotype_Vector; // A Genotype Position Vectors to store all the missing genotype position information.

    float calcuate_distance_two_indivudals(char *data_individual_1, char *data_individual_2, int marker_len);
    char  Impute_one(int missing_individual_pos, int missing_marker_pos); // Impute one specific missing genotype at position of (marker_pos, individual_pos)
    void  Impute_all();
    void  calculate_pairwise_corr_sum(float *Lpldbin_marker_sum_corr, int method_corr);
    void  Synthesize(int syn_mode);  // 0: not to be synthesized; 1: norm2 synthesized; 2: to select the optimal marker

    protected:
    private:
};

#endif // CIMPUTINGMISSING_H
