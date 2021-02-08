/*Author: Wenchao Zhang, Noble Research Institute
  Written on 02/24/2020                          */
#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <vector>     // std vector
#include <algorithm>  // std:: sort
#define Min(x, y) (((x) < (y)) ? (x) : (y))
#define Max(x, y) (((x) < (y)) ? (y) : (x))
#define Abs(x)    (((x) > 0) ? (x) : (-x))
#define Sum(x,y)  (x+y)
#define ABS_CORR  1

using namespace std;  // ignoring the std::
enum METTHOD_CORR {Pairwise_PCC=0, Pairwise_LD =1};
//enum METTHOD_SYNTHESIS {Not_Synthesis=0, Norm_Synthesis=1, Representative_Synthesis =2};
enum METTHOD_SYNTHESIS {Norm_Synthesis=1, Representative_Synthesis =2};


struct LDBin_Map_Feature
{
    int   LDBin_Start;             // Index recorded the start of the constructed LDBin, the order index in the orignal raw SNP
    int   LDBin_End;               // Index recorded the end of the constructed LDBin, the order index in the orignal raw SNP
    int   LDBin_Pos;               // Index of the representative position of the constructed LDBin, the order index in the orignal raw SNP
    float R2_Breakpoint;           // The R2 Vaule at the break point
};

struct LDBin_Object
{
    int   LDBin_Start;             // Index recorded the start of the constructed LDBin, the order index in the orignal raw SNP
    int   LDBin_End;               // Index recorded the end of the constructed LDBin, the order index in the orignal raw SNP
    int   LDBin_Pos;               // Index of the representative position of the constructed LDBin, the order index in the orignal raw SNP
    float *LDBinData;
};
bool compareLDBin(LDBin_Object val, LDBin_Object p);
float Calculate_Pairwise_LD(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num, int ABS_Mode);
float Calculate_Pairwise_LD(float *Genotype_Array1, float *Genotype_Array2, int Individual_Num, int ABS_Mode);
float Calculate_Pairwise_PCC(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num, int ABS_Mode);   // Use the overide function to accept the representative char type genotype as 0,1,2
float Calculate_Pairwise_PCC(float *Genotype_Array1, float *Genotype_Array2, int Individual_Num, int ABS_Mode); // Use the overide function to accept the synthesised float type genotype

float Calculate_r2_one_bin_to_another_bin(float *bin_one, float *bin_two, int Individual_Num, int method_corr);
void calculate_pairwise_corr_sum(float *pLDBin_sum_corr, float *pdata_Buffer, int Individual_Num, int Bin_Num, int method_corr);

void Synthesis_Process(float *pdata_Buffer, float* pdata_Synthesised, int Individual_Num, int Bin_Num, int *Pos, int LDBin_Corr_Method, int LDBin_Synthesis_Method); //

#endif // COMMON_H_INCLUDED

