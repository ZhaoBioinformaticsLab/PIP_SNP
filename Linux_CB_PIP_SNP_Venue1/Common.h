/*Author: Wenchao Zhang, Noble Research Institute
  Written on 12/13/2018                          */

#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#define ADD_CONSTANT 1
#include <vector>     // std vector
#include <algorithm>  // std:: sort
using namespace std;  // ignoring the std::


enum METTHOD_CORR {Pairwise_PCC=0, Pairwise_LD =1};                  // Correlation  methods to measure the similarities of two SNP markers, either as (1) LD based D value , or (2) Pearson Cross-Correlation
enum METHOD_DETECTION {Right_Breakthorugh =0, Left_Breakthorugh =1, RandL_Breakthrough =2, RorL_Breakthrough =3};  // Detection methods to define a LD Bin is ready. Right_Breakthorugh (1) to compare R2(SNP_Breakthrough,SNP_Right) with the user configured parameter; Left_Breakthrough (2) to compare R2(SNP_Breakthrough,SNP_Left)

/*Define a Genotype Position structure to represent the position */
struct Genotype_Positon
{
    int Marker_Pos;
    int Individual_Pos;
};

struct Individual_Pair_Distance_KnownGenotype
{
    int   X_Individual;
    int   Y_Individual;
    float Distance;
    char  KnownGenotype;
};

struct LDBin_Map_Feature
{
    int   LDBin_Start;             // Index recorded the start of the constructed LDBin, the order index in the orignal raw SNP
    int   LDBin_End;               // Index recorded the end of the constructed LDBin, the order index in the orignal raw SNP
    int   LDBin_Pos;               // Index of the representative position of the constructed LDBin, the order index in the orignal raw SNP
    float R2_LD_Left2;             // Recorded the R2 value between the most left two SNP markers
    float R2_LD_Right2;            // Record the R2 Value between the two most-right SNP markers. R2(SNP_Right-1,SNP_Right)
    float R2_LD_Left_Right;        // Recorded the R2 value between the most left and left SNP markers.
    float R2_LD_Average_Forward;   // Recorded the average of R2 values of the start SNP Marker with all other SNP marker at the detected LD Bin
    float R2_LD_Average_Neighbor2; // Recorded the average of R2 values of all the two neighbor SNP marker at the detected LD Bin
    float R2_Right_Breakthrough;   // Recorded the R2 Value of at the  boundary between the breakthrough and the right point
    float R2_Left_Breakthrough;    // Recorded the R2 Value between the breakthrough and the left point
};

bool Comp_Individual_Pair_by_Distance(Individual_Pair_Distance_KnownGenotype val, Individual_Pair_Distance_KnownGenotype p);

void destroy_grouplines(std::vector <char *> &group_marker_lines); //void destroy_grouplines(std::vector <char *marker_line> &group_marker_lines);

float Calculate_r2_one_line_to_another_line(char *line_one, char *line_two, int Individual_Num, int Method_Corr);
/* Calculate the R2 values between the current line to the group marker line */
void Calculate_r2_one_line_to_grouplines(vector <char *> group_marker_lines, char *line_one, int Individual_Num, int Method_Corr, float *R2_Average, float *R2_Boundary);

float Calculate_Pairwise_LD(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num);
float Calculate_Pairwise_PCC(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num);

#endif // COMMON_H_INCLUDED
