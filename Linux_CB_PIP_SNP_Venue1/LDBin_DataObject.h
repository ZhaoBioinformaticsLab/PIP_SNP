/*Author: Wenchao Zhang, Noble Research Institute
  Written on 12/13/2018                          */
#ifndef LDBIN_DATAOBJECT_H_INCLUDED
#define LDBIN_DATAOBJECT_H_INCLUDED
#include <vector>
#include  "Common.h"
using namespace std;
class CLDBin_DataObject
{

public:
    int   LDBin_Individual_Num;    // Initialized for the Individual Number
    int   LDBin_Method_Corr;       // Initialized for Correlation Method
    int   LDBin_Start;
    int   LDBin_Current;
    float R2_LD_Left2;             // Record the R2 Value between the two most-left SNP markers. R2(SNP1,SNP2)
    float R2_LD_Right2;            // Record the R2 Value between the two most-right SNP markers. R2(SNP_Right-1,SNP_Right)
    float Sum_R2_Left_To_Current;  // Record the R2 Values between the most-left SNP marker to all other SNP markers. R2(SNP1,SNP2)+R2(SNP1,SNP3)+...+R2(SNP1,SNP_Right)
    float Sum_R2_Neighbor2;        // Record the R2 Values between the two neighbor SNP markers of the LDBin. R2(SNP1,SNP2)+R2(SNP2,SNP3)+...+R2(SNP_Right-1,SNP_Right)
    float R2_LD_Left_Right;        // Record the R2 Value between the most-left SNP marker and the most right SNP marker. R2(SNP1,SNP_Right)

    vector <char *> LDBin_marker_lines;

    int   Get_markerlines_number();
    void  Reset_LDBin_DataObject(int Line_Count, char * &markerline);
    void  Add_one_markerline(int Line_Count, char * &markerline, float R2_Left_Breakthrough, float R2_Right_Breakthorugh);
    float Calculate_R2_Left_one_markerline(char *markerline);   //Calculate the R2 value between the external markerline and the Left marker line of the LDMarker Group.
    float Calculate_R2_Right_one_markerline(char *markerline);  //Calculate the R2 value between the external markerline and the Right marker line of the LDMarker Group.

    void Clean_LDBin_marker_lines();
    CLDBin_DataObject(int Individual_Num, int Method_Corr, int Line_Count);

    ~CLDBin_DataObject(void);

};

#endif // LDBIN_DATAOBJECT_H_INCLUDED
