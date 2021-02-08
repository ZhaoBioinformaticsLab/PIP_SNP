/*Author: Wenchao Zhang, Noble Research Institute
  Written on 12/13/2018                          */

#include <vector>   // For using Vector
#include <math.h>
#include  "LDBin_DataObject.h"
#include <iostream>  // Cin/ Cout/Cerr
using namespace std;
int CLDBin_DataObject::Get_markerlines_number()
{
    return LDBin_marker_lines.size();
}


void CLDBin_DataObject::Reset_LDBin_DataObject(int Line_Count, char * &markerline)
{
    /* Clean the data buffer*/
    while(!LDBin_marker_lines.empty())
    {

       //delete LDBin_marker_lines.back();
       char *current_line= LDBin_marker_lines.back();
       if(current_line==NULL )
       {
           cerr <<"Null Pointer at Reset LDBin" <<endl;
       }

       delete current_line;
       LDBin_marker_lines.pop_back();
    }

    /* Reset the object property variables*/
    if(markerline !=NULL)
    {
        LDBin_marker_lines.push_back(markerline);
    }

    LDBin_Start             =Line_Count;
    LDBin_Current           =Line_Count;

    R2_LD_Left2             =0.0;
    R2_LD_Right2            =0.0;
    R2_LD_Left_Right        =0.0;
    Sum_R2_Left_To_Current  =0.0;
    Sum_R2_Neighbor2        =0.0;

}

void CLDBin_DataObject::Add_one_markerline(int Line_Count, char * &markerline, float R2_Left_Breakthrough, float R2_Right_Breakthrough)
{
    LDBin_Current           =Line_Count;
    if(LDBin_Current <LDBin_Start)
    {
        cerr<< "LD Bin Start can not greater than the LD Current, Error at: LDBin_Start=" << LDBin_Start << "LDBin_Current =" << Line_Count<<endl;
    }
    else
    {
        if(LDBin_Current >LDBin_Start) //Indicate A new marker line has been added, make update the R2 sum value
        {

            R2_LD_Left_Right                       = R2_Left_Breakthrough;
            Sum_R2_Left_To_Current                += R2_Left_Breakthrough;

            if((LDBin_Current -LDBin_Start) ==1)
            {
               R2_LD_Left2 = R2_Left_Breakthrough;
            }
            R2_LD_Right2                           =R2_Right_Breakthrough;
            Sum_R2_Neighbor2                      +=R2_Right_Breakthrough;

        }
    }

    if(markerline !=NULL)
    {
        LDBin_marker_lines.push_back(markerline);
    }


}

//Calculate the R2 value between the external markerline and the Left marker line of the LDMarker Group.
float CLDBin_DataObject::Calculate_R2_Left_one_markerline(char *markerline)
{
    vector<char *>::iterator it_firsit_line = LDBin_marker_lines.begin();
    float R2                                = Calculate_r2_one_line_to_another_line(*it_firsit_line, markerline, LDBin_Individual_Num, LDBin_Method_Corr);
    return R2;
}

//Calculate the R2 value between the external markerline and the Right marker line of the LDMarker Group.
float CLDBin_DataObject::Calculate_R2_Right_one_markerline(char *markerline)
{
    char * Right_line                       = LDBin_marker_lines.back();
    float R2                                = Calculate_r2_one_line_to_another_line(Right_line, markerline, LDBin_Individual_Num, LDBin_Method_Corr);
    return R2;
}

void CLDBin_DataObject::Clean_LDBin_marker_lines()
{
    while(!LDBin_marker_lines.empty())
    {
       delete LDBin_marker_lines.back();
       LDBin_marker_lines.pop_back();
    }
}

CLDBin_DataObject::CLDBin_DataObject(int Individual_Num, int Method_Corr, int Line_Count)
{
    LDBin_Individual_Num    =Individual_Num;
    LDBin_Method_Corr       =Method_Corr;
    LDBin_Start             =Line_Count;
    LDBin_Current           =Line_Count;

    R2_LD_Left2             =0.0;
    R2_LD_Right2            =0.0;
    R2_LD_Left_Right        =0.0;
    Sum_R2_Left_To_Current  =0.0;
    Sum_R2_Neighbor2        =0.0;

}

CLDBin_DataObject::~CLDBin_DataObject(void)
{
    while(!LDBin_marker_lines.empty())
    {
       delete LDBin_marker_lines.back();
       LDBin_marker_lines.pop_back();
    }
}
