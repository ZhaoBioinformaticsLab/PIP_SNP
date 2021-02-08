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


void CLDBin_DataObject::Reset_LDBin_DataObject(int Line_Strat, int Line_End)
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

    LDBin_Start             =Line_Strat;
    LDBin_Current           =Line_End;

}

void CLDBin_DataObject::Add_one_markerline(int Line_Count, char * &markerline)
{
    LDBin_Current           =Line_Count;
    if(LDBin_Current <LDBin_Start)
    {
        cerr<< "LD Bin Start can not greater than the LD Current, Error at: LDBin_Start=" << LDBin_Start << "LDBin_Current =" << Line_Count<<endl;
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

CLDBin_DataObject::CLDBin_DataObject(int Individual_Num, int Method_Corr, int Line_Strat, int Line_End)
{
    LDBin_Individual_Num    =Individual_Num;
    LDBin_Method_Corr       =Method_Corr;
    LDBin_Start             =Line_Strat;
    LDBin_Current           =Line_End;

}

CLDBin_DataObject::~CLDBin_DataObject(void)
{
    while(!LDBin_marker_lines.empty())
    {
       delete LDBin_marker_lines.back();
       LDBin_marker_lines.pop_back();
    }
}
