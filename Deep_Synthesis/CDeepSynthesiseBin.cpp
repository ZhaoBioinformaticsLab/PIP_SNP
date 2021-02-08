/*Author: Wenchao Zhang, Noble Research Institute
  Written on 03/18/2020                          */
#include "CDeepSynthesiseBin.h"
#include <vector>   // For using Vector
#include <assert.h>
#include <math.h>
#include <iostream>  // Cin/ Cout
#include <cfloat>

using namespace std;

CDeepSynthesiseBin::CDeepSynthesiseBin(int Individual_Num, int Corr_Method, int Synthesis_Method, float R2Th)
{
    LDBin_Individual_Num  = Individual_Num;
    LDBin_Corr_Method     = Corr_Method;
    LDBin_Synthesis_Method= Synthesis_Method;
    LDBin_R2Th            = R2Th;
    //ctor
}

/* This function is used to load the read bin line by line to processing bin buffer until it reach to 3*/
void CDeepSynthesiseBin::Load_One_Bin(LDBin_Object Read_Bin)
{
    if(Processing_Bin_Buffer.size()>=3)
    {
       cerr <<"Processing Bin Bufffer is full! There are logical error"<<endl;
       return;
    }
    else
    {
       Processing_Bin_Buffer.push_back(Read_Bin);
    }
    return;
}

/* This function is based on the full loaded (3) processing bin buffer, and used to process the neighbor bins, and
finally update the processing bin buffer with new bin */
void CDeepSynthesiseBin::Process_ReadBin(LDBin_Object Read_Bin)
{
    if(Processing_Bin_Buffer.size()!=3)
    {
        cerr << "The processing bin buffer is not full(3), logical error!"<<endl;
        return;
    }
    else
    {
        float R2_01, R2_02; //Declare Two variables
        LDBin_Object m_Bin0 = Processing_Bin_Buffer.front();
        LDBin_Object m_Bin1 = Processing_Bin_Buffer.at(1);
        LDBin_Object m_Bin2 = Processing_Bin_Buffer.back();
        switch(LDBin_Corr_Method)
        {
           case Pairwise_PCC:  // 0:
               R2_01 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin1.LDBinData, LDBin_Individual_Num, ABS_CORR);
               R2_02 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin2.LDBinData, LDBin_Individual_Num, ABS_CORR);
               break;
           case Pairwise_LD:  //1:
               R2_01 = Calculate_Pairwise_LD(m_Bin0.LDBinData, m_Bin1.LDBinData, LDBin_Individual_Num, ABS_CORR);
               R2_02 = Calculate_Pairwise_LD(m_Bin0.LDBinData, m_Bin2.LDBinData, LDBin_Individual_Num, ABS_CORR);
               break;
           default:
               R2_01 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin1.LDBinData, LDBin_Individual_Num, ABS_CORR);
               R2_02 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin2.LDBinData, LDBin_Individual_Num, ABS_CORR);
               break;
        }

        if(R2_01 >=LDBin_R2Th)
        {
           /* Indicate that the first two bins are correlated to each other*/
           Continuous_Bin_Buffer.push_back(m_Bin0);  // Back up the first bin to continuous buffer
           Processing_Bin_Buffer.erase(Processing_Bin_Buffer.begin()); // Remove the first
        }
        else
        {
            if(R2_02>=LDBin_R2Th)
            {
                /* Indicate that the first and thrid bins are correlated to each other*/
                Uncontinuous_Bin_Buffer.push_back(m_Bin1); // Back up the second bin to continuous buffer
                Processing_Bin_Buffer.erase(Processing_Bin_Buffer.begin()+1); //Remove the second
            }
            else
            {
                /*Indicate thhat the first bin is not correlated with its two neighboring bins*/
                Continuous_Bin_Buffer.push_back(m_Bin0);  // Back up the first bin to continuous buffer
                Processing_Bin_Buffer.erase(Processing_Bin_Buffer.begin()); // Remove the first
                /*Need to process the continuous and uncontinuous buffers through synthesising and combination*/
                Synthesis_Combine_Process();
            }
        }

        assert(Processing_Bin_Buffer.size()==2);
        Processing_Bin_Buffer.push_back(Read_Bin);
    }

    return;
}

/*This function is used to update process the continuous and uncontinuous bin buffers*/
void CDeepSynthesiseBin::Synthesis_Combine_Process()
{
    /*all the bins in continuous buffer will be synthesisied into only one synthesised bin*/
    assert(Continuous_Bin_Buffer.size()>0);
    int Bin_Num =Continuous_Bin_Buffer.size();
    LDBin_Object Continuous_Synthesis_Bin;
    float *pdata_ContinuousBuffer= new float[LDBin_Individual_Num*Bin_Num];
    float *pdata_Synthesised     = new float[LDBin_Individual_Num];
    int i_bin =0;
    int Synthesised_Bin_Start=DBL_MAX;
    int Synthesised_Bin_End  =0;
    int Bin_Pos_Array[Bin_Num];
    for (vector <LDBin_Object>::iterator it_cont_bin =Continuous_Bin_Buffer.begin(); it_cont_bin !=Continuous_Bin_Buffer.end(); it_cont_bin++)
    {
        for (int i_indiviual=0; i_indiviual< LDBin_Individual_Num; i_indiviual++)
        {
             float value = *((*it_cont_bin).LDBinData +i_indiviual);
             *(pdata_ContinuousBuffer+ i_bin*LDBin_Individual_Num+ i_indiviual) =value;
        }

        if((*it_cont_bin).LDBinData!=NULL)
        {
            delete (*it_cont_bin).LDBinData; //we free the linebin here, and it was allocated in the main() function
            (*it_cont_bin).LDBinData =NULL;
        }

        if(Synthesised_Bin_Start > (*it_cont_bin).LDBin_Start)
            Synthesised_Bin_Start = (*it_cont_bin).LDBin_Start;

        if(Synthesised_Bin_End < (*it_cont_bin).LDBin_End)
           Synthesised_Bin_End = (*it_cont_bin).LDBin_End;

        Bin_Pos_Array[i_bin]   = (*it_cont_bin).LDBin_Pos;

        i_bin++;
    }
    Continuous_Bin_Buffer.clear();

    int Pos_Index=0;
    Synthesis_Process(pdata_ContinuousBuffer, pdata_Synthesised, LDBin_Individual_Num, Bin_Num, &Pos_Index, LDBin_Corr_Method, LDBin_Synthesis_Method);

    Continuous_Synthesis_Bin.LDBinData  =pdata_Synthesised;
    Continuous_Synthesis_Bin.LDBin_Start=Synthesised_Bin_Start;
    Continuous_Synthesis_Bin.LDBin_End  =Synthesised_Bin_End;
    Continuous_Synthesis_Bin.LDBin_Pos  =Bin_Pos_Array[Pos_Index];
    DS_Processed_Bin_Buffer.push_back(Continuous_Synthesis_Bin);

    if(pdata_ContinuousBuffer!=NULL)
    {
       delete pdata_ContinuousBuffer;
       pdata_ContinuousBuffer =NULL;
    }

    /*all the bins in the uncontinuous buffer will be syntesisied into multiple bins*/
    if(Uncontinuous_Bin_Buffer.size()>0)
    {
        while(!Uncontinuous_Bin_Buffer.empty())
        {
            LDBin_Object Uncon_Bin1= Uncontinuous_Bin_Buffer.front();
            int UB_Size=Uncontinuous_Bin_Buffer.size();
            if(UB_Size==1) // only 1 bin
            {
                DS_Processed_Bin_Buffer.push_back(Uncon_Bin1);
                Uncontinuous_Bin_Buffer.erase(Uncontinuous_Bin_Buffer.begin()); // vector do not have member function as pop_front() remove the only one bin
            }
            else // 2 more bins
            {
                float R2_Neighbor2;
                float *pdata_bin1;
                float *pdata_bin2;

                Uncontinuous_Bin_Buffer.erase(Uncontinuous_Bin_Buffer.begin()); //pop_front();
                LDBin_Object Uncon_Bin2 = Uncontinuous_Bin_Buffer.front();
                pdata_bin1= Uncon_Bin1.LDBinData;
                pdata_bin2= Uncon_Bin2.LDBinData;

                R2_Neighbor2 = Calculate_r2_one_bin_to_another_bin(pdata_bin1, pdata_bin2, LDBin_Individual_Num, LDBin_Corr_Method);

                if(R2_Neighbor2>LDBin_R2Th) //merge the two neighnor bins into one due to their high correlation
                {
                   int bin1_start, bin2_start, bin1_end, bin2_end, bin1_pos, bin2_pos, synthesis_start, synthesis_end, synthesis_pos;
                   bin1_start= Uncon_Bin1.LDBin_Start;
                   bin2_start= Uncon_Bin2.LDBin_Start;
                   bin1_end  = Uncon_Bin1.LDBin_End;
                   bin2_end  = Uncon_Bin2.LDBin_Start;
                   bin1_pos  = Uncon_Bin1.LDBin_Pos;
                   bin2_pos  = Uncon_Bin2.LDBin_Pos;
                   synthesis_start = (bin1_start< bin2_start)? bin1_start:bin2_start;
                   synthesis_end = (bin1_end < bin2_end)? bin2_end:bin1_end;
                   synthesis_pos = (bin1_pos+ bin2_pos)/2;
                   Uncon_Bin1.LDBin_Start = synthesis_start;
                   Uncon_Bin1.LDBin_End   = synthesis_end;
                   if(LDBin_Synthesis_Method ==Norm_Synthesis)
                   {
                      Uncon_Bin1.LDBin_Pos = synthesis_pos;
                      for(int i_individual =0; i_individual< LDBin_Individual_Num; i_individual++)
                      {
                          float temp1=*(pdata_bin1+i_individual);
                          float temp2=*(pdata_bin2+i_individual);
                          *(pdata_bin1+i_individual)=sqrt(temp1*temp1+temp2*temp2);
                      }
                  }

                  if(pdata_bin2!=NULL)
                  {
                      delete pdata_bin2; //note, we free a linebin here, and it was allocated in the main() function
                      pdata_bin2=NULL;
                  }

                  DS_Processed_Bin_Buffer.push_back(Uncon_Bin1);
                  Uncontinuous_Bin_Buffer.erase(Uncontinuous_Bin_Buffer.begin()); // pop_front() remove the second bin
              }
              else // the two neighbor bins are not correlated
              {
                  DS_Processed_Bin_Buffer.push_back(Uncon_Bin1);
              }
          }
       }
    }

    /*the synthesised bin from continuous buffer and the synthesised bins from uncontinuous buffer will be combined or */
    sort(DS_Processed_Bin_Buffer.begin(), DS_Processed_Bin_Buffer.end(), compareLDBin);

    return;
}

void CDeepSynthesiseBin::Flush_process()
{
   float R2_01, R2_02; //Declare Two variables
   LDBin_Object m_Bin0 = Processing_Bin_Buffer.front();
   LDBin_Object m_Bin1 = Processing_Bin_Buffer.at(1);
   LDBin_Object m_Bin2 = Processing_Bin_Buffer.back();
   switch(LDBin_Corr_Method)
   {
      case Pairwise_PCC:  // 0:
          R2_01 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin1.LDBinData, LDBin_Individual_Num, ABS_CORR);
          R2_02 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin2.LDBinData, LDBin_Individual_Num, ABS_CORR);
          break;
      case Pairwise_LD:  //1:
          R2_01 = Calculate_Pairwise_LD(m_Bin0.LDBinData, m_Bin1.LDBinData, LDBin_Individual_Num, ABS_CORR);
          R2_02 = Calculate_Pairwise_LD(m_Bin0.LDBinData, m_Bin2.LDBinData, LDBin_Individual_Num, ABS_CORR);
          break;
      default:
          R2_01 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin1.LDBinData, LDBin_Individual_Num, ABS_CORR);
          R2_02 = Calculate_Pairwise_PCC(m_Bin0.LDBinData, m_Bin2.LDBinData, LDBin_Individual_Num, ABS_CORR);
          break;
     }

     Continuous_Bin_Buffer.push_back(m_Bin0);
     if(R2_01>=LDBin_R2Th)
     {
        Continuous_Bin_Buffer.push_back(m_Bin1);
     }
     else
     {
        Uncontinuous_Bin_Buffer.push_back(m_Bin1);
     }

     if(R2_02>=LDBin_R2Th)
     {
        Continuous_Bin_Buffer.push_back(m_Bin2);
     }
     else
     {
        Uncontinuous_Bin_Buffer.push_back(m_Bin2);
     }

     Processing_Bin_Buffer.clear();
     /*After Flush the Processing Buffer' contents, we need to synthesis & combine processing to the continuous and uncontinuous buffers*/
     Synthesis_Combine_Process();
     return;

}

CDeepSynthesiseBin::~CDeepSynthesiseBin()
{
    while(!Processing_Bin_Buffer.empty())
    {
       LDBin_Object bin_ob =Processing_Bin_Buffer.back();
       float *pdata= bin_ob.LDBinData ;
       if(pdata!=NULL)
       {
          delete pdata;
          pdata =NULL;
       }
       Processing_Bin_Buffer.pop_back();
    }

    while(!Continuous_Bin_Buffer.empty())
    {
       LDBin_Object bin_ob =Continuous_Bin_Buffer.back();
       float *pdata= bin_ob.LDBinData ;
       if(pdata!=NULL)
       {
          delete pdata;
          pdata =NULL;
       }

       Continuous_Bin_Buffer.pop_back();
    }

    while(!Uncontinuous_Bin_Buffer.empty())
    {
       LDBin_Object bin_ob =Uncontinuous_Bin_Buffer.back();
       float *pdata= bin_ob.LDBinData ;
       if(pdata!=NULL)
       {
          delete pdata;
          pdata =NULL;
       }
       Uncontinuous_Bin_Buffer.pop_back();
    }

    while(!DS_Processed_Bin_Buffer.empty())
    {
       LDBin_Object bin_ob =DS_Processed_Bin_Buffer.back();
       float *pdata= bin_ob.LDBinData ;
       if(pdata!=NULL)
       {
          delete pdata;
          pdata =NULL;
       }

       DS_Processed_Bin_Buffer.pop_back();
    }
    //dtor
}
