/*Author: Wenchao Zhang, Noble Research Institute
  Written on 03/18/2020                          */
#ifndef CDEEPSYNTHESISEBIN_H
#define CDEEPSYNTHESISEBIN_H
#include <vector>   // For using Vector
#include "Common.h"

/*This Class is defined for a deep synthesis processing for the initial constructed Bins, which will be finally directly output*/
class CDeepSynthesiseBin
{
    public:
    int   LDBin_Individual_Num;
    int   LDBin_Corr_Method;
    int   LDBin_Synthesis_Method;
    float LDBin_R2Th;

    vector <LDBin_Object> Processing_Bin_Buffer;    // Used to dynamically store the processing Bins, which are connected witht he readin Bin lines
    vector <LDBin_Object> Continuous_Bin_Buffer;    // Used to dynamically store the continuous (Clustered Block) Bins
    vector <LDBin_Object> Uncontinuous_Bin_Buffer;  // Used to dynamically store the uncontinuous Jumping Bins
    vector <LDBin_Object> DS_Processed_Bin_Buffer;  // Used to dynamically store the post-processed DS Bins, which can be output to a file

    void Load_One_Bin(LDBin_Object Read_Bin);       // use to add one line (together with LD Bin information) into the Processing_Bin Buffer until it's full.
    void Process_ReadBin(LDBin_Object Read_Bin);   // Process how to deal with new read line bin and the process_bin_buffer
    void Synthesis_Combine_Process();                 // Process the continuous and uncontinuous buffers
    void Flush_process();                             // Used to flush process when there is not new readline.
        CDeepSynthesiseBin(int Individual_Num, int Corr_Method, int Synthesis_Method, float R2Th);
        virtual ~CDeepSynthesiseBin();
    protected:
    private:
};

#endif // CDEEPSYNTHESISEBIN_H
