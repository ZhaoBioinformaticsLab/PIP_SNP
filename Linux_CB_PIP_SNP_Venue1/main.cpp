/*Written by Wenchao Zhang, Noble Research Institute
  01/08/2020 */
#include <iostream>  // Cin/ Cout
#include <fstream>  // For Read/Write File
#include <string>
#include <sstream>  // using sstream
#include <vector>   // For using Vector
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstddef>  // std::size_t
#include "Common.h"
#include "LDBin_DataObject.h"
#include "CImputingMissing.h"

using namespace std;
int parse_cmd_line(int argc, char* argv[], string & G_File_Name, string & LD_BinMap_FileName, int & Individual_Num, float & LD_R2_Th, int &Method_Corr, int &Method_Detection, int &KNN, int &Method_Syn, string &Out_File_Name);
int LDBin_Detection(int Method_Detection, float R2_Right_Breakthrough, float R2_Left_Breakthrough, float LD_R2_Th);
void Fill_LDBin_Feature(CLDBin_DataObject &LDBin_DataObj, float R2_Left_Breakthrough, float R2_Right_Breakthrough,  LDBin_Map_Feature &LDBin_Feature);

int parse_cmd_line(int argc, char* argv[], string &G_File_Name, string &LD_BinMap_FileName, int &Individual_Num, float &LD_R2_Th, int &Method_Corr, int &Method_Detection, int &KNN, int &Method_Syn, string &Out_File_Name)
{
	if((argc!=2)&&(argc <19))
	{
		cerr <<"Please specific the correct parameters, or use parameter -u for user manuals!" <<endl;
		return 1;
	}
	else
	{
	   for (int i=1; i<argc;i++)
       {
	       string arg=argv[i];
	       if((arg=="-g")||(arg=="-G"))
	       {
             if(i+1<argc)
		     {
		       G_File_Name =argv[i+1];
		     }
		     else
		     {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the full file name of genotype matrix !" <<endl;
			   return 1; // Parse cmd_line fail
		      }
		    }

		   if((arg=="-l")||(arg=="-L"))
	       {
             if(i+1<argc)
		     {
		       LD_BinMap_FileName =argv[i+1];
		     }
		     else
		     {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the full file name of the output LD Bin Mapping Result!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
		    }

            if((arg=="-i")||(arg=="-I"))
	        {
              if(i+1<argc)
		      {
		           Individual_Num =atoi(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the Individual Number!" <<endl;
			   return 1; // Parse cmd_line failofstream
		      }
	        }

	        if((arg=="-r")||(arg=="-R"))
	        {
              if(i+1<argc)
		      {
		           LD_R2_Th =atof(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the Pairwise LD R2 Threshold!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }


	        if((arg=="-c")||(arg=="-C"))
	        {
              if(i+1<argc)
		      {
                   Method_Corr =atoi(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the  Correlation Method for a Pairwise Genotypic Markers!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }

            if((arg=="-d")||(arg=="-D"))
	        {
              if(i+1<argc)
		      {
                   Method_Detection =atoi(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the  Detecting Method for a Marker group(LD Bin) !" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }

            if((arg=="-k")||(arg=="-K"))
	        {
              if(i+1<argc)
		      {
                   KNN =atoi(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the KNN (K nearest neighbors) individuals for imputing a missing genotype form a group(LD Bin) !" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }

	        if((arg=="-s")||(arg=="-S"))
	        {
              if(i+1<argc)
		      {
                   Method_Syn =atoi(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the integer value for how to generate or not generate the synthesizied (binned) marker for the detected group(LD Bin) !" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }

           if((arg=="-o")||(arg=="-O"))
           {
              if(i+1<argc)
		      {
                   Out_File_Name = argv[i+1];
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the output file name for storing the result after imputing the missing genotypes or/and synthesizing multiple Genotype SNP markers of the LD blocks!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
           }

		   if((arg=="-u")||(arg=="-U"))
	       {
			     cout<<"Welcome to use this program to do LD Marker Bin Detecting" <<endl;
				 cout << "The usuage of input parameter arguments are listed as followings:" <<endl;
				 cout << "-u or -U: Output this help usuage message" <<endl;
				 cout << "-g or -G: The full name of Genotype file"<<endl;
				 cout << "-l or -L: The full name of LD Bin Mapping Result file"<<endl;
				 cout << "-i or -I: The Individual number" <<endl;
				 cout << "-r or -R: The Threshold for the pairwise LD R2" <<endl;
				 cout << "-c or -C: the  Correlation Method for a Pairwise Genotypic Markers,0:Pearson_Correlation_R2; 1: LD_D_R2; Default(>2): Pearson_Correlation_R2 " <<endl;
				 cout << "-d or -D: the  Detection Method for a Marker Group(LD Block), 0: Right_Breakthrough; 1: Left_Breakthrough; 2: Left and Right Breakthrough; 3: Left or Right Breakthrough; Default(>3) : Right_Breakthrough" <<endl;
			     cout << "-k or -K: the  K Nearest Neighbor individuals in The KNN method " <<endl;
			     cout << "-s or -S: the method how to generate the syntesized(binned) genotypic marker, 0: not to synthesize; 1: norm integration of the multiple markers' genotype values; 2: select the optimal one; Default(>2) :norm integration of the multiple markers' genotype values  " <<endl;
			     cout << "-o or -O: the  full name of Imputing results" <<endl;
			     return 1; // Parse cmd_line fail
		    }
	    }
	}
    return 0;
}

int LDBin_Detection(int Method_Detection, float R2_Right_Breakthrough, float R2_Left_Breakthrough, float LD_R2_Th )
{
    int ret =0;
    switch (Method_Detection)
    {
        case Right_Breakthorugh:
            if(fabs(R2_Right_Breakthrough)< LD_R2_Th)  ret =1;
            break;
        case Left_Breakthorugh:
            if(fabs(R2_Left_Breakthrough)< LD_R2_Th)  ret =1;
            break;
        case RandL_Breakthrough:
            if((fabs(R2_Right_Breakthrough)< LD_R2_Th)&&(fabs(R2_Left_Breakthrough)< LD_R2_Th))  ret =1;
            break;
        case RorL_Breakthrough:
            if((fabs(R2_Right_Breakthrough)< LD_R2_Th)||(fabs(R2_Left_Breakthrough)< LD_R2_Th))  ret =1;
            break;
        default:
            if(fabs(R2_Right_Breakthrough)< LD_R2_Th)  ret =1;
            break;

    }
    return ret;
}

/*Once we detect a LDBin Object, we need to acquire and fill the LDBin_Map_Feature */
void Fill_LDBin_Feature(CLDBin_DataObject &LDBin_DataObj, float R2_Left_Breakthrough, float R2_Right_Breakthrough,  LDBin_Map_Feature &LDBin_Feature)
{

    int R2_PairNum = LDBin_DataObj.Get_markerlines_number()-1;

    LDBin_Feature.LDBin_End               =LDBin_DataObj.LDBin_Current;
    LDBin_Feature.LDBin_Start             =LDBin_DataObj.LDBin_Start;
    LDBin_Feature.LDBin_Pos               =(LDBin_Feature.LDBin_End -LDBin_Feature.LDBin_Start)/2;

    LDBin_Feature.R2_LD_Average_Forward   =0.0;
    LDBin_Feature.R2_LD_Average_Neighbor2 =0.0;
    if(R2_PairNum>0)
    {
       LDBin_Feature.R2_LD_Average_Forward   =LDBin_DataObj.Sum_R2_Left_To_Current/R2_PairNum;
       LDBin_Feature.R2_LD_Average_Neighbor2 =LDBin_DataObj.Sum_R2_Neighbor2/R2_PairNum;
    }

    LDBin_Feature.R2_LD_Left2                =LDBin_DataObj.R2_LD_Left2;
    LDBin_Feature.R2_LD_Right2               =LDBin_DataObj.R2_LD_Right2;
    LDBin_Feature.R2_LD_Left_Right           =LDBin_DataObj.R2_LD_Left_Right;

    LDBin_Feature.R2_Left_Breakthrough       =R2_Left_Breakthrough;
    LDBin_Feature.R2_Right_Breakthrough      =R2_Right_Breakthrough;

    return;

}

int main(int argc, char* argv[])
{
    cout << "Hello, Dear User! SNP Processing Pipeline for LD Bin Detecting, Missing Genotype Imputing, and or Synthesizing is Starting!" << endl;

    string G_File_Name, LD_BinMap_FileName, Out_File_Name;
    int Individual_Num;
    int Method_Corr;
    int Method_Detection;
    float LD_R2_Th;
    int KNN;
    int Method_Syn;

     // Parse the command line for the inputting
    if(1==parse_cmd_line(argc, argv, G_File_Name, LD_BinMap_FileName, Individual_Num, LD_R2_Th, Method_Corr, Method_Detection, KNN, Method_Syn, Out_File_Name))
        return 1;

    // Begin to read
    //Open and Read the Geneotype Matrix Marker_Block one by one.
	ifstream G_File(G_File_Name.c_str(), ios::in);
	ofstream Binmap_File(LD_BinMap_FileName.c_str(), ios::out);
	ofstream O_File(Out_File_Name.c_str(), ios::out);


	char out_delim =',';

    if(!G_File.is_open())
	{
		cerr<< G_File_Name <<" Can't be accessed!"<<endl;
		return 1;
	}

	if(G_File)
    {
       string sLine;
       int Line_Count=0;
       char *marker_line;

       vector <LDBin_Map_Feature> Vector_LDBin_MapFeatures;
       CLDBin_DataObject LDBin_Data(Individual_Num,Method_Corr, Line_Count);

	   while(getline(G_File, sLine))
	   {
		  if(sLine.empty()) ; // Ignore empty lines
		  else
		  {
            stringstream ss(sLine);
			vector <string> s_v;
			string item;
			char delim1 =',';
			char delim2 ='\t';

			while(getline(ss, item, delim1))
			{
				s_v.push_back(item);
			}

			int Col_Num=s_v.size();
			if (Col_Num < Individual_Num)
			{
				s_v.clear();
				while(getline(ss, item, delim2))// try delim2;
			    {
				   s_v.push_back(item);
			    }
			    Col_Num=s_v.size();
			}

			if(Col_Num != Individual_Num)
			{
				cerr<< G_File_Name <<"File Format is not right, Error at Line=" << Line_Count<<endl;
		        return 1;
			}

			marker_line = new char [Individual_Num];
			for (int i=0; i<Col_Num;i++)
			{
			   char value= atoi(s_v.at(i).c_str());
			     marker_line[i] = value;
			}

			// Read one marker line and need to make a decision how to process it.
			float R2_Left_Breakthrough =0.0;
            float R2_Right_Breakthrough=0.0;

			if (LDBin_Data.Get_markerlines_number()==0)
            {
                // Indicate the LDBin Object is empty, just add it
                LDBin_Data.Add_one_markerline(Line_Count,marker_line, R2_Left_Breakthrough,R2_Right_Breakthrough);
            }
            else
            {
                // Indicate the LDBin Object have some marker already, we need to make some decisions
                R2_Left_Breakthrough = LDBin_Data.Calculate_R2_Left_one_markerline(marker_line);
                R2_Right_Breakthrough= LDBin_Data.Calculate_R2_Right_one_markerline(marker_line);
                if(LDBin_Detection(Method_Detection, R2_Right_Breakthrough, R2_Left_Breakthrough, LD_R2_Th )==1)
                {
                    // Indicate a new LDBin has been detected.
                    // We need to get the LDBin Map Feature Information from the LDBin Object
                    LDBin_Map_Feature LDBin_Feature;

                    Fill_LDBin_Feature(LDBin_Data, R2_Left_Breakthrough, R2_Right_Breakthrough, LDBin_Feature);
                    /*The codes for imputing the missing genotype values that located in the detected LD Bin Regions*/
                    // Initialized the Cimputing Object with detected LDBin Data and KNN (10);
                    CImputingMissing LDBinData_Impute(LDBin_Data, KNN);
                    LDBinData_Impute.Impute_all();

                    /*Finished for imputing the missing genotypes located in the detected LD Bin region */

                    if(Method_Syn>0)
                    {
                      /*And Need to synthesis/bin the SNP markers of a LDBin into one representative marker*/
                      LDBinData_Impute.Synthesize(Method_Syn);
                      int i_individual;
                      for(i_individual=0; i_individual<LDBinData_Impute.LDBin_Individual_Num-1; i_individual++)
                      {
                         O_File<< LDBinData_Impute.Lpldbin_Synthesized[i_individual]<< out_delim;
                      }
                      O_File<< LDBinData_Impute.Lpldbin_Synthesized[i_individual]<< endl;

                      LDBin_Feature.LDBin_Pos = LDBinData_Impute.LDBin_Pos; // Update the Position of Bin with the represenatative pos
                    }
                    else
                    {
                        /*Need to output the imputation results to a new out file*/
                        for (int i_marker=0; i_marker<LDBinData_Impute.LDBin_Marker_Num; i_marker++)
                        {
                           int i_individual;
                           for(i_individual=0; i_individual<LDBinData_Impute.LDBin_Individual_Num-1; i_individual++)
                           {
                              O_File<< int(LDBinData_Impute.Lpldbin_ImputedData[i_marker*LDBinData_Impute.LDBin_Individual_Num+i_individual])<<out_delim;
                           }
                           O_File<< int(LDBinData_Impute.Lpldbin_ImputedData[i_marker*LDBinData_Impute.LDBin_Individual_Num+i_individual])<< endl;
                        }
                    }

                    LDBin_Feature.LDBin_Pos += LDBin_Feature.LDBin_Start;
                    Vector_LDBin_MapFeatures.push_back(LDBin_Feature);
                    // Need to reset the LDBin Object for next LDBin
                    LDBin_Data.Reset_LDBin_DataObject(Line_Count, marker_line);
                }
                else
                {
                    //Indicate there is not a new LDBin that has been detected, still use the old LDBin
                    LDBin_Data.Add_one_markerline(Line_Count,marker_line, R2_Left_Breakthrough,R2_Right_Breakthrough);

                }

            }

          } // Finish one line;
		  Line_Count++;
	   }

	   /* remember to process the remain marker lines that stored in the LDBin Object */
	   LDBin_Map_Feature LDBin_Feature;
       Fill_LDBin_Feature(LDBin_Data, 0.0, 0.0, LDBin_Feature);

       /*The codes for imputing the missing genotype values that located in the detected LD Bin Regions*/
       // int KNN =10;
        CImputingMissing LDBinData_Impute(LDBin_Data, KNN);
        LDBinData_Impute.Impute_all();
        /*Finished for imputing the missing genotypes located in the detected LD Bin region */
        if(Method_Syn>0)
        {
            /*And Need to synthesis/bin the SNP markers of a LDBin into one representative marker*/
            LDBinData_Impute.Synthesize(Method_Syn);
            int i_individual;
            for(i_individual=0; i_individual<LDBinData_Impute.LDBin_Individual_Num-1; i_individual++)
            {
                O_File<< LDBinData_Impute.Lpldbin_Synthesized[i_individual]<< out_delim;
            }
            O_File<< LDBinData_Impute.Lpldbin_Synthesized[i_individual]<< endl;

            LDBin_Feature.LDBin_Pos = LDBinData_Impute.LDBin_Pos; // Update the Position of Bin with the represenatative pos

        }
        else
        {
            /*Finished for imputing the missing genotypes located in the detected LD Bin region */
            /*And Need to output the imputation results to a new out file*/
            for (int i_marker=0; i_marker<LDBinData_Impute.LDBin_Marker_Num; i_marker++)
            {
               int i_individual;
               for(i_individual=0; i_individual<LDBinData_Impute.LDBin_Individual_Num-1; i_individual++)
               {
                 O_File<< int(LDBinData_Impute.Lpldbin_ImputedData[i_marker*LDBinData_Impute.LDBin_Individual_Num+i_individual])<<out_delim;
               }
               O_File<< int(LDBinData_Impute.Lpldbin_ImputedData[i_marker*LDBinData_Impute.LDBin_Individual_Num+i_individual])<< endl;
            }
        }

        LDBin_Feature.LDBin_Pos += LDBin_Feature.LDBin_Start;
        Vector_LDBin_MapFeatures.push_back(LDBin_Feature);// Store the current detected LD_Bin Information.

	    LDBin_Data.Clean_LDBin_marker_lines();

       /* output the detected LD Bin Map index  */

       Binmap_File<< "Bin_Start" << out_delim;
       Binmap_File<< "Bin_End"<< out_delim;
       Binmap_File<< "Bin_Pos"<< out_delim;
       Binmap_File<< "R2_LD_Left2" << out_delim;
       Binmap_File<< "R2_LD_Right2" << out_delim;
       Binmap_File<< "R2_LD_Left_Right" << out_delim;
       Binmap_File<< "R2_LD_Average_Forward"<<out_delim;
       Binmap_File<< "R2_LD_Average_Neighbor2"<<out_delim;
       Binmap_File<< "R2_Left_Breakthrough"<<out_delim;
       Binmap_File<< "R2_Right_Breakthrough"<<  endl;

       for (vector <LDBin_Map_Feature>::iterator it_LD= Vector_LDBin_MapFeatures.begin(); it_LD!= Vector_LDBin_MapFeatures.end();  it_LD++ )
       {

           Binmap_File<< (*it_LD).LDBin_Start << out_delim;
           Binmap_File<< (*it_LD).LDBin_End<< out_delim;
           Binmap_File<< (*it_LD).LDBin_Pos<< out_delim;
           Binmap_File<< (*it_LD).R2_LD_Left2<< out_delim;      // Recorded the R2 value between the most left two SNP markers
           Binmap_File<< (*it_LD).R2_LD_Right2<< out_delim;     // Recorded the R2 value between the most right two SNP markers
           Binmap_File<< (*it_LD).R2_LD_Left_Right<< out_delim; // Recorded the R2 value between the most left and right SNP markers.
           Binmap_File<< (*it_LD).R2_LD_Average_Forward<<out_delim;   // Recorded the average of R2 values of the start SNP Marker with all other SNP marker at the detected LD Bin
           Binmap_File<< (*it_LD).R2_LD_Average_Neighbor2<<out_delim; // Recorded the average of R2 values of all the two neighbor SNP marker at the detected LD Bin
           Binmap_File<< (*it_LD).R2_Left_Breakthrough<<out_delim;    // Recorded the R2 Value of between the breakthrough and the left point
           Binmap_File<< (*it_LD).R2_Right_Breakthrough<<  endl;      // Recorded the R2 Value of between the breakthrough and the right point
       }

       Vector_LDBin_MapFeatures.clear();
    }

    G_File.close();
    Binmap_File.close();
    O_File.close();

    cout << "A Series of Tentative SNP Data Processing has been Finished!  " <<endl;

    return 0;
}
