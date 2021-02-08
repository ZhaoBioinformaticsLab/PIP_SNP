/*Written by Wenchao Zhang, Noble Research Institute
  01/08/2019 */
#include <iostream>  // Cin/ Cout
#include <fstream>  // For Read/Write File
#include <string>
#include <sstream>  // using sstream
#include <vector>   // For using Vector
#include <list>     // For Using List
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "Common.h"
#include "LDBin_DataObject.h"
#include "CImputingMissing.h"

using namespace std;

int parse_cmd_line(int argc, char* argv[], string &G_File_Name, string &LD_BinMap_FileName, int &Individual_Num, int &Method_Corr, int &KNN, int &Method_Syn, string &Out_File_Name, string &Out_LD_Map_FileName)
{
	if((argc!=2)&&(argc <17))
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

           if((arg=="-m")||(arg=="-M"))
           {
              if(i+1<argc)
		      {
                   Out_LD_Map_FileName = argv[i+1];
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the output file name for LD block Map after synthesizing!" <<endl;
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
				 cout << "-c or -C: the  Correlation Method for a Pairwise Genotypic Markers, 0: Pearson_Correlation_R2; 1: LD_D_R2; Default(>=2): Pearson_Correlation_R2" <<endl;
                 cout << "-k or -K: the  K Nearest Neighbor individuals in The KNN method " <<endl;
			     cout << "-s or -S: the method how to generate the syntesized(binned) genotypic marker, 0: not to synthesize; 1: norm integration of the multiple markers' genotype values; 2: select the optimal one; default(>2): norm integration   " <<endl;
			     cout << "-o or -O: the  full name of Imputing results" <<endl;
			     cout << "-m or -M: the  full name of LD Mapping results after synthesing" <<endl;
			     return 1; // Parse cmd_line fail
		    }
	    }
	}
    return 0;
}


int main(int argc, char* argv[])
{
    cout << "Hello, Dear User! Welcome to use SNP Processing Pipeline for LD Bin Filling, Missing Genotype Imputing, and/or Marker Synthesizing!" << endl;

    string G_File_Name, LD_BinMap_FileName, Out_File_Name, Out_LD_Map_FileName;
    int Individual_Num;
    int Method_Corr;
    int KNN;
    int Method_Syn;

     // Parse the command line for the inputting
    if(1==parse_cmd_line(argc, argv, G_File_Name, LD_BinMap_FileName, Individual_Num, Method_Corr, KNN, Method_Syn, Out_File_Name, Out_LD_Map_FileName))
        return 1;

    // Begin to read
    //Open and Read the Geneotype Matrix Marker_Block one by one.
	ifstream G_File(G_File_Name.c_str(), ios::in);
	ifstream Binmap_File(LD_BinMap_FileName.c_str(), ios::in);
	ofstream O_File(Out_File_Name.c_str(), ios::out);
	ofstream O_Binmap_File(Out_LD_Map_FileName.c_str(), ios::out);

	char out_delim =',';

    if(!G_File.is_open())
	{
		cerr<< G_File_Name <<" Can't be accessed!"<<endl;
		return 1;
	}

	if(!Binmap_File.is_open())
	{
		cerr<< LD_BinMap_FileName <<" Can't be accessed!"<<endl;
		return 1;
	}

	if(G_File && Binmap_File)
    {
       string sLine;
       int Line_Count=0;
       char *marker_line;

       list <LDBin_Map_Feature> List_LDBin_MapFeatures;

       while(getline(Binmap_File, sLine)) // Read the LD Bin Map File
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
			if(Col_Num < 2)
			{
				s_v.clear();
				while(getline(ss, item, delim2))// try delim2;
			    {
				   s_v.push_back(item);
			    }
			    Col_Num=s_v.size();
			}

			if(Col_Num <2)
			{
				cerr<< Binmap_File <<"File Format is not right(should contain at least two int number), Error at Line=" << Line_Count<<endl;
		        return 1;
			}

            int line_start= atoi(s_v.at(0).c_str());
            int line_end  = atoi(s_v.at(1).c_str());
            if((line_start==0)&&(line_end==0)) continue ; // If meet the head information, by pass

            LDBin_Map_Feature bin_map;
            bin_map.LDBin_Start = line_start;
            bin_map.LDBin_End   = line_end;
            List_LDBin_MapFeatures.push_back(bin_map);
            Line_Count ++;
		  }
	   }

       LDBin_Map_Feature bin_map =List_LDBin_MapFeatures.front();
       Line_Count                =0;
       CLDBin_DataObject LDBin_Data(Individual_Num,Method_Corr,bin_map.LDBin_Start,bin_map.LDBin_End);
       if(bin_map.LDBin_Start!= 0)
       {
           cerr << "The LD Bin Map File is not correctly recorded, the first LDBin_Start should start from 0"<< endl;
           return 1;
       }

       O_Binmap_File << "SynthesisBin_Start" << out_delim << "SynthesisBin_End" << out_delim << "SynthesisBin_Pos" <<endl;

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

			// After Read one marker line

			if(Line_Count > bin_map.LDBin_End)
            {
                cerr << "LD Bin Map File and Genotype Marker File are not matched! LDBin_End= "<<  bin_map.LDBin_End << "But, The next recorded Line Count=" << Line_Count<< endl;
                return 1;
            }
            else
            {
                // LD Bin buffer is not reach to the LD Mapping indicator, Add the line buffer into the LDBin Buffer
                LDBin_Data.Add_one_markerline(Line_Count,marker_line);

                if(Line_Count == bin_map.LDBin_End)	// LDBin Buffer
                {
                     // One LD Bin is ready, we need to process it
                     /*The codes for imputing the missing genotype values that located in the detected LD Bin Regions*/
                     // Initialized the Cimputing Object with detected LDBin Data and KNN (10);
                     CImputingMissing LDBinData_Impute(LDBin_Data, KNN);
                     LDBinData_Impute.Impute_all();

                     /*Finished for imputing the missing genotypes located in the detected LD Bin region */

                     if(Method_Syn>0)
                     {
                         /*And Need to synthesis/bin the SNP markers of a LDBin into one representative marker*/
                         LDBinData_Impute.Synthesize(Method_Syn);

                         LDBinData_Impute.LDBin_Pos+= LDBinData_Impute.LDBin_Start;

                         O_Binmap_File << LDBinData_Impute.LDBin_Start << out_delim << LDBinData_Impute.LDBin_End << out_delim << LDBinData_Impute.LDBin_Pos <<endl;

                         int i_individual;
                         for(i_individual=0; i_individual<LDBinData_Impute.LDBin_Individual_Num-1; i_individual++)
                         {
                            O_File<< LDBinData_Impute.Lpldbin_Synthesized[i_individual]<< out_delim;
                         }
                         O_File<< LDBinData_Impute.Lpldbin_Synthesized[i_individual]<< endl;
                      }
                      else
                      {

                         LDBinData_Impute.LDBin_Pos+= LDBinData_Impute.LDBin_Start;
                         O_Binmap_File << LDBinData_Impute.LDBin_Start << out_delim << LDBinData_Impute.LDBin_End << out_delim << LDBinData_Impute.LDBin_Pos <<endl;

                         /*Need to output the imputation results to a new out file*/
                         for(int i_marker=0; i_marker<LDBinData_Impute.LDBin_Marker_Num; i_marker++)
                         {
                             int i_individual;
                             for(i_individual=0; i_individual<LDBinData_Impute.LDBin_Individual_Num-1; i_individual++)
                             {
                                 O_File<< int(LDBinData_Impute.Lpldbin_ImputedData[i_marker*LDBinData_Impute.LDBin_Individual_Num+i_individual])<<out_delim;
                             }
                             O_File<< int(LDBinData_Impute.Lpldbin_ImputedData[i_marker*LDBinData_Impute.LDBin_Individual_Num+i_individual])<< endl;
                          }
                       }

                       // After processing, we need to prepare another LD Bin
                       List_LDBin_MapFeatures.pop_front();
                       if(!List_LDBin_MapFeatures.empty())
                       {
                           bin_map =List_LDBin_MapFeatures.front();
			               LDBin_Data.Reset_LDBin_DataObject(bin_map.LDBin_Start, bin_map.LDBin_End);
			               if(bin_map.LDBin_Start!= (Line_Count+1))
			               {
			                  cerr << "LD Bin Map File and Genotype Marker File are not matched! LDBin_Start= "<<  bin_map.LDBin_Start << "But, The next recorded Line Count=" << (Line_Count+1)<< endl;
			                  return 1;
			               }
                       }
                }
            }

          Line_Count++;
          } // Finish one line;

	   }

	   LDBin_Data.Clean_LDBin_marker_lines();

    }

    G_File.close();
    Binmap_File.close();
    O_File.close();
    O_Binmap_File.close();

    cout << "SNP Data Processing Pipeline Finished!  " <<endl;

    return 0;
}
