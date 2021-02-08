/*Written by Wenchao Zhang, Noble Research Institute
  02/21/2020 */
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
#include <assert.h>
#include "Common.h"
#include "CDeepSynthesiseBin.h"

using namespace std;
int parse_cmd_line(int argc, char* argv[], string & G_File_Name, string & LD_BinMap_FileName, int & Individual_Num, int & Corr_Method, float &R2_Th, int &Syntheis_Method, string &Out_G_File_Name, string &Out_LDBin_File_Name);
int parse_cmd_line(int argc, char* argv[], string & G_File_Name, string & LD_BinMap_FileName, int & Individual_Num, int & Corr_Method, float &R2_Th, int &Syntheis_Method, string &Out_G_File_Name, string &Out_LDBin_File_Name)
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
			   cerr<<"Parse cmd_line fail, Need to clearly specific the full file name of the input LD Bin Mapping!" <<endl;
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
                   Corr_Method =atoi(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the Correlation Method!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }

	        if((arg=="-r")||(arg=="-R"))
	        {
              if(i+1<argc)
		      {
                   R2_Th =atof(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the Pairwise LD R2 Threshold!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }

	        if((arg=="-s")||(arg=="-S"))
	        {
              if(i+1<argc)
		      {
                   Syntheis_Method =atoi(argv[i+1]);
		      }
		      else
		      {
			   cerr<<"Parse cmd_line fail, Need to clearly specific the Deep Synthesis Method!" <<endl;
			   return 1; // Parse cmd_line fail
		      }
	        }

           if((arg=="-o")||(arg=="-O"))
           {
              if(i+1<argc)
		      {
                   Out_G_File_Name = argv[i+1];
		      }
		      else
		      {
			    cerr<<"Parse cmd_line fail, Need to clearly specific the output file name for storing the final processed SNP result!" <<endl;
			    return 1; // Parse cmd_line fail
		      }
           }

           if((arg=="-m")||(arg=="-M"))
           {
              if(i+1<argc)
		      {
                   Out_LDBin_File_Name = argv[i+1];
		      }
		      else
		      {
			    cerr<<"Parse cmd_line fail, Need to clearly specific the output file name for final LD Bin Mapping!" <<endl;
			    return 1; // Parse cmd_line fail
		      }
           }

           if((arg=="-u")||(arg=="-U"))
	       {
			     cout<<"Welcome to use this program to do Deep Synthesising" <<endl;
				 cout << "The usuage of input parameter arguments are listed as followings:" <<endl;
				 cout << "-u or -U: Output this help usuage message" <<endl;
				 cout << "-g or -G: The full name of inputing Synthesised SNP file"<<endl;
				 cout << "-l or -L: The full name of inputting  LD Bin file Mapping Synthesising SNPs  "<<endl;
				 cout << "-i or -I: The Individual number" <<endl;
				 cout << "-c or -C: The Method for the pairwise R2 Correlation, 0:Pearson_Correlation_R2; 1: LD_D_R2; Default(>=2): Pearson_Correlation_R2" <<endl;
				 cout << "-r or -R: The Threshold for the pairwise R2 Threshold, a float value at (0~1.0) " <<endl;
				 cout << "-s or -S: The Method for the multiple SNP bin Synthesising. 1: norm integration of the multiple markers' genotype values; 2: select the optimal one; Default(0, or >2) :norm integration of the multiple markers' genotype values" <<endl;
			     cout << "-o or -O: The full name of Deep Synthesised SNP results" <<endl;
			     cout << "-m or -M: The full LD Bin mapping file name after deep synthesising" <<endl;
			     return 1; // Parse cmd_line fail
		    }
	     }
	}
	return 0;
}

int main(int argc, char* argv[])
{
    cout << "Hello, Dear User! To furrther  reduce the SNP size, we develop such specific program to deep synthesis SNP Bins by considering the correlationship of a SNP and its neighbor jumped SNP Bins!" << endl;

    string G_File_Name, LD_BinMap_FileName, Out_G_File_Name, Out_LDBin_File_Name;
    int Individual_Num, Method_Corr, Method_Synthesis;
    float R2_Th;

    if(1==parse_cmd_line( argc, argv, G_File_Name, LD_BinMap_FileName, Individual_Num, Method_Corr, R2_Th, Method_Synthesis, Out_G_File_Name, Out_LDBin_File_Name))
      return 0;

    ifstream In_G_File(G_File_Name.c_str(), ios::in);
	ifstream In_BinMap_File(LD_BinMap_FileName.c_str(), ios::in);
	ofstream O_G_File(Out_G_File_Name.c_str(), ios::out);
	ofstream O_BinMap_File(Out_LDBin_File_Name.c_str(), ios::out);

    char in_delim =',';
    char out_delim =',';

    if(!In_G_File.is_open())
	{
		cerr<< G_File_Name <<" Can't be accessed!"<<endl;
		return 1;
	}

	if(!In_BinMap_File.is_open())
	{
		cerr<< LD_BinMap_FileName <<" Can't be accessed!"<<endl;
		return 1;
	}

    int Bin_Count_Map =0;
    vector <LDBin_Map_Feature> Input_LDBin_Map_Vector;
    if(In_BinMap_File)
    {
       /*Try to load the Bin Map Information*/
       string sLine;
       getline(In_BinMap_File, sLine); // Read the head information line
       while(getline(In_BinMap_File, sLine))
	   {
	      if(sLine.empty()) ; // Ignore empty lines
		  else
		  {
		     stringstream ss(sLine);
             /* Some code to extract the Bin Map Information*/
             string item;
             LDBin_Map_Feature bin_feature;
             getline(ss, item, in_delim); // Read the first column as the Bin_Start
             bin_feature.LDBin_Start = atoi(item.c_str());
             getline(ss, item, in_delim); // Read the second column as the Bin_End
             bin_feature.LDBin_End   = atoi(item.c_str());
             getline(ss, item, in_delim); // Read the third column as the Bin_Pos
             bin_feature.LDBin_Pos   = atoi(item.c_str());
             bin_feature.R2_Breakpoint=0.0;
             Input_LDBin_Map_Vector.push_back(bin_feature);
             Bin_Count_Map++;
		  }
	   }
    }

    if(In_G_File)
    {
       O_BinMap_File << "DS_Bin_Start" << out_delim << "DS_Bin_End" << out_delim<< "DS_Bin_Pos" <<endl;
       string sLine;
       CDeepSynthesiseBin DeepSynthesis(Individual_Num, Method_Corr, Method_Synthesis, R2_Th);

       int Line_Count=0;
       int DS_Line_Count =0;
       while(getline(In_G_File, sLine))
	   {
          if(sLine.empty()) ; // Ignore empty lines
		  else
		  {
             float *read_line;
             LDBin_Object read_ldbin;
             stringstream ss(sLine);
			 vector <string> s_v;
			 string item;
			 while(getline(ss, item, in_delim))
             {
                s_v.push_back(item);
             }

             if(s_v.size()!= Individual_Num)
             {
                cerr << G_File_Name <<"File Format is not right, Error at Line=" << Line_Count<<endl;
		        return 1;
             }

             read_line  = new float[Individual_Num];
             for (int i=0; i<Individual_Num;i++)
             {
                 float value= atof(s_v.at(i).c_str());
                 read_line[i] = value;
             }

             LDBin_Map_Feature temp_bin_feature= Input_LDBin_Map_Vector.at(Line_Count);
             read_ldbin.LDBin_Start            = temp_bin_feature.LDBin_Start;
             read_ldbin.LDBin_End              = temp_bin_feature.LDBin_End;
             read_ldbin.LDBin_Pos              = temp_bin_feature.LDBin_Pos;
             read_ldbin.LDBinData              = read_line;

             // Load one linebin, then, we need to use the DeepSynthesis Object to implement the deep synthesis processing
             if(DeepSynthesis.Processing_Bin_Buffer.size()<3)
             {
                DeepSynthesis.Load_One_Bin(read_ldbin);
             }
             else
             {
                DeepSynthesis.Process_ReadBin(read_ldbin);
             }


             if(DeepSynthesis.DS_Processed_Bin_Buffer.size()>0)
             {
                 int DSBin_Start, DSBin_End, DSBin_Pos;
                 float *pdata;
                 // Need to output the DS processed LDbins to output files
                 //for (vector <LDBin_Object>::iterator it_DSbin = DeepSynthesis.DS_Processed_Bin_Buffer.begin(); it_DSbin !=DeepSynthesis.DS_Processed_Bin_Buffer.end(); it_DSbin++)
                 while(!DeepSynthesis.DS_Processed_Bin_Buffer.empty())
                 {
                     LDBin_Object DSbin= DeepSynthesis.DS_Processed_Bin_Buffer.back();
                     DSBin_Start= DSbin.LDBin_Start;
                     DSBin_End  = DSbin.LDBin_End;
                     DSBin_Pos  = DSbin.LDBin_Pos;
                     pdata      = DSbin.LDBinData;

                     O_BinMap_File << DSBin_Start << out_delim << DSBin_End << out_delim << DSBin_Pos <<endl;
                     int i_individual;
                     for (i_individual =0; i_individual<Individual_Num-1; i_individual++)
                     {
                        O_G_File<< *(pdata+i_individual)<< out_delim;
                     }
                     O_G_File<< *(pdata+i_individual)<< endl;

                     //delete the data pointer, which is allocated during the input reading
                     if(pdata!=NULL)
                     {
                        delete pdata;
                        pdata=NULL;
                     }

                     DeepSynthesis.DS_Processed_Bin_Buffer.pop_back();
                     DS_Line_Count ++;
                 }

                 DeepSynthesis.DS_Processed_Bin_Buffer.clear();
             }

		  }

		  cout <<"#Process Line= " << Line_Count <<" and #DS Line= " << DS_Line_Count <<endl;

		  Line_Count ++;
	   }

	   assert(Bin_Count_Map==Line_Count);

	   /*After the Reading processing finished, there are some linebins stored in DS object buffers, we need a flush processing*/
	   DeepSynthesis.Flush_process();
	   /*Finally, we flush all the intermediated linebins to the DS_Processed_Bin Buffers, we need to output them to the output files*/
       if(DeepSynthesis.DS_Processed_Bin_Buffer.size()>0)
       {
            int DSBin_Start, DSBin_End, DSBin_Pos;
            float *pdata;
            // Need to output the DS processed LDbins to output files
           // for (vector <LDBin_Object>::iterator it_DSbin = DeepSynthesis.DS_Processed_Bin_Buffer.begin(); it_DSbin !=DeepSynthesis.DS_Processed_Bin_Buffer.end(); it_DSbin++)
            while(!DeepSynthesis.DS_Processed_Bin_Buffer.empty())
            {
                LDBin_Object DSbin= DeepSynthesis.DS_Processed_Bin_Buffer.back();
                DSBin_Start= DSbin.LDBin_Start;
                DSBin_End  = DSbin.LDBin_End;
                DSBin_Pos  = DSbin.LDBin_Pos;
                pdata =DSbin.LDBinData;
                O_BinMap_File << DSBin_Start << out_delim << DSBin_End << out_delim << DSBin_Pos <<endl;
                int i_individual;
                for (i_individual =0; i_individual<Individual_Num-1; i_individual++)
                {
                    O_G_File<< *(pdata+i_individual)<< out_delim;
                }
                O_G_File<< *(pdata+i_individual)<< endl;

               //delete the data pointer, which is allocated during the input reading
                if(pdata!=NULL)
                {
                   delete pdata;
                   pdata=NULL;
                }

                DeepSynthesis.DS_Processed_Bin_Buffer.pop_back();
                DS_Line_Count ++;
             }
        }
        DeepSynthesis.DS_Processed_Bin_Buffer.clear();
        cout <<"#DS Line= " << DS_Line_Count <<endl;
    }

    return 0;

}
