/*Author: Wenchao Zhang, Noble Research Institute
  Written on 12/31/2018                          */
#include "CImputingMissing.h"
#include <math.h>
#include <float.h>
#include <iostream>  // Cin/ Cout

CImputingMissing::CImputingMissing(CLDBin_DataObject &ldbin_obj, int K)
{
    K_NN                          =K;
    Method_Corr                   =ldbin_obj.LDBin_Method_Corr;
    LDBin_Start                   =ldbin_obj.LDBin_Start;
    LDBin_End                     =ldbin_obj.LDBin_Current;
    LDBin_Marker_Num              =LDBin_End-LDBin_Start+1; // ldbin_obj.Get_markerlines_number();
    LDBin_Pos                     =(LDBin_Marker_Num-1)/2;
    LDBin_Individual_Num          =ldbin_obj.LDBin_Individual_Num;

    Lpldbin_OriginalData_Transpose= new char [LDBin_Marker_Num*LDBin_Individual_Num];
    Lpldbin_ImputedData           = new char [LDBin_Marker_Num*LDBin_Individual_Num];
    Lpldbin_Synthesized           = NULL;

    int marker_line_count =0;
    for (vector<char *>::iterator it_markerline =ldbin_obj.LDBin_marker_lines.begin(); it_markerline !=ldbin_obj.LDBin_marker_lines.end(); it_markerline++)
    {
        char *current_markerline = *(it_markerline);
        for (int individual_count =0; individual_count< LDBin_Individual_Num; individual_count++)
        {
            char current_genotype = *(current_markerline+individual_count);
            Lpldbin_OriginalData_Transpose[individual_count*LDBin_Marker_Num+ marker_line_count] = current_genotype;
            Lpldbin_ImputedData[marker_line_count*LDBin_Individual_Num + individual_count]       = current_genotype;

            if(current_genotype<0) // 0,1,2 means the meaningful known genotype while -1 corresponds to missing genotype value, and need to be imputed
            {
               Genotype_Positon missing_genotype_pos;
               missing_genotype_pos.Marker_Pos     =marker_line_count;
               missing_genotype_pos.Individual_Pos =individual_count;
               Missing_Genotype_Vector.push_back(missing_genotype_pos);
            }
        }

        marker_line_count++;
    }

   return;
}

/* Can be implemented by GPU in parallel*/
float CImputingMissing::calcuate_distance_two_indivudals(char *data_individual_1, char *data_individual_2, int marker_len)
{
     float sum_dis =0.0;
     int count_num =0;
     for (int i_marker= 0; i_marker < marker_len; i_marker++ )
     {
         char gen1= data_individual_1[i_marker];
         char gen2= data_individual_2[i_marker];
         if((gen1!=-1)&&(gen2 !=-1))
         {
            sum_dis +=fabs(gen1-gen2);
            count_num++;
         }
     }

     if(count_num==0)
     {
         return DBL_MAX;
     }
     else
     {
         return (sum_dis+ ADD_CONSTANT)/count_num;
     }
}

char CImputingMissing::Impute_one(int missing_individual_pos, int missing_marker_pos)
{
     char ret_imputed;

     char *data_individual_1, *data_individual_2;

     data_individual_1 = Lpldbin_OriginalData_Transpose+missing_individual_pos*LDBin_Marker_Num;
     vector <Individual_Pair_Distance_KnownGenotype> Individual_pair_distance_vector;
     for (int individual_pos=0; individual_pos< LDBin_Individual_Num; individual_pos++)
     {
         char Knowngenotype        = Lpldbin_OriginalData_Transpose[individual_pos*LDBin_Marker_Num+missing_marker_pos];
         if ((individual_pos != missing_individual_pos)&&(Knowngenotype>-1))
         {
             Individual_Pair_Distance_KnownGenotype pair_distance;
             pair_distance.X_Individual = missing_individual_pos;
             pair_distance.Y_Individual = individual_pos;
             data_individual_2          = Lpldbin_OriginalData_Transpose + individual_pos*LDBin_Marker_Num;
             pair_distance.Distance     = calcuate_distance_two_indivudals( data_individual_1, data_individual_2, LDBin_Marker_Num);
             pair_distance.KnownGenotype= Knowngenotype;

             Individual_pair_distance_vector.push_back(pair_distance);
         }
     }

     sort(Individual_pair_distance_vector.begin(), Individual_pair_distance_vector.end(), Comp_Individual_Pair_by_Distance);

     int f_k =0;   //used to count the efficient individual fop imputing.
     float weight[3];  // used to stored the weighed sum of three biological meaningful genotype values: 0,1,2
     weight[0] =0.0;
     weight[1] =0.0;
     weight[2] =0.0;
     for ( vector <Individual_Pair_Distance_KnownGenotype>::iterator It_pair_distance = Individual_pair_distance_vector.begin(); (It_pair_distance!=Individual_pair_distance_vector.end())&&(f_k< K_NN); It_pair_distance++)
     {
        /*cout << "f_k=" << f_k <<"The Corresponding Distance is" << (*It_pair_distance).Distance << "1.0/Distance is" << 1.0/(*It_pair_distance).Distance <<endl;
          cout << "weight[0]=" << weight[0] <<"weight[1]=" << weight[1] << "weight[2]=" << weight[2] <<endl; */
          switch((*It_pair_distance).KnownGenotype)
          {
              case 0:
                weight[0] += 1.0/((*It_pair_distance).Distance);
                break;
              case 1:
                weight[1] += 1.0/((*It_pair_distance).Distance);
                break;
              case 2:
                weight[2] += 1.0/((*It_pair_distance).Distance);
                break;
             default:
                break;
          }

         f_k++;
     }

    if((weight[0]>=weight[1])&&(weight[0]>=weight[2]))
    {
       ret_imputed =0;
    }
    else if(weight[1]>=weight[2])
    {
        ret_imputed =1;
    }
    else
    {
        ret_imputed =2;
    }

    return ret_imputed;

}

void CImputingMissing::Impute_all()
{
    if(!Missing_Genotype_Vector.empty())
    {
        for( vector<Genotype_Positon>::iterator It_MissingGenotype= Missing_Genotype_Vector.begin(); It_MissingGenotype!= Missing_Genotype_Vector.end(); It_MissingGenotype++)
        {
           int  individual_pos                                                         = (*It_MissingGenotype).Individual_Pos;
           int  marker_pos                                                             = (*It_MissingGenotype).Marker_Pos;
           char imputed_genotype                                                       = Impute_one(individual_pos, marker_pos);
    //     cout << "Missing genotype at individual_pos=" << individual_pos << "marker_pos=" << marker_pos << "imputed genotype=" << +imputed_genotype<< endl;
           Lpldbin_ImputedData[marker_pos*LDBin_Individual_Num + individual_pos]       = imputed_genotype;
        }
    }
}

void CImputingMissing::calculate_pairwise_corr_sum(float *Lpldbin_marker_sum_corr, int method_corr)
{
    for (int i_line=0; i_line< LDBin_Marker_Num; i_line++)
    {
        char *line_one =  Lpldbin_ImputedData+i_line*LDBin_Individual_Num;
        for (int j_line=i_line+1; j_line< LDBin_Marker_Num; j_line++)
        {
           char *line_two = Lpldbin_ImputedData+j_line*LDBin_Individual_Num;
           float temp_r2  = Calculate_r2_one_line_to_another_line(line_one, line_two, LDBin_Individual_Num, method_corr);
           Lpldbin_marker_sum_corr[i_line] += temp_r2;
           Lpldbin_marker_sum_corr[j_line] += temp_r2;
        }
    }
}

void CImputingMissing::Synthesize(int syn_mode)
{
   if (syn_mode>0)
   {
       /*Indicate that it's necessary to synthesize the detected LDBin */
       Lpldbin_Synthesized =new float[LDBin_Individual_Num];
       if(syn_mode==2) // Indicate need to find the optimal SNP marker as the representative synthesized marker of the detected LDBin
       {
           float *Lpldbin_marker_sum_corr= new float [LDBin_Marker_Num]; // (1xm ) To store the sum of correlation coefficients of a marker to other markers
           for (int i_marker=0; i_marker < LDBin_Marker_Num; i_marker++ )
           {
               /* Intialize each sum as 0.0 */
               Lpldbin_marker_sum_corr[i_marker] =0.0;
           }
           calculate_pairwise_corr_sum(Lpldbin_marker_sum_corr, Method_Corr);
           int opt_marker   =0;
           float opt_cor_sum=0.0;
           for (int j_marker=0; j_marker < LDBin_Marker_Num; j_marker++ )
           {
               if (opt_cor_sum < Lpldbin_marker_sum_corr[j_marker])
               {
                   opt_cor_sum = Lpldbin_marker_sum_corr[j_marker];
                   opt_marker  = j_marker;
               }
           }

           LDBin_Pos = opt_marker;  // Recorded the position index of the LD Bin

           delete Lpldbin_marker_sum_corr;
           for (int i_individual=0; i_individual<LDBin_Individual_Num; i_individual++)
           Lpldbin_Synthesized[i_individual] = *(Lpldbin_ImputedData+opt_marker*LDBin_Individual_Num+i_individual);
       }
       else    // Indicate need to calculate the norm of multiple SNP marker values as the synthesized marker of the detected LDBin
       {
           LDBin_Pos = (LDBin_Marker_Num -1)/2;  // Recorded the position index of the LD Bin

           for(int i_individual=0; i_individual<LDBin_Individual_Num; i_individual++)
           {
                float p_synthesis_marker = 0.0;
                for (int i_marker=0; i_marker< LDBin_Marker_Num; i_marker++)
                {
                    char genotype_value=Lpldbin_ImputedData[i_marker*LDBin_Individual_Num+i_individual];
                    p_synthesis_marker += genotype_value*genotype_value;
                }
                Lpldbin_Synthesized[i_individual] =sqrt(p_synthesis_marker);
            }
       }
   }

}


CImputingMissing::~CImputingMissing()
{
    if(Lpldbin_OriginalData_Transpose!=NULL)
    {
        delete Lpldbin_OriginalData_Transpose;
        Lpldbin_OriginalData_Transpose=NULL;
    }

    if(Lpldbin_ImputedData!=NULL)
    {
        delete Lpldbin_ImputedData;
        Lpldbin_ImputedData =NULL;
    }

    if(Lpldbin_Synthesized!=NULL)
    {
        delete  Lpldbin_Synthesized;
        Lpldbin_Synthesized =NULL;
    }
    return;
}
