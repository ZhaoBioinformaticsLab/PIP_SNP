/*Author: Wenchao Zhang, Noble Research Institute
  Written on 02/24/2020                          */

#include <vector>   // For using Vector
#include <math.h>
#include "Common.h"

bool compareLDBin(LDBin_Object val, LDBin_Object p)
{
   return val.LDBin_Pos > p.LDBin_Pos;
}

/* Calculation of Linkage disequilibrium measure D;
 D equals to p11*p22- p12*p21, which can be prove it by solving the four equations of Genetic Theory.
 Analysis of linkage disequilbrium (LD) between polymorphic sites in a locus identifies "clusters" of
 highly correlated sites based on the r2LD statistic. Sites are binned into sets of highly informative
  markers to minimize redundant data. This data is useful for the development of a minimal set of SNPs
  which can be used for large-scale genotyping of similar sample populations  */
float Calculate_Pairwise_LD(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num, int ABS_Mode)
{
     float D;
     int count_arry[3][3] ={0};
     int c= 0;
     for (int i=0; i< Individual_Num; i++)
     {
         char g1 = Genotype_Array1[i];
         char g2 = Genotype_Array2[i];
         if(g1 >=0 && g2>=0)
         {
           count_arry[g1][g2]++;
           c++;
         }
     }

     int tota = count_arry[1][0] + count_arry[1][1] + count_arry[1][2] + 2 * (count_arry[2][0] + count_arry[2][1] + count_arry[2][2]);
     double meana = 1.0*tota/(1.0*c);

     int totb = count_arry[0][1] + count_arry[1][1] + count_arry[2][1] + 2 * (count_arry[0][2] + count_arry[1][2] + count_arry[2][2]);
     double meanb = 1.0*totb/(1.0*c);

     double xy = 0.0;
     double xx = 0.0;
     double yy = 0.0;

     for (int i = 0; i < 3; i++)
     {
        for (int j = 0; j < 3; j++)
        {
            xy += count_arry[i][j] * (i - meana) * (j - meanb);
            xx += count_arry[i][j] * (i - meana) * (i - meana);
            yy += count_arry[i][j] * (j - meanb) * (j - meanb);
        }
     }

      if((xx == 0.0) || (yy == 0.0))
          D =0.0;
      else
          D = (xy * xy) / (xx * yy);

      if(ABS_Mode==0)
         return D;
      else
         return fabs(D);
}
/*Override function for float data type*/
float Calculate_Pairwise_LD(float *Genotype_Array1, float *Genotype_Array2, int Individual_Num, int ABS_Mode)
{
     float D;
     int count_arry[3][3] ={0};
     int c= 0;
     for (int i=0; i< Individual_Num; i++)
     {
         char g1 = Genotype_Array1[i];
         char g2 = Genotype_Array2[i];
         if(g1 >=0 && g2>=0)
         {
           count_arry[g1][g2]++;
           c++;
         }
     }

     int tota = count_arry[1][0] + count_arry[1][1] + count_arry[1][2] + 2 * (count_arry[2][0] + count_arry[2][1] + count_arry[2][2]);
     double meana = 1.0*tota/(1.0*c);

     int totb = count_arry[0][1] + count_arry[1][1] + count_arry[2][1] + 2 * (count_arry[0][2] + count_arry[1][2] + count_arry[2][2]);
     double meanb = 1.0*totb/(1.0*c);

     double xy = 0.0;
     double xx = 0.0;
     double yy = 0.0;

     for (int i = 0; i < 3; i++)
     {
        for (int j = 0; j < 3; j++)
        {
            xy += count_arry[i][j] * (i - meana) * (j - meanb);
            xx += count_arry[i][j] * (i - meana) * (i - meana);
            yy += count_arry[i][j] * (j - meanb) * (j - meanb);
        }
     }

      if((xx == 0.0) || (yy == 0.0))
          D =0.0;
      else
          D = (xy * xy) / (xx * yy);

      if(ABS_Mode==0)
         return D;
      else
         return fabs(D);
}

// Use the overide function to accept the representative char type genotype as 0,1,2
float Calculate_Pairwise_PCC(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num, int ABS_Mode)
{
    int i;
    int c=0;
    float corr;
    float mean_a =0.0;
    float mean_b =0.0;
    for (i=0; i< Individual_Num ; i++)
    {
         char g1 = Genotype_Array1[i];
         char g2 = Genotype_Array2[i];
         if(g1 >=0 && g2>=0)
         {
             mean_a += g1;
             mean_b += g2;
             c++;
         }
    }

    mean_a =mean_a/c;
    mean_b =mean_b/c;

    float covab =0.0;
    float var_a =0.0;
    float var_b =0.0;

    for ( i=0; i< Individual_Num ; i++)
    {
        char g1 = Genotype_Array1[i];
        char g2 = Genotype_Array2[i];
        if(g1 >=0 && g2>=0)
        {
            covab +=  (g1- mean_a) *(g2- mean_b);
            var_a +=  (g1- mean_a) *(g1- mean_a);
            var_b +=  (g2- mean_b) *(g2- mean_b);
        }
    }

    if ((var_a ==0) ||(var_b ==0))
         corr =0.0;
    else
         corr = covab/ sqrt(var_a * var_b);

    if(ABS_Mode==0)
      return corr;
    else
      return fabs(corr);

}

// Use the overide function to accept the synthesised float type genotype
float Calculate_Pairwise_PCC(float *Genotype_Array1, float *Genotype_Array2, int Individual_Num, int ABS_Mode)
{
    int i;
    int c=0;
    float corr;
    float mean_a =0.0;
    float mean_b =0.0;
    for (i=0; i< Individual_Num ; i++)
    {
         float g1 = Genotype_Array1[i];
         float g2 = Genotype_Array2[i];

         mean_a += g1;
         mean_b += g2;
         c++;

    }

    mean_a =mean_a/c;
    mean_b =mean_b/c;

    float covab =0.0;
    float var_a =0.0;
    float var_b =0.0;

    for ( i=0; i< Individual_Num ; i++)
    {
        float g1 = Genotype_Array1[i];
        float g2 = Genotype_Array2[i];

        covab +=  (g1- mean_a) *(g2- mean_b);
        var_a +=  (g1- mean_a) *(g1- mean_a);
        var_b +=  (g2- mean_b) *(g2- mean_b);

    }

    if ((var_a ==0) ||(var_b ==0))
         corr =0.0;
    else
         corr = covab/ sqrt(var_a * var_b);

    if(ABS_Mode==0)
      return corr;
    else
      return fabs(corr);

}

void Synthesis_Process(float *pdata_Buffer, float* pdata_Synthesised, int Individual_Num, int Bin_Num, int *Pos, int LDBin_Corr_Method, int LDBin_Synthesis_Method)
{
   switch (LDBin_Synthesis_Method)
   {
      case Representative_Synthesis:// 2 Indicate need to find the optimal bin as the representative
         {
           float *Lpldbin_sum_corr= new float [Bin_Num]; // (1xm ) To store the sum of correlation coefficients of a bin to other bin
           for (int i_bin=0; i_bin < Bin_Num; i_bin++ )
           {
               /* Intialize each sum as 0.0 */
               Lpldbin_sum_corr[i_bin] =0.0;
           }
           calculate_pairwise_corr_sum(Lpldbin_sum_corr, pdata_Buffer, Individual_Num, Bin_Num, LDBin_Corr_Method);
           int opt_bin   =0;
           float opt_cor_sum=0.0;
           for (int j_bin=0; j_bin < Bin_Num; j_bin++ )
           {
               if (opt_cor_sum < Lpldbin_sum_corr[j_bin])
               {
                   opt_cor_sum = Lpldbin_sum_corr[j_bin];
                   opt_bin  = j_bin;
               }
           }

           *Pos = opt_bin;  // Recorded the position index of the Bin

           delete Lpldbin_sum_corr;
           for (int i_individual=0; i_individual<Individual_Num; i_individual++)
           {
              pdata_Synthesised[i_individual] = *(pdata_Buffer+opt_bin*Individual_Num+i_individual);
           }

        }
        break;
      case Norm_Synthesis: // 1. Indicate need to calculate the norm of multiple SNP marker values as the synthesized marker
        {
           for(int i_individual=0; i_individual<Individual_Num; i_individual++)
           {
                float value_synthesis = 0.0;
                for (int i_bin=0; i_bin< Bin_Num; i_bin++)
                {
                    float temp_value=pdata_Buffer[i_bin*Individual_Num+i_individual];
                    value_synthesis += temp_value*temp_value;
                }
                pdata_Synthesised[i_individual] =sqrt(value_synthesis);
           }
           *Pos = Bin_Num/2;
        }
        break;
      default:
        {
           for(int i_individual=0; i_individual<Individual_Num; i_individual++)
           {
                float value_synthesis = 0.0;
                for (int i_bin=0; i_bin< Bin_Num; i_bin++)
                {
                    float temp_value=pdata_Buffer[i_bin*Individual_Num+i_individual];
                    value_synthesis += temp_value*temp_value;
                }
                pdata_Synthesised[i_individual] =sqrt(value_synthesis);
           }
           *Pos = Bin_Num/2;
        }
        break;
   }

   return;
}

void calculate_pairwise_corr_sum(float *pLDBin_sum_corr, float *pdata_Buffer, int Individual_Num, int Bin_Num, int method_corr)
{
    for (int i_bin=0; i_bin< Bin_Num; i_bin++)
    {
        float *bin_one =  pdata_Buffer+i_bin*Individual_Num;
        for (int j_bin=i_bin+1; j_bin< Bin_Num; j_bin++)
        {
           float *bin_two = pdata_Buffer+j_bin*Individual_Num;
           float temp_r2  = Calculate_r2_one_bin_to_another_bin(bin_one, bin_two, Individual_Num, method_corr);
           pLDBin_sum_corr[i_bin] += fabs(temp_r2);
           pLDBin_sum_corr[j_bin] += fabs(temp_r2);
        }
    }
    return;
}

float Calculate_r2_one_bin_to_another_bin(float *bin_one, float *bin_two, int Individual_Num, int method_corr)
{
    float r2_value =0.0;
    switch (method_corr)
    {
        case  Pairwise_LD:
            r2_value+= Calculate_Pairwise_LD(bin_two, bin_one, Individual_Num, ABS_CORR);
            break;
        case Pairwise_PCC:
            r2_value+= Calculate_Pairwise_PCC(bin_two, bin_one, Individual_Num, ABS_CORR);
            break;
        default :
            r2_value+= Calculate_Pairwise_PCC(bin_two, bin_one, Individual_Num, ABS_CORR);
            break;
    }
    return r2_value;
}
