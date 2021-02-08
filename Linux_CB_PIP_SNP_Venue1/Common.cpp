/*Author: Wenchao Zhang, Noble Research Institute
  Written on 12/13/2018                          */

#include <vector>   // For using Vector
#include <math.h>
#include  "Common.h"

/* Free all the pointers that stored in the group_marker_lines */
void destroy_grouplines(vector <char *> &group_marker_lines)
{
    while(!group_marker_lines.empty())
    {
        delete group_marker_lines.back();
        group_marker_lines.pop_back();
    }
}

/* Calculate the R2 values between line one to another marker line*/
float Calculate_r2_one_line_to_another_line(char *line_one, char *line_two, int Individual_Num, int Method_Corr)
{
    float r2_value =0.0;
    switch (Method_Corr)
    {
        case  Pairwise_LD:
            r2_value+= Calculate_Pairwise_LD(line_two, line_one, Individual_Num);
            break;
        case Pairwise_PCC:
            r2_value+= Calculate_Pairwise_PCC(line_two, line_one, Individual_Num);
            break;
        default :
            r2_value+= Calculate_Pairwise_PCC(line_two, line_one, Individual_Num);
            break;
    }
    return r2_value;
}

/* Calculate the R2 values between the current line to the group marker line */
void Calculate_r2_one_line_to_grouplines(vector <char *> group_marker_lines, char *line_one, int Individual_Num, int Method_Corr, float *R2_Average, float *R2_Boundary)
{
    float r2_value;
    float r2_sum     =0.0;
    int marker_num   =group_marker_lines.size();
    int marker_count =0;
    for ( vector<char *>::iterator it_marker_line=  group_marker_lines.begin(); it_marker_line != group_marker_lines.end(); it_marker_line++)
    {
        switch (Method_Corr)
        {
           case  Pairwise_LD:
              r2_value = Calculate_Pairwise_LD(*it_marker_line, line_one, Individual_Num);
              break;
           case Pairwise_PCC:
              r2_value = Calculate_Pairwise_PCC(*it_marker_line, line_one, Individual_Num);
              break;
           default :
              r2_value = Calculate_Pairwise_PCC(*it_marker_line, line_one, Individual_Num);
              break;
        }

        r2_sum  += r2_value;

        marker_count ++;
        if(marker_count ==marker_num) *R2_Boundary = r2_value;

    }

    *R2_Average =r2_sum/marker_num;
    return;
}

/* Calculation of Linkage disequilibrium measure D;
 D equals to p11*p22- p12*p21, which can be prove it by solving the four equations of Genetic Theory.
 Analysis of linkage disequilbrium (LD) between polymorphic sites in a locus identifies "clusters" of
 highly correlated sites based on the r2LD statistic. Sites are binned into sets of highly informative
  markers to minimize redundant data. This data is useful for the development of a minimal set of SNPs
  which can be used for large-scale genotyping of similar sample populations  */

float Calculate_Pairwise_LD(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num)
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

      return D;

}

/*Calculate the Pearson correlation coefficient */
float Calculate_Pairwise_PCC(char *Genotype_Array1, char *Genotype_Array2, int Individual_Num)
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

    return corr;

}

bool Comp_Individual_Pair_by_Distance(Individual_Pair_Distance_KnownGenotype val, Individual_Pair_Distance_KnownGenotype p)
{
    return val.Distance < p.Distance;
}
