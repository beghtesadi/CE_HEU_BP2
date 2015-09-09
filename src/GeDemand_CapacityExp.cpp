#include <iostream>
#include <math.h>
#include <vector>

#include "GeDemand_CapacityExp.h"

using namespace std;

void GenerateDemand (vector<int> &x , vector<double> &y, double mean, double K){
	int a, b;
        a = CalMin(mean,ceil((K-1)/2));
        if (mean > ceil ((K-1)/2))
           b= floor((K-1)/2);
        else
           b= K-mean-1;
        double prob = 0;
        for (int i= 0; i< mean-a +1 ; i++){
          prob = prob +CalProb (i , mean);
          }
        double probb = prob;
        x.push_back(mean-a);
        //y.push_back(CalProb ( mean-a, mean));
        y.push_back(probb);
	for (int i=1 ; i<a+b ; i++){
           y.push_back(CalProb(mean-a+i , mean)); 
	       x.push_back(mean-a+i);
           probb = probb +y[i];}
        x.push_back (mean+b);
        //y.push_back(CalProb ( mean+b, mean));
        y.push_back(1-probb);
}
double CalProb (int n, double mean) 
 {
  double prob=exp(-1*mean);
  if (n > 0)
  {
    for (int i=1 ; i<n+1 ; i++)
    {
     prob = prob * mean / i ;
    }
   } 
  return (prob);
 }

int CalMin (int  a, int b)
{
   int min;
   if (a < b)
     min=a;
   else
     min=b;
   return (min);
}

