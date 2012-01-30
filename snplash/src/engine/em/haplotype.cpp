# include "haplotype.h"
# include <iostream>
# include <iomanip>
# include <fstream>
  
haplotype::haplotype()
{
   descendents = NULL;
   numDescendents = 0;
   haplotypeIndex = 0;
}
 
haplotype::haplotype(int *currAlleles, int numMarkers, int *numAlleles)
{
   descendents = NULL;
   int imark;
   for(imark = 0; imark < numMarkers; imark++){
      if((currAlleles[imark] < 128) && (currAlleles[imark] > 0)){
         alleles.push_back(static_cast<char>(currAlleles[imark]));
      }
      else{
         throw alleleNumOutsideRangeEx();
      }
   }
   haplotypeIndex = alleles2Index(numAlleles);
   fillinDescendents(numMarkers, numAlleles);
}

haplotype::haplotype(int *currAlleles, vector<int> &numAlleles)
{
   descendents = NULL;
   int imark;
   for(imark = 0; imark < static_cast<int>(numAlleles.size()); imark++){
      if((currAlleles[imark] < 128) && (currAlleles[imark] > 0)){
         alleles.push_back(static_cast<char>(currAlleles[imark]));
      }
      else{
         throw alleleNumOutsideRangeEx();
      }
   }
   haplotypeIndex = alleles2Index(numAlleles);
   fillinDescendents(numAlleles);
}

haplotype::haplotype(int haplotypeIndex, int numMarkers, int *numAlleles)
          :haplotypeIndex(haplotypeIndex)
{
   descendents = NULL;
   int imark;
   for(imark = 0; imark < numMarkers; imark++){
      alleles.push_back(static_cast<char>(index2Allele(imark, numMarkers, numAlleles)));
   }
   fillinDescendents(numMarkers, numAlleles);
}
    
haplotype::haplotype(int haplotypeIndex, vector<int> &numAlleles)
          :haplotypeIndex(haplotypeIndex)
{
   descendents = NULL;
   int imark;
   for(imark = 0; imark < static_cast<int>(numAlleles.size()); imark++){
      alleles.push_back(static_cast<char>(index2Allele(imark, numAlleles)));
   }
   fillinDescendents(numAlleles);
}
    
haplotype::~haplotype()
{
   delete [] descendents;
}

void haplotype::init(const vector<int> &currAlleles)
{
   
   alleles.reserve(currAlleles.size());
   for (unsigned int imark = 0; imark < currAlleles.size(); ++imark){
      alleles.push_back(static_cast<char>(currAlleles[imark]));
   }
}
    
void haplotype::setHaplotypeIndex(int *numAlleles)
{
   haplotypeIndex = alleles2Index(numAlleles);
}

void haplotype::setHaplotypeIndex(vector<int> &numAlleles)
{
   haplotypeIndex = alleles2Index(numAlleles);
}

int haplotype::alleles2Index(int *numAlleles) const
{
   int imark;
   int hapIndex = 0;

   for (imark = 0; imark < static_cast<int>(alleles.size()); ++imark){
      hapIndex = (numAlleles[imark] + 1)*hapIndex + static_cast<int>(alleles[imark]);
   }
   return hapIndex;
}

int haplotype::index2Allele(int imark, int numMarkers, int *numAlleles) const
{
   int jmark;
   int product = 1;

   for (jmark = numMarkers - 1 ; jmark > imark ; --jmark){
      product *= (numAlleles[jmark]+1);
   }
   
   return static_cast<int>((haplotypeIndex/product) % (numAlleles[imark] + 1));
}

int haplotype::alleles2Index(vector<int> &numAlleles) const
{
   int imark;
   int hapIndex = 0;

   for (imark = 0; imark < static_cast<int>(alleles.size()); ++imark){
      hapIndex = (numAlleles[imark] + 1)*hapIndex + static_cast<int>(alleles[imark]);
   }
   return hapIndex;
}

int haplotype::index2Allele(int imark, vector<int> &numAlleles) const
{
   int jmark;
   int product = 1;

   for (jmark = static_cast<int>(numAlleles.size()) - 1 ; jmark > imark ; --jmark){
      product *= (numAlleles[jmark]+1);
   }
   return static_cast<int>((haplotypeIndex/product) % (numAlleles[imark] + 1));
}

void haplotype::fillinDescendents(int numMarkers, int *numAlleles)
{
   int firstZero = numMarkers;
   int imark;
   for(imark = 0; imark < numMarkers; ++imark){
      if(static_cast<int>(alleles[imark]) == 0){
         firstZero = imark;
         break;
      }
   }
   if(firstZero == numMarkers){
      numDescendents = 0;
   }
   else{
      numDescendents = numAlleles[firstZero];
      descendents = new int[numDescendents];

      int idescend;
      for (idescend = 0; idescend < numDescendents; ++idescend){
         alleles[firstZero] = static_cast<char>(idescend + 1);
         descendents[idescend] = alleles2Index(numAlleles);
      }
      alleles[firstZero] = 0;
   }
}

void haplotype::fillinDescendents(vector<int> &numAlleles)
{
   int firstZero = static_cast<int>(numAlleles.size());
   int imark;
   for(imark = 0; imark < static_cast<int>(numAlleles.size()); ++imark){
      if(static_cast<int>(alleles[imark]) == 0){
         firstZero = imark;
         break;
      }
   }
   if(firstZero == static_cast<int>(numAlleles.size())){
      numDescendents = 0;
   }
   else{
      numDescendents = numAlleles[firstZero];
      descendents = new int[numDescendents];

      int idescend;
      for (idescend = 0; idescend < numDescendents; ++idescend){
         alleles[firstZero] = static_cast<char>(idescend + 1);
         descendents[idescend] = alleles2Index(numAlleles);
      }
      alleles[firstZero] = 0;
   }
}

