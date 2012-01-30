/*
 * Imported directly from Dprime by the WFUBMC group.
 * 
 * There were no comments when I got this, so any changes
 * are marked by the presence of comments. 
 * 
 */

# ifndef HAPLOTYPE_H
# define HAPLOTYPE_H

# include <vector>
# include <iomanip>

using namespace std;

class haplotype
{
	public:
		haplotype(); 
		haplotype(int *currAlleles, int numMarkers, int *numAlleles);
		haplotype(int *currAlleles, vector<int> &numAlleles);
		haplotype(int haplotypeIndex, int numMarkers, int *numAlleles);
		haplotype(int haplotypeIndex, vector<int> &numAlleles);
		~haplotype();
		
		void init(const vector<int> &hapAlleles);
		void setHaplotypeIndex(int *numAlleles);
		void setHaplotypeIndex(vector<int> &numAlleles);
		int alleles2Index(int *numAlleles) const;
		int alleles2Index(vector<int> &numAlleles) const;
		int index2Allele(int imark, int numMarkers, int *numAlleles) const;
		int index2Allele(int imark, vector<int> &numAlleles) const;
		
		int getNumDescendents() const { return numDescendents; }
		int getDescendent(int idesc) const { return descendents[idesc]; }
		int getHaplotypeIndex() const { return haplotypeIndex; }
		
		int getAllele(int imark) const { return alleles[imark]; }
		
		int getNumMarkers() const { return static_cast<int>(alleles.size()); }
		
		void setAllele(int markerNum, int allele){ alleles[markerNum] = static_cast<char>(allele); }
		
		class alleleNumOutsideRangeEx {};
		
	private:
		void fillinDescendents(int numMarkers, int *numAlleles);
		void fillinDescendents(vector<int> &numAlleles);
		
		int haplotypeIndex;
		int numDescendents;
		int *descendents;
		vector<char> alleles;
};
# endif
