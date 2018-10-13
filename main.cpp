#include <fstream>
#include <vector>
#include <bitset>
#include <utility>
#include <map>
#include <exception>

#include "sequence.h"
#include "common.h"
#include "utils.h"
#include "solver.hxx"
#include "bloSum62.hxx"

using namespace std;

int main(int argc, char* argv[])
{
	// argument judgement.
	if(argc < 3)
	{
		cout << "usage: program sequence_file blosum62_file" << endl;
		cout << "sequence_file format: " << endl;
		cout << "seqname species	---	[1 line]" << endl;
		cout << "sequence	---	[multiple line]" << endl << endl;
		cout << "notice: each seqname at least appear twice to be compared." << endl;
		cout << "blosum62_file format: see the file blosum62.hxx" << endl;

		return -1;
	}

	// get file names.
	std::string seqFile = std::string(argv[1]);
	std::string bloFile = std::string(argv[2]);
	// check the two parameters.
	if(seqFile.empty() || bloFile.empty())
		return -2;

	// read protein sequence file.
	ComparePairList cmpList;
	GetProteinSequencePairs(seqFile, cmpList);
	cout << "Compare List Size: " << cmpList.size() << endl;	// output the size of sequence pair that need to be analysis.
	
	// read protein BloSum62 file.
	BloSum62<int> bloSum62(bloFile);

	//--------------------------------------------------//
	// affine-gap local alignment.						//
	//--------------------------------------------------//

	// align each sequence pair.
	for(const CompareSeqPair& seqPair : cmpList)
	{
		// the analysis process can not stop when feeding sequence 'CN1A' & '5016' into the solver.
		// maybe there are some bug in analysis process.
		if(seqPair.first.GetSequenceName() != "CN1A" && 
		   seqPair.first.GetSequenceName() != "5016")
		{
			// construct a solver to solve this sequence pair.
			Solver<int> solver(seqPair.first, seqPair.second, bloSum62);
			// update the three temporary matrix in the solver, and get the final score between the two sequences.
			int score = solver.Update();
			// construct the alignment result of the two sequences.
			SeqPair result = solver.Construct();
			// output the alignment result.
			cout << seqPair.first.GetSequenceName() << ": length = " << seqPair.first.GetSequence().length() 
				<< ", score = "<< score << endl;
			cout << "seq1: " << result.first << endl;
			cout << "seq2: " << result.second << endl;
		}
	}

	return 0;
}