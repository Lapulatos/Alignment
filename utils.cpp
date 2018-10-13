#include "utils.h"

void GetProteinSequencePairs(const std::string& seqfile, ComparePairList& list)
{
	// seqfile format:
	//	seqname species		---		[1 line]
	//	sequence			---		[multiple line]
	// notice: in seqfile, there are multiple species' sequences.

	// open the file.
	std::fstream fs(seqfile);

	if(fs.is_open())
	{
		// 
		SequenceMap seqMap;
		SequenceMap::iterator search;
		// 
		std::string line;
		std::string name = "", species = "", sequence = "";
		size_t spaceIndex = 0;

		// read the file.
		while(std::getline(fs, line))
		{
			// analys this line.
			if((spaceIndex = line.find(' ')) != std::string::npos)	// a new sequence.
			{
				// process the old sequence.
				if(!name.empty() && !species.empty() && !sequence.empty())
				{
					// the new sequence is already in search map, get the old one & construct the new one, and put them in the list.
					if((search = seqMap.find(name)) != seqMap.end())
					{
						list.push_back(CompareSeqPair(search->second, Sequence(name, species, sequence)));
					} 
					else 
					{
						// the new sequence is not in the search map, just insert it in the search map.
						seqMap.insert(SearchPair(name, Sequence(name, species, sequence)));
					}
					// clear the old sequence.
					sequence = "";
				}

				// process the new sequence.
				// get the information of this sequence.
				name = line.substr(0, spaceIndex);
				species = line.substr(spaceIndex + 1, line.length() - spaceIndex - 1);
			} 
			else		// just a protein sequence, append it to the sequence.
			{
				sequence += line;
			}
			// the previous sequence could be constructed.
		}

		// process the last sequence.
		if((search = seqMap.find(name)) != seqMap.end())
		{
			list.push_back(CompareSeqPair(search->second, Sequence(name, species, sequence)));
			sequence = "";
		}

		// close the file.
		fs.close();
	}
}
