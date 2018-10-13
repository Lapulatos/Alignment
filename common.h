#ifndef __COMMON_H__
#define __COMMON_H__

#include <memory>
#include <limits>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <map>
#include <memory.h>

#include "sequence.h"

// type defination.
typedef std::map<std::string, Sequence> SequenceMap;
typedef std::pair<Sequence, Sequence> CompareSeqPair;
typedef std::pair<std::string, Sequence> SearchPair;
typedef std::vector<CompareSeqPair> ComparePairList;
typedef std::pair<std::string, std::string> SeqPair;


//constexpr int NEGINF = INT_MIN;
constexpr int NEGINF = -10000;

#endif		// __COMMON_H__
