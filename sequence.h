#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <string>
#include <iostream>

class Sequence
{
public:
	Sequence();
	Sequence(const std::string& seqName_, const std::string& species_, const std::string& sequence_);

	const std::string& GetSequenceName() const;
	const std::string& GetSpecies() const;
	const std::string& GetSequence() const;
	void SetSequenceName(const std::string& seqName_);
	void SetSpecies(const std::string& species_);
	void SetSequence(const std::string& sequence_);

	const void Print() const;

	const bool operator== (const Sequence& other_) const;

private:
	std::string m_SeqName_;
	std::string m_Species_;
	std::string m_Sequence_;
};

#endif		// __SEQUENCE_H__
