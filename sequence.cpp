#include "sequence.h"

Sequence::Sequence() :
	m_SeqName_(""), m_Species_(""), m_Sequence_("")
{}

Sequence::Sequence(const std::string& seqName_, const std::string& species_, const std::string& sequence_) :
	m_SeqName_(seqName_), m_Species_(species_), m_Sequence_(sequence_)
{}

const std::string& Sequence::GetSequenceName() const
{
	return m_SeqName_;
}

const std::string& Sequence::GetSpecies() const
{
	return m_Species_;
}

const std::string& Sequence::GetSequence() const
{
	return m_Sequence_;
}

void Sequence::SetSequenceName(const std::string& seqName_)
{
	m_SeqName_ = seqName_;
}

void Sequence::SetSpecies(const std::string& species_)
{
	m_Species_ = species_;
}

void Sequence::SetSequence(const std::string& sequence_)
{
	m_Sequence_ = sequence_;
}

const bool Sequence::operator== (const Sequence& other_) const
{
	if(m_SeqName_ == other_.m_SeqName_ && m_Species_ == other_.m_Species_)
		return true;
	else
		return false;
}

const void Sequence::Print() const
{
	std::cout << "SeqName: " << m_SeqName_ << "  Species: " << m_Species_ << "  Length: " << m_Sequence_.length() << std::endl;
	std::cout << "Sequence: " << m_Sequence_ << std::endl;
}
