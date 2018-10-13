#ifndef __BLOSUM62_HXX__
#define __BLOSUM62_HXX__

#include "common.h"

template <typename T>
class BloSum62
{
public:
	BloSum62() :
		m_Matrix_(nullptr), m_Proteins_("")
	{}

	BloSum62(const std::string& blofile) :
		m_Proteins_("")
	{
		// blofile format:
		//		A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V
		//	A	4	-1	-2	-2	0	-1	-1	0	-2	-1	-1	-1	-1	-2	-1	1	0	-3	-2	0
		//	R	-1	.
		//	N	-2		.
		//	D	-2			.
		//	C	0				.
		//	Q	-1					.
		//	E	-1						.
		//	G	0							.		
		//	H	-2								.
		//	I	-1									.
		//	L	-1										.
		//	K	-1											.
		//	M	-1												.
		//	F	-2													.
		//	P	-1														.
		//	S	1															.
		//	T	0																.
		//	W	-3																	.
		//	Y	-2																		.
		//	V	-0																			.
		std::fstream fs(blofile);

		if(fs.is_open())
		{
			// allocate space.
			AllocateSpace();

			// read the file.
			std::string line;
			std::getline(fs, line);
			// read the protein name.
			for(char ch : line)
			{
				if(ch != ' ')
				{
					m_Proteins_ += std::string(1, ch);
				}
			}
			// read the matrix.
			char protein = ' ';
			T value = 0;
			size_t proteinLine = 0, count = 0;
			while(count < s_ProteinNumber_)
			{
				// read a line.
				std::getline(fs, line);
				std::stringstream ss(line);
				// read the protein name and it's line number.
				char protein;
				ss >> protein;
				proteinLine = m_Proteins_.find(protein);
				// read the values correspond to this protein.
				for(size_t i = 0; i < s_ProteinNumber_; ++i)
				{
					ss >> value;
					m_Matrix_[proteinLine][i] = value;
				}
				++count;
			}

			fs.close();
		}
	}

	BloSum62(const T **matrix_, const std::string& proteins_) :
		m_Proteins_(proteins_)
	{
		// allocate new space.
		AllocateSpace();
		// copy data.
		for(size_t i = 0; i < s_ProteinNumber_; ++i)
		{
			memcpy(m_Matrix_[i], matrix_, sizeof(T)*s_ProteinNumber_);
		}
	}

	BloSum62(const BloSum62& other_)
	{
		if(this != &other_)
		{
			// copy protein information.
			this->m_Proteins_ = other_.m_Proteins_;
			// allocate new space.
			this->AllocateSpace();
			// copy the matrix.
			for(size_t i = 0; i < s_ProteinNumber_; ++i)
			{
				memcpy(this->m_Matrix_[i], other_.m_Matrix_[i], sizeof(int)*s_ProteinNumber_);
			}
		}
	}

	~BloSum62()
	{
		ReleaseSpace();
	}

	int** GetBloSum62Matrix() const
	{
		return m_Matrix_;
	}

	std::string GetProteins() const
	{
		return m_Proteins_;
	}

	int GetValue(char p1, char p2) const
	{
		int p1Idx = m_Proteins_.find(p1), p2Idx = m_Proteins_.find(p2);
		assert(p1Idx != std::string::npos && p2Idx != std::string::npos);
		return m_Matrix_[p1Idx][p2Idx];
	}

	void Print() const
	{
		// print header.
		for(size_t i = 0; i < s_ProteinNumber_; ++i)
			std::cout << "\t" << m_Proteins_[i];
		std::cout << std::endl;

		// print body.
		for(size_t i = 0; i < s_ProteinNumber_; ++i)
		{
			// print the name of this protein.
			std::cout << m_Proteins_[i];
			// print the value.
			for(size_t j = 0; j < s_ProteinNumber_; ++j)
			{
				std::cout << "\t" << m_Matrix_[i][j];
			}
			std::cout << std::endl;
		}
	}

private:
	void AllocateSpace()
	{
		// release space.
		ReleaseSpace();
		// allocate new space.
		m_Matrix_ = new int*[s_ProteinNumber_];
		for(size_t i = 0; i < s_ProteinNumber_; ++i)
		{
			m_Matrix_[i] = new int[s_ProteinNumber_];
			memset(m_Matrix_[i], (T)0, sizeof(T)*s_ProteinNumber_);
		}
	}

	void ReleaseSpace()
	{
		if(m_Matrix_ != nullptr)
		{
			for(size_t i = 0; i < s_ProteinNumber_; ++i)
			{
				if(m_Matrix_[i] != nullptr)
				{
					delete[] m_Matrix_[i];
				}
			}
			delete[] m_Matrix_, m_Matrix_ = nullptr;
		}
	}

public:
	static const size_t s_ProteinNumber_;
private:
	T **m_Matrix_;
	std::string m_Proteins_;
};

template <typename T>
const size_t BloSum62<T>::s_ProteinNumber_ = 20;

#endif	// __BLOSUM62_HXX__
