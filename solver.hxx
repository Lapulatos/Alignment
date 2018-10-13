#ifndef __SOLVER_HXX__
#define __SOLVER_HXX__

#include "common.h"
#include "blosum62.hxx"

typedef enum
{
	None,		// from none
	FromA,		// from the matrix A
	FromB,		// from the matrix B
	FromC		// from the matrix C
} eTrace;

template <typename T>
class Space
{
public:
	Space(size_t rows_, size_t cols_) :
		m_Rows_(rows_), m_Cols_(cols_)
	{
		// allocate value space.
		m_Value_ = new T*[m_Rows_];
		for(size_t i = 0; i < m_Rows_; ++i)
		{
			m_Value_[i] = new T[m_Cols_];
			memset(m_Value_[i], (T)0, sizeof(T)*m_Cols_);
		}
		// allocate trace space.
		m_Trace_ = new eTrace*[m_Rows_];
		for(size_t i = 0; i < m_Rows_; ++i)
		{
			m_Trace_[i] = new eTrace[m_Cols_];
			memset(m_Trace_[i], (int)eTrace::None, sizeof(eTrace)*m_Cols_);
		}
	}
	Space(size_t rows_, size_t cols_, const std::initializer_list<T>& values_, const std::initializer_list<eTrace>& traces_) :
		m_Rows_(rows_), m_Cols_(cols_)
	{
		//
		assert((values_.size() == rows_*cols_) && (traces_.size() == rows_*cols_));
		// allocate value space.
		m_Value_ = new T*[m_Rows_];
		for(size_t i = 0; i < m_Rows_; ++i)
		{
			m_Value_[i] = new T[m_Cols_];
			memcpy(m_Value_[i], values_.begin() + m_Cols_*i, sizeof(T)*m_Cols_);
		}
		// allocate trace space.
		m_Trace_ = new eTrace*[m_Rows_];
		for(size_t i = 0; i < m_Rows_; ++i)
		{
			m_Trace_[i] = new eTrace[m_Cols_];
			memcpy(m_Trace_[i], traces_.begin() + m_Cols_*i, sizeof(eTrace)*m_Cols_);
		}
	}
	~Space()
	{
		if(m_Value_ != nullptr)
		{
			for(size_t i = 0; i < m_Rows_; ++i)
			{
				delete[] m_Value_[i];
			}
			delete[] m_Value_;
			m_Value_ = nullptr;
		}
		if(m_Trace_ != nullptr)
		{
			for(size_t i = 0; i < m_Rows_; ++i)
			{
				delete[] m_Trace_[i];
			}
			delete[] m_Trace_;
			m_Trace_ = nullptr;
		}
	}

	void PrintValueMatrix() const
	{
		for(size_t i = 0; i < m_Rows_; ++i)
		{
			for(size_t j = 0; j < m_Cols_; ++j)
			{
				std::cout << m_Value_[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}
	void PrintTraceMatrix() const
	{
		for(size_t i = 0; i < m_Rows_; ++i)
		{
			for(size_t j = 0; j < m_Cols_; ++j)
			{
				if(m_Trace_[i][j] == eTrace::None)
					std::cout << "N\t";
				else if(m_Trace_[i][j] == eTrace::FromA)
					std::cout << "A\t";
				else if(m_Trace_[i][j] == eTrace::FromB)
					std::cout << "B\t";
				else 
					std::cout << "C\t";
			}
			std::cout << std::endl;
		}
	}

	inline size_t GetRows() const
	{
		return m_Rows_;
	}
	inline size_t GetCols() const
	{
		return m_Cols_;
	}

	inline void SetValue(const T& value_, size_t i_, size_t j_)
	{
		assert((0 <= i_ && i_ < m_Rows_) && (0 <= j_ && j_ < m_Cols_));
		m_Value_[i_][j_] = value_;
	}
	inline void SetTrace(const eTrace& trace_, size_t i_, size_t j_)
	{
		assert((0 <= i_ && i_ < m_Rows_) && (0 <= j_ && j_ < m_Cols_));
		m_Trace_[i_][j_] = trace_;
	}
	inline T& GetValue(size_t i_, size_t j_) const
	{
		assert((0 <= i_ && i_ < m_Rows_) && (0 <= j_ && j_ < m_Cols_));
		return m_Value_[i_][j_];
	}
	inline eTrace& GetTrace(size_t i_, size_t j_) const
	{
		assert((0 <= i_ && i_ < m_Rows_) && (0 <= j_ && j_ < m_Cols_));
		return m_Trace_[i_][j_];
	}

	inline T** GetValue() const
	{
		return m_Value_;
	}
	inline eTrace** GetTrace() const
	{
		return m_Trace_;
	}
	inline T**& GetValue()
	{
		return m_Value_;
	}
	inline eTrace**& GetTrace()
	{
		return m_Trace_;
	}

private:
	size_t m_Rows_, m_Cols_;
	T **m_Value_;
	eTrace **m_Trace_;
};

template <typename T>
class Solver
{
public:
	Solver(const Sequence& seq1_, const Sequence& seq2_, const BloSum62<T>& blosum62_) :
//	Solver(const Sequence& seq1_, const Sequence& seq2_) :
		m_Seq1_(seq1_), m_Seq2_(seq2_), m_BloSum62_(blosum62_),
//		m_Seq1_(seq1_), m_Seq2_(seq2_),
		m_Seq1Len_(seq1_.GetSequence().length()), m_Seq2Len_(seq2_.GetSequence().length()),
		m_MtxA_(m_Seq1Len_ + 1, m_Seq2Len_ + 1), m_MtxB_(m_Seq1Len_ + 1, m_Seq2Len_ + 1), m_MtxC_(m_Seq1Len_ + 1, m_Seq2Len_ + 1)
	{
		InitializeSpace();
	}

	inline Space<T>& GetSpaceA()
	{
		return m_MtxA_;
	}
	inline Space<T>& GetSpaceB()
	{
		return m_MtxB_;
	}
	inline Space<T>& GetSpaceC()
	{
		return m_MtxC_;
	}

	T Update()
	{
		// get these two sequences.
		std::string seq1 = m_Seq1_.GetSequence(), seq2 = m_Seq2_.GetSequence();
		// update matrix A, B and C.
		T valA = (T)0, valB = (T)0, valC = (T)0;
		for(size_t i = 1; i <= m_Seq1Len_; ++i)
		{
			for(size_t j = 1; j <= m_Seq2Len_; ++j)
			{
				/// update the matrix A
				T sigma = m_BloSum62_.GetValue(seq1[i - 1], seq2[j - 1]);
				valA = sigma + m_MtxA_.GetValue(i - 1, j - 1);
				valB = sigma + m_MtxB_.GetValue(i - 1, j - 1);
				valC = sigma + m_MtxC_.GetValue(i - 1, j - 1);
				if((valA >= valB) && (valA >= valC))			// valA is the max value.
				{
					m_MtxA_.SetValue(valA, i, j);
					m_MtxA_.SetTrace(eTrace::FromA, i, j);
				} 
				else if((valB >= valA) && (valB >= valC))		// valB is the max value.
				{
					m_MtxA_.SetValue(valB, i, j);
					m_MtxA_.SetTrace(eTrace::FromB, i, j);
				} else if((valC >= valA) && (valC >= valB))		// valC is the max value.
				{
					m_MtxA_.SetValue(valC, i, j);
					m_MtxA_.SetTrace(eTrace::FromC, i, j);
				}

				/// update the matrix B
				valA = m_MtxA_.GetValue(i, j - 1) - (g_Wg + g_Ws); 
				valB = m_MtxB_.GetValue(i, j - 1) - g_Ws; 
				valC = m_MtxC_.GetValue(i, j - 1) - (g_Wg + g_Ws);
				if((valA >= valB) && (valA >= valC))
				{
					m_MtxB_.SetValue(valA, i, j);
					m_MtxB_.SetTrace(eTrace::FromA, i, j);
				} else if((valB >= valA) && (valB >= valC))	
				{
					m_MtxB_.SetValue(valB, i, j);
					m_MtxB_.SetTrace(eTrace::FromB, i, j);
				} else if((valC >= valA) && (valC >= valB))
				{
					m_MtxB_.SetValue(valC, i, j);
					m_MtxB_.SetTrace(eTrace::FromC, i, j);
				}

				/// update the matrix C
				valA = m_MtxA_.GetValue(i - 1, j) - (g_Wg + g_Ws);
				valB = m_MtxB_.GetValue(i - 1, j) - (g_Wg + g_Ws);
				valC = m_MtxC_.GetValue(i - 1, j) - g_Ws;
				if((valA >= valB) && (valA >= valC))
				{
					m_MtxC_.SetValue(valA, i, j);
					m_MtxC_.SetTrace(eTrace::FromA, i, j);
				} else if((valB >= valA) && (valB >= valC))
				{
					m_MtxC_.SetValue(valB, i, j);
					m_MtxC_.SetTrace(eTrace::FromB, i, j);
				} else if((valC >= valA) && (valC >= valB))
				{
					m_MtxC_.SetValue(valC, i, j);
					m_MtxC_.SetTrace(eTrace::FromC, i, j);
				}
			}
		}

		// get the last score of each value matrix.
		valA = m_MtxA_.GetValue(m_Seq1Len_, m_Seq2Len_);
		valB = m_MtxB_.GetValue(m_Seq1Len_, m_Seq2Len_);
		valC = m_MtxB_.GetValue(m_Seq1Len_, m_Seq2Len_);
		// return the max score.
		return (valA >= valB) ? ((valA >= valC) ? valA : valC) : ((valB >= valC) ? valB : valC);
	}

	SeqPair Construct() const
	{
		std::string origSeq1 = m_Seq1_.GetSequence(), origSeq2 = m_Seq2_.GetSequence();
		std::string resSeq1 = "", resSeq2 = "";

		// initialize the current trace & some other parameters.
		T valA = (T)0, valB = (T)0, valC = (T)0;
		eTrace current = eTrace::None;
		size_t i = m_Seq1Len_, j = m_Seq2Len_;

		// get the value of each matrix.
		valA = m_MtxA_.GetValue(i, j); valB = m_MtxB_.GetValue(i, j); valC = m_MtxC_.GetValue(i, j);

		// get the start matrix.
		if((valA >= valB) && (valA >= valC))
			current = eTrace::FromA;
		else if((valB >= valA) && (valB >= valC))
			current = eTrace::FromB;
		else if((valC >= valA) && (valC >= valB))
			current = eTrace::FromC;

		// start trace back.
		while((i > 0) || (j > 0))
		{
			if(current == eTrace::FromA)
			{
				// get the trace.
				current = m_MtxA_.GetTrace(i, j);
				// trace back.
				resSeq1 = std::string(1, origSeq1[i - 1]) + resSeq1;	// there is a offset 1 between origSeq1 & (A or B or C).
				resSeq2 = std::string(1, origSeq2[j - 1]) + resSeq2;
				--i; --j;
			} 
			else if(current == eTrace::FromB)
			{
				// get the trace.
				current = m_MtxB_.GetTrace(i, j);
				// trace back.
				resSeq1 = "-" + resSeq1;	
				resSeq2 = std::string(1, origSeq2[j - 1]) + resSeq2;	// there is a offset 1 between origSeq1 & (A or B or C).
				--j;
			} 
			else if(current == eTrace::FromC)
			{
				// get the trace.
				current = m_MtxC_.GetTrace(i, j);
				// trace back.
				resSeq1 = std::string(1, origSeq1[i - 1]) + resSeq1;	// there is a offset 1 between origSeq1 & (A or B or C).
				resSeq2 = "-" + resSeq2;
				--i;
			}
			// debug
//			std::cout << "[" << i << ", " << j << "]" << std::endl;
		}

		return SeqPair(resSeq1, resSeq2);
	}

public:
	static const T g_Wg;
	static const T g_Ws;

private:
	void InitializeSpace()
	{
		T** ptrValue = nullptr;

		// initialize matrix A
		ptrValue = m_MtxA_.GetValue();

		ptrValue[0][0] = 0;							// A[0][0] = 0
		for(size_t k = 1; k <= m_Seq1Len_; ++k)		// A[1 ... end][0] = -inf
			ptrValue[k][0] = (T)NEGINF;
		for(size_t k = 1; k <= m_Seq2Len_; ++k)		// A[0][1 ... end] = -inf
			ptrValue[0][k] = (T)NEGINF;

		// initialize matrix B
		ptrValue = m_MtxB_.GetValue();

		for(size_t k = 0; k <= m_Seq1Len_; ++k)		// B[0 ... end][0] = -inf
			ptrValue[k][0] = (T)NEGINF;
		for(size_t k = 1; k <= m_Seq2Len_; ++k)		// B[0][1 ... end] = -(Wg + k*Ws)
			ptrValue[0][k] = (T)(-(g_Wg + ((T)k)*g_Ws));

		// initialize matrix C
		ptrValue = m_MtxC_.GetValue();

		for(size_t k = 1; k <= m_Seq1Len_; ++k)		// B[1 ... end][0] = -(Wg + k*Ws)
			ptrValue[k][0] = (T)(-(g_Wg + ((T)k)*g_Ws));
		for(size_t k = 0; k <= m_Seq2Len_; ++k)		// B[0 ... end][0] = -inf
			ptrValue[0][k] = (T)NEGINF;
	}

private:
	BloSum62<T> m_BloSum62_;
	Sequence m_Seq1_, m_Seq2_;
	size_t m_Seq1Len_, m_Seq2Len_;

	Space<T> m_MtxA_;
	Space<T> m_MtxB_;
	Space<T> m_MtxC_;
};

template <typename T>
const T Solver<T>::g_Wg = 10;

template <typename T>
const T Solver<T>::g_Ws = 2;

#endif	// __SOLVER_HXX__
