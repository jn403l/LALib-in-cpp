// ２次元行列に関する基本的な計算機能
#ifndef QBMATRIX2_H
#define QBMATRIX2_H

#include <stdexcept>
/**
standard two dimensional array:
matrixArray = new double[3][3]
however, we can use matrixArray = new double[nElements]
where nElements = nRows * nCols;
and linearIndex = (row * numCols) + col;
 **/

template<class T>
class qbMatrix2 {
public:
	// define the various constructors
	qbMatrix2();
	qbMatrix2(int nRows, int nCols);
	qbMatrix2(int nRows, int nCols, const T *inputData);
	qbMatrix2(const qbMatrix2<T>& inputMatrix);
  qbMatrix2<T>& operator=(const qbMatrix2<T>& otherMatrix);

	// define the destructor
	~qbMatrix2();

	// configuration methods
	bool resize(int numRows, int numCols);

	// element access methods
	T GetElement(int row, int col);
	bool SetElement(int row, int col, T elementValue);
	int GetNumRows();
	int GetNumCols();

	// overload == operator
	bool operator== (const qbMatrix2<T>& rhs) const;

	// overload +, - and * operators (friends)
  /* (1) 行列と行列の足し算 */
	template<class U> friend qbMatrix2<U> operator+ (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
  /* (2) 行列とスカラーの足し算 */
	template<class U> friend qbMatrix2<U> operator+ (const U& lhs, const qbMatrix2<U>& rhs);
  /* (3) スカラーと行列の足し算 */
	template<class U> friend qbMatrix2<U> operator+ (const qbMatrix2<U>&lhs, const U& rhs);

	template<class U> friend qbMatrix2<U> operator- (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
	template<class U> friend qbMatrix2<U> operator- (const U& lhs, const qbMatrix2<U>& rhs);
	template<class U> friend qbMatrix2<U> operator- (const qbMatrix2<U>&lhs, const U& rhs);

	template<class U> friend qbMatrix2<U> operator* (const qbMatrix2<U>& lhs, const qbMatrix2<U>& rhs);
	template<class U> friend qbMatrix2<U> operator* (const U& lhs, const qbMatrix2<U>& rhs);
	template<class U> friend qbMatrix2<U> operator* (const qbMatrix2<U>&lhs, const U& rhs);

private:
	int Sub2Ind(int row, int col);

private:
  T *m_matrixData;
  int m_nRows, m_nCols, m_nElements;
};

/***************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTION
**************************************************/

// default constructor
// TODO: replace literal 0.0 with value-initialization T{}?
template <class T>
qbMatrix2<T>::qbMatrix2() {
  m_nRows = 1;
  m_nCols = 1;
  m_nElements = 1;
  m_matrixData = new T[m_nElements];
  m_matrixData[0] = 0.0;
}

// construct 零行列 (すべての要素は0)
template<class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols) {
  m_nRows = nRows;
  m_nCols = nCols;
  m_nElements = m_nRows * m_nCols;
  m_matrixData = new T[m_nElements];
  for (int i = 0; i < m_nElements; i++)
    m_matrixData[i] = 0.0;
}

// construct from const 配列
template<class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols, const T *inputData) {
  m_nRows = nRows;
  m_nCols = nCols;
  m_nElements = m_nRows * m_nCols;
  m_matrixData = new T[m_nElements];
  for (int i = 0; i < m_nElements; i++)
    m_matrixData[i] = inputData[i];
}

// copy constructor
template<class T>
qbMatrix2<T>::qbMatrix2(const qbMatrix2<T>& inputMatrix) {
  m_nRows = inputMatrix.m_nRows;
  m_nCols = inputMatrix.m_nCols;
  m_nElements = m_nRows * m_nCols;
  m_matrixData = new T[m_nElements];
  for (int i = 0; i < m_nElements; i++)
    m_matrixData[i] = inputMatrix.m_matrixData[i];
}

// copy assignment operator
template<class T>
qbMatrix2<T>& qbMatrix2<T>::operator=(const qbMatrix2<T>& otherMatrix) {
  if (this == &otherMatrix) return *this;
  m_nRows = otherMatrix.m_nRows;
  m_nCols = otherMatrix.m_nCols;
  m_nElements = m_nRows * m_nCols;
  if (this != &otherMatrix) {
    T* m_newMatrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
      m_newMatrixData[i] = otherMatrix.m_matrixData[i];
    delete[] m_matrixData;
    m_matrixData = m_newMatrixData;
  }
  return *this;
}

// destructor
template<class T>
qbMatrix2<T>::~qbMatrix2() {
  if (m_matrixData != nullptr)
    delete[] m_matrixData;
}


/***************************************************
CONFIGURATION FUNCTIONS
**************************************************/
template<class T>
bool qbMatrix2<T>::resize(int numRows, int numCols) {
  m_nRows = numRows;
  m_nCols = numCols;
  m_nElements = (m_nRows * m_nCols);
  delete[] m_matrixData;
  m_matrixData = new T[m_nElements];
  if (m_matrixData != nullptr) {
    for (int i = 0; i < m_nElements; i++)
      m_matrixData[i] = 0.0;

    return true;
  } else {
    return false;
  }
}

/***************************************************
ELEMENT FUNCTIONS
 **************************************************/
// TODO: throw on bounds error?
template<class T>
T qbMatrix2<T>::GetElement(int row, int col) {
  int linearIndex = Sub2Ind(row, col);
  if (linearIndex >= 0)
    return m_matrixData[linearIndex];
  else
    return 0.0;
}

template<class T>
bool qbMatrix2<T>::SetElement(int row, int col, T elementValue) {
  int linearIndex = Sub2Ind(row, col);
  if (linearIndex >= 0) {
    m_matrixData[linearIndex] = elementValue;
    return true;
  } else  {
    return false;
  }
}

template<class T>
int qbMatrix2<T>::GetNumRows() {
  return m_nRows;
}

template<class T>
int qbMatrix2<T>::GetNumCols() {
  return m_nCols;
}

/***************************************************
OVERLOADED OPERATOR FUNCTIONS
**************************************************/

/***************************************************
THE + OPERATOR
**************************************************/
// 行列 + 行列
template<class T>
qbMatrix2<T> operator+ (const qbMatrix2<T>& lhs, const qbMatrix2<T>& rhs) {
  if ((lhs.m_nCols != rhs.m_nCols) && (lhs.m_nRows != rhs.m_nRows))
    throw std::invalid_argument("operator+: shape mismatch");
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; i++)
    tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

// スカラー + 行列
template<class T>
qbMatrix2<T> operator+ (const T& lhs, const qbMatrix2<T>& rhs) {
  int numRows = rhs.m_nRows;
  int numCols = rhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; ++i)
    tempResult[i] = lhs + rhs.m_matrixData[i];

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

// 行列 + スカラー
template<class T>
qbMatrix2<T> operator+ (const qbMatrix2<T>& lhs, const T& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; ++i)
    tempResult[i] = lhs.m_matrixData[i] + rhs;

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

/***************************************************
THE - OPERATOR
**************************************************/
// 行列 - 行列
template<class T>
qbMatrix2<T> operator- (const qbMatrix2<T>& lhs, const qbMatrix2<T>& rhs) {
  if ((lhs.m_nCols != rhs.m_nCols) && (lhs.m_nRows != rhs.m_nRows))
    throw std::invalid_argument("operator-: shape mismatch");
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; i++)
    tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

// スカラー - 行列
template<class T>
qbMatrix2<T> operator- (const T& lhs, const qbMatrix2<T>& rhs) {
  int numRows = rhs.m_nRows;
  int numCols = rhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; ++i)
    tempResult[i] = lhs - rhs.m_matrixData[i];

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

// 行列 - スカラー
template<class T>
qbMatrix2<T> operator- (const qbMatrix2<T>& lhs, const T& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; ++i)
    tempResult[i] = lhs.m_matrixData[i] - rhs;

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

/***************************************************
THE * OPERATOR
**************************************************/
// スカラー * 行列
template<class T>
qbMatrix2<T> operator* (const T& lhs, const qbMatrix2<T>& rhs) {
  int numRows = rhs.m_nRows;
  int numCols = rhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; ++i)
    tempResult[i] = lhs * rhs.m_matrixData[i];

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

// 行列 - スカラー
template<class T>
qbMatrix2<T> operator* (const qbMatrix2<T>& lhs, const T& rhs) {
  int numRows = lhs.m_nRows;
  int numCols = lhs.m_nCols;
  int numElements = numRows * numCols;
  T *tempResult = new T[numElements];
  for (int i = 0; i < numElements; ++i)
    tempResult[i] = lhs.m_matrixData[i] * rhs;

  qbMatrix2<T> result(numRows, numCols, tempResult);
  delete[] tempResult;
  return result;
}

// 行列 * 行列
template<class T>
qbMatrix2<T> operator* (const qbMatrix2<T>& lhs, const qbMatrix2<T>& rhs) {
  int r_numRows = rhs.m_nRows;
  int r_numCols = rhs.m_nCols;
  int l_numRows = lhs.m_nRows;
  int l_numCols = lhs.m_nCols;

  if (l_numCols == r_numRows) {
    // this is the standard matrix multiplication condition.
    // the output will be the same size as the RHS
    T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];

    // loop through each row of the LHS
    for (int  lhsRow = 0; lhsRow < l_numRows; lhsRow++) {
      // loop through each column on the RHS
      for (int rhsCol = 0; rhsCol < r_numCols; rhsCol++) {
        T elementResult = 0.0;
        // loop through each element of this LHS row
        for (int lhsCol = 0; lhsCol < l_numCols; lhsCol++) {
          // compute the LHS linear index
          int lhsLinearIndex = (lhsRow + l_numCols) + lhsCol;

          // compute the RHS linear index (based on LHS col)
          // rhs row number equal to lhs colum number
          int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;

          // perform the calculation on these elements
          elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
        }

        // store the result
        int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
        tempResult[resultLinearIndex] = elementResult;
      }
    }
    qbMatrix2<T> result(l_numRows, r_numCols, tempResult);
    delete[] tempResult;
    return result;
  } else {
    qbMatrix2<T> result(1, 1);
    return result;
  }
}


/***************************************************
THE == OPERATOR
**************************************************/
template<class T>
bool qbMatrix2<T>::operator== (const qbMatrix2<T>& rhs) const {
  // TODO: check if && should be ||
  // check if the matrices are the same size, if not return false
  if ((this->m_nRows != rhs.m_nRows) && (this->m_nCols != rhs.m_nCols))
    return false;

  // check if the elements are equal
  bool flag = true;
  for (int i = 0; i < this->m_nElements; ++i) {
    if (this->m_matrixData[i] != rhs.m_matrixData[i])
      flag = false;
  }
  return flag;
}

/***************************************************
PRIVATE FUNCTIONS
**************************************************/
template<class T>
int qbMatrix2<T>::Sub2Ind(int row, int col) {
  if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
    return (row * m_nCols) + col;
  else
    return -1;
}

#endif
