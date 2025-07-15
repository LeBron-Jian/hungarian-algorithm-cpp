///////////////////////////////////////////////////////////////////////////////
/*
 Hungarian.h: Header file for Class HungarianAlgorithm.

 This is a C++ wrapper with slight modification of a hungarian algorithm implementation by Markus Buehren.
 The original implementation is a few mex-functions for use in MATLAB, found here:
 http://www.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem

 Both this code and the orignal code are published under the BSD license.
 by Cong Ma, 2016

    HungarianAlogorithm类 是一个典型的，基于“算法步骤“的匈牙利算法实现。
    这个实现将整个算法过程分解为一系列独立的，可按顺序调用的步骤函数
*/

#ifndef HUNGARIAN_H
#define HUNGARIAN_H

#include <iostream>
#include <vector>

using namespace std;


class HungarianAlgorithm
{
public:
	HungarianAlgorithm();  // 构造函数，通常用于初始化类的成员变量
	~HungarianAlgorithm();  // 析构函数，用于清理在对象声明周期中分配的资源，例如动态分配的数组
	// 主函数入口：接受两个引用参数，DistMatrix：成本矩阵（或距离矩阵），Assignment：函数将最终的匹配结果写入这个vector，返回最终的最小总成本
	double Solve(vector <vector<double> >& DistMatrix, vector<int>& Assignment);

private:
	void assignmentoptimal(int *assignment, double *cost, double *distMatrix, int nOfRows, int nOfColumns);
	void buildassignmentvector(int *assignment, bool *starMatrix, int nOfRows, int nOfColumns);
	void computeassignmentcost(int *assignment, double *cost, double *distMatrix, int nOfRows);
	void step2a(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step2b(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step3(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
	void step4(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim, int row, int col);
	void step5(int *assignment, double *distMatrix, bool *starMatrix, bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns, bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
};


#endif
