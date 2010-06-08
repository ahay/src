#ifndef _VECMATOP_HPP_
#define _VECMATOP_HPP_

#include "nummat.hh"

using std::vector;
//typedef long int int;

//--------------------------------------------------
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C);
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C);

int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y);
int dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y);

//--------------------------------------------------
int zgemm(cpx alpha, const CpxNumMat& A, const CpxNumMat& B, cpx beta, CpxNumMat& C);
int zgemm(int m, int n, int k, cpx alpha, cpx* A, cpx* B, cpx beta, cpx* C);

int zgemv(cpx alpha, const CpxNumMat& A, const CpxNumVec& X, cpx beta, CpxNumVec& Y);
int zgemv(int m, int n, cpx alpha, cpx* A, cpx* X, cpx beta, cpx* Y);

//--------------------------------------------------
int dgmres(int (*A)(const DblNumVec&, DblNumVec&), const DblNumVec& b, const DblNumVec& x0,
	   int restart, double tol, int maxit, int print,
	   DblNumVec& x, int& flag, double& relres, int& iter, vector<double>& resvec);

int zgmres(int (*A)(const CpxNumVec&, CpxNumVec&), const CpxNumVec& b, const CpxNumVec& x0,
	   int restart, double tol, int maxit, int print,
	   CpxNumVec& x, int& flag, double& relres, int& iter, vector<double>& resvec);

//--------------------------------------------------
int pinv(const DblNumMat& M, double eps, DblNumMat& R);

//--------------------------------------------------
int lowrank(int m, int n, int (*sample)(vector<int>&, vector<int>&, DblNumMat&), double eps, int npk,
	    vector<int>& cidx, vector<int>& ridx, DblNumMat& mid);

#endif

