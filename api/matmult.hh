#ifndef _matmult_hh
#define _matmult_hh

// Matrix Multiplication
class MatMult { 
public:
    // constructor
    MatMult(int nx, int ny, float **matrix);
    // linear operators
    void forward(bool add, const float *x, float *y);
    void adjoint(bool add, float *x, const float *y);
private:
    float **matrix_;
    int nx_, ny_;
};

#endif
