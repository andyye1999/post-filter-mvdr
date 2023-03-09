/* Created on 2016-08-29
 * Author: Zhang Binbin
 * About: Matrix implemention in cc 
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

// Here we just implement serveral simple matrix function

class Matrix {
public:
    Matrix(int r, int c): num_row_(r), num_col_(c) {
        data_ = new float[r * c](); 
    }

    Matrix(): num_row_(0), num_col_(0), data_(NULL) {}

    ~Matrix() {
        if (data_ != NULL) delete [] data_;
    }

    int NumRow() const { return num_row_; }
    int NumCol() const { return num_col_; }

    const float *Data() const {
        return data_;
    }
    
    float operator () (int r, int c) const {
        assert(r < num_row_);
        assert(c < num_col_);
        return *(data_ + r * num_col_ + c);
    }

    float & operator () (int r, int c) {
        assert(r < num_row_);
        assert(c < num_col_);
        return *(data_ + r * num_col_ + c);
    }

    void Resize(int r, int c) {
        if (r != num_row_ || c != num_col_) {
            if (data_ != NULL) delete [] data_;
            num_row_ = r;
            num_col_ = c;
            data_ = new float[r * c]();
        }
        else {
            memset(data_, 0, r * c * sizeof(float));
        }
    }

    void Copy(const Matrix &mat) {
        this->Resize(mat.NumRow(), mat.NumCol());
        memcpy(data_, mat.Data(), num_row_ * num_col_ * sizeof(float));
    }

    void Zeros() {
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                (*this)(i, j) = 0;
            }
        }
    }

    // *this = beta * (*this) + alpha * (mat)
    void Add(const Matrix &mat, float alpha = 1.0, float beta = 1.0) {
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                (*this)(i, j) = alpha * mat(i, j) + beta * (*this)(i, j); 
            }
        }
    }
    
    // Sum all element
    float Sum() {
        float sum = 0.0;
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                sum += (*this)(i, j);
            }
        }
        return sum;
    }

    // *this = mat1 .* mat2
    void MulElement(const Matrix &mat1, const Matrix &mat2) {
        assert(mat1.NumRow() == mat2.NumRow());
        assert(mat1.NumCol() == mat2.NumCol());
        this->Resize(mat1.NumRow(), mat1.NumCol());
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                (*this)(i, j) = mat1(i, j) * mat2(i, j);
            }
        }
    }

    // *this = mat1 * mat2
    void Mul(const Matrix &mat1, const Matrix &mat2) {
        assert(mat1.NumCol() == mat2.NumRow());
        this->Resize(mat1.NumRow(), mat2.NumCol());
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                for (int k = 0; k < mat1.NumCol(); k++) {
                    (*this)(i, j) += mat1(i, k) * mat2(k, j);        
                }
            }
        }
    }
    /*
    这段代码实现了一个矩阵求逆的函数，其实现方法是高斯-约旦消元法（Gaussian-Jordan elimination），
    基本思路是通过初等矩阵的左乘，将原矩阵变为一个单位矩阵，这个过程中左乘的初等矩阵相当于把原矩阵的每一行都变为对角线元素为 $1$，其余元素为 $0$ 的矩阵。
    具体实现步骤：
    首先判断输入矩阵是否是方阵，如果不是则无法求逆。
    创建一个和输入矩阵大小一致的矩阵，用于存储求逆结果。同时将这个矩阵初始化为单位矩阵，即对角线元素为 $1$，其余元素为 $0$。
    为了避免修改原始矩阵，将输入矩阵拷贝到临时矩阵中。
    对于每一行 $i$，选取对角线元素的绝对值最大的列 $c$，然后把第 $c$ 列所在的行变为对角线元素为 $1$，其余元素为 $0$ 的矩阵，
    同时在输出矩阵中对应的行做相应的变换，这个过程通过消元和代入得到，详见代码。
    如果有一行 $i$ 的对角线元素的绝对值为 $0$，则说明输入矩阵是奇异矩阵，无法求逆，返回 false。
    重复上述过程直到处理完所有行。
    返回 true，表示求逆成功，输出矩阵即为原矩阵的逆。
    */
    // *this = inv(mat);
    // if singular return false
    bool Inverse(const Matrix &mat) {
        assert(mat.NumRow() == mat.NumCol());
        Resize(mat.NumRow(), mat.NumCol());
        Matrix tmp_mat;
        tmp_mat.Copy(mat);
        for (int i = 0; i < num_row_; i++) {
            (*this)(i, i) = 1.0;
        }
        int *table = (int *)calloc(sizeof(int), num_row_);
        for (int i = 0; i < num_row_; i++) {
            // select max diag element
            int c = -1;
            float max = 0.0;
            for (int j = 0; j < num_row_; j++) {
                if (table[j] == 0 && fabs(tmp_mat(i, j)) > max) {
                    c = j;
                    max = fabs(tmp_mat(i, j));
                }
            }
            if (max == 0.0) {
                printf("singular matrix\n");
                return false;
            }
            
            table[c] = 1;
            float den = tmp_mat(c, c);
            for (int j = 0; j < num_row_; j ++) {
                if (j != c) { // other lines, not selected
                    float ratio = tmp_mat(j, c) / den;
                    for (int k = 0; k < num_col_; k++) {
                        tmp_mat(j, k) -= ratio * tmp_mat(c, k);
                        (*this)(j, k) -= ratio * (*this)(c, k);
                    }
                }
            }
            // c row, just scale
            for (int k = 0; k < num_row_; k++) {
                tmp_mat(c, k) /= den;
                (*this)(c, k) /= den;
            }
        }

        free(table);
        return true;
    }

    void Print() {
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                printf("%f\t", (*this)(i, j));
            }
            printf("\n");
        }
    }

private:
    int num_row_, num_col_;
    float *data_;
    // Disallow assign and copy
    Matrix(const Matrix &mat);
    Matrix & operator = (const Matrix &mat);
};


struct Complex {
    Complex(float r = 0, float i = 0): real(r), img(i) {}
    float real, img;
};

// Complex Matrix
class ComplexMatrix {
public:
    ComplexMatrix(int r, int c): num_row_(r), num_col_(c) {
        data_ = new Complex[r * c]();
    }

    ComplexMatrix(): num_row_(0), num_col_(0), data_(0) {}

    ~ComplexMatrix() {
        if (data_ != NULL) delete [] data_;
    }

    int NumRow() const { return num_row_; }
    int NumCol() const { return num_col_; }

    const Complex *Data() const {
        return data_;
    }
    
    Complex operator () (int r, int c) const {
        return data_[r * num_col_ + c];
    }

    Complex & operator () (int r, int c) {
        return data_[r * num_col_ + c];
    }

    void Resize(int r, int c) {
        if (r != num_row_ || c != num_col_) {
            if (data_ != NULL) delete [] data_;
            num_row_ = r;
            num_col_ = c;
            data_ = new Complex[r * c]();
        } else {
            memset(data_, 0, r * c * sizeof(Complex));
        }
    }
    
    void Copy(const ComplexMatrix &mat) {
        this->Resize(mat.NumRow(), mat.NumCol());
        memcpy(data_, mat.Data(), num_row_ * num_col_ * sizeof(Complex));
    }

    // *this = conj(mat);
    void Conj(const ComplexMatrix &mat) {
        this->Resize(mat.NumRow(), mat.NumCol());
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                (*this)(i, j).real = mat(i, j).real;
                (*this)(i, j).img = -mat(i, j).img;
            }
        }
    }

    // *this = beta * (*this) + alpha * (mat)
    void Add(const ComplexMatrix &mat, float alpha = 1.0, float beta = 1.0) {
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                ((*this)(i, j)).real = alpha * mat(i, j).real + 
                                       beta * ((*this)(i, j)).real; 
                ((*this)(i, j)).img = alpha * mat(i, j).img + 
                                       beta * ((*this)(i, j)).img; 
            }
        }
    }

    // *this = mat1 * mat2
    void Mul(const ComplexMatrix &mat1, const ComplexMatrix &mat2) {
        assert(mat1.NumCol() == mat2.NumRow());
        this->Resize(mat1.NumRow(), mat2.NumCol());
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                for (int k = 0; k < mat1.NumCol(); k++) {
                    ((*this)(i, j)).real += mat1(i, k).real * mat2(k, j).real - 
                                            mat1(i, k).img * mat2(k, j).img;
                    ((*this)(i, j)).img += mat1(i, k).real * mat2(k, j).img + 
                                            mat1(i, k).img * mat2(k, j).real;
                }
            }
        }
    }

    // *this = inv(mat);
    // if singular return false
    // C = A + iB
    // inv(C) = inv(A+B*inv(A)*B) - i inv(A)*B*inv(A+B*inv(A)*B)
    bool Inverse(const ComplexMatrix &mat) {
        assert(mat.NumRow() == mat.NumCol());
        this->Resize(mat.NumRow(), mat.NumCol());
        // assign real_mat(A), img_mat(B) C = A + iB
        Matrix a_mat(num_row_, num_col_), b_mat(num_row_, num_col_);
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                a_mat(i, j) = (mat(i,j)).real;
                b_mat(i, j) = (mat(i,j)).img;
            }
        }
        
        bool flag;
        Matrix inv_a_mat;
        flag = inv_a_mat.Inverse(a_mat);
        if (!flag) return false;

        Matrix mat1, mat2, real_mat, img_mat;
        mat1.Mul(b_mat, inv_a_mat);
        mat2.Mul(mat1, b_mat);
        mat2.Add(a_mat);
        flag = real_mat.Inverse(mat2);
        if (!flag) return false;

        mat1.Mul(inv_a_mat, b_mat);
        img_mat.Mul(mat1, real_mat);
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                ((*this)(i, j)).real = real_mat(i, j);
                ((*this)(i, j)).img = -img_mat(i, j);
            }
        }
        return true;
    }

    void Print() {
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                printf("%f+%fi\t", ((*this)(i, j)).real, ((*this)(i, j)).img);
            }
            printf("\n");
        }
    }

    // Trace sum of diag elem
    float Trace() {
        assert(num_row_ == num_col_);
        float sum = 0.0;
        for (int i = 0; i < num_row_; i++) {
            sum += (*this)(i, i).real;
        }
        return sum;
    }

    // *this = alpha * (*this)
    void Scale(float alpha) {
        for (int i = 0; i < num_row_; i++) {
            for (int j = 0; j < num_col_; j++) {
                (*this)(i, j).real *= alpha; 
                (*this)(i, j).img *= alpha; 
            }
        }
    }

    void ApplyDiagCeil(float val) {
        assert(num_row_ == num_col_);
        for (int i = 0; i < num_row_; i++) {
            (*this)(i, i).real += val;
            //if (fabs((*this)(i, i).real < val)) {
            //    (*this)(i, i).real = val;
            //}
        }
    }

private:
    int num_row_, num_col_;
    Complex *data_;
    // Disallow assign and copy
    ComplexMatrix(const ComplexMatrix &mat);
    ComplexMatrix & operator = (const ComplexMatrix &mat);
};

#endif
