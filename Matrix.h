#pragma once
#include <vector>
#include <complex>
#include <ostream>
#include <stdexcept>
#include <cmath>
#include <utility>
#include <cstddef>

/* 템플릿을 사용하는 경우 선언과 정의 모두 하나의 파일에 정리하는게 컴파일러의 입장에서 낫다고 한다.
 아니면 에러가 나는 경우가 생긴다.
 이걸 해결하기 위해서 Matrix.cpp를 Matrix.tpp로 바꾸고 include 해줘야만 했다.
 그리고 그냥 Matrix.cpp는 결국 Matrix_impl.hpp로 바꼈다.*/
// 결국 통합하여 정리함

// 혹시나 다른 헤더 등에서 함수명이 겹칠 일이 생길까봐 익명 namespace를 활용함
namespace{
template <typename U>
U conjugate_if_needed(const U& x) {
    return x;
}

template <typename U>
std::complex<U> conjugate_if_needed(const std::complex<U>& x) {
    return std::conj(x);
}
}

template <typename T>
class Matrix{
private:
    std::size_t rows_, cols_;
    std::vector<std::vector<T>> contents_;
    
public:
    Matrix() : rows_{0}, cols_{0}, contents_{}{}
    Matrix(std::size_t ro, std::size_t co) : rows_{ro}, cols_{co}, contents_{ro, std::vector<T>(co)}{}
    Matrix(std::vector<std::vector<T>> vec);
    std::size_t rows() const{return rows_;}
    std::size_t cols() const{return cols_;}
    /*
    void set_row(std::size_t a){rows_ = a;}
    void set_col(std::size_t b){cols_ = b;}
    */
    T& operator()(std::size_t i, std::size_t j){return contents_[i][j];}
    const T& operator()(std::size_t i, std::size_t j) const{return contents_[i][j];}
    void setdata(std::size_t row, std::size_t col, T data){contents_[row][col] = data;}
    
    Matrix& operator+=(const Matrix<T>& other);
    Matrix& operator-=(const Matrix& other);
    // 대입연산자는 멤버함수로만 정의할 수 있으며, 클래스 외부함수로 정의할 수는 없다. => C++문법 규칙이라고..
    Matrix& operator=(const Matrix& other); // 직접 정의하지 않아도 컴파일러가 대개 자동으로 정의해주긴 한다고...
    Matrix& operator*=(const Matrix& other);
    Matrix& operator*=(const T& scalar);
    Matrix& operator/=(const T& scalar);
    
    bool operator==(const Matrix& other); // 구현해야 함
    
    Matrix transpose() const;
    Matrix dagger() const;
    Matrix adjoint() const{return dagger();}
};

template <typename T>
class SquareMatrix : public Matrix<T>{
private:
    
public:
    explicit SquareMatrix(size_t n) : Matrix<T>(n,n) {}
    SquareMatrix(std::vector<std::vector<T>> vec);
    SquareMatrix() : Matrix<T>(){}
    size_t order() const{return this->rows();}
    T trace() const;
    SquareMatrix inverse() const;
    
    T determinant() const;
    T det() const{return determinant();}
    
    static SquareMatrix identity(std::size_t n){
        SquareMatrix ansmat(n);
        for(std::size_t i=0; i < n; ++i){
            ansmat(i, i) = T{1};
        }
        return ansmat;
    }
    
    // 이하 구현해야 함
    SquareMatrix& operator+=(const SquareMatrix<T>& other);
    SquareMatrix& operator-=(const SquareMatrix<T>& other);
    SquareMatrix& operator*=(const SquareMatrix<T>& other);
    SquareMatrix& operator*=(const T& scalar);
    SquareMatrix& operator/=(const T& scalar);
    
    SquareMatrix transpose() const;
    SquareMatrix dagger() const;
    SquareMatrix adjoint() const{return dagger();}
    
    // 가장 큰 eigenvalue를 구하는 함수
    T dominant_eigenvalue(int max_iter = 1000, double EPS = 1.0e-10) const;
};

// 추가 연산 설계도

template <typename T>
Matrix<T> operator+(Matrix<T>, const Matrix<T>&);
template <typename t>
Matrix<t> operator-(Matrix<t>, const Matrix<t>&);
template <typename T>
Matrix<T> operator*(const Matrix<T>&, const Matrix<T>&);
template <typename T>
Matrix<T> operator*(const T&, const Matrix<T>&);
template <typename T>
Matrix<T> operator*(const Matrix<T>&, const T&);
template <typename T>
Matrix<T> operator/(Matrix<T>, const T&);

template <typename T>
SquareMatrix<T> operator+(SquareMatrix<T>&, const SquareMatrix<T>&);
template <typename t>
SquareMatrix<t> operator-(SquareMatrix<t>, const SquareMatrix<t>&);
template <typename T>
SquareMatrix<T> operator*(const SquareMatrix<T>&, const SquareMatrix<T>&);
template <typename T>
SquareMatrix<T> operator*(const T&, SquareMatrix<T>&);
template <typename T>
SquareMatrix<T> operator*(SquareMatrix<T>&, const T&);
template <typename T>
SquareMatrix<T> operator/(SquareMatrix<T>, const T&);

template <typename T>
std::ostream& operator<<(std::ostream&, const Matrix<T>&);


// 여기부터 정의들

// Matrix member

template <typename T>
Matrix<T>::Matrix(std::vector<std::vector<T>> vec) : rows_{vec.size()}, cols_{vec.empty()?0:vec[0].size()}, contents_{std::move(vec)}{
    for(const auto& row : contents_){
        if(row.size() != cols_)
            throw std::invalid_argument("Rows have differnt sizes.");
    }
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other){
    if(rows_!=other.rows()||cols_!=other.cols())
        throw std::invalid_argument("Two Matrix have different sizes.");
    for(std::size_t i=0; i<rows_; ++i){
        for(std::size_t j = 0; j < cols_; ++j){
            contents_[i][j]+=other(i,j);
        }
    }
    return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix& other){
    if(rows_!=other.rows()||cols_!=other.cols())
        throw std::invalid_argument("Two Matrix have different sizes.");
    for(std::size_t i=0; i<rows_; ++i){
        for(std::size_t j = 0; j < cols_; ++j){
            contents_[i][j]-=other(i,j);
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& other){
    rows_ = other.rows_; // 멤버함수라 굳이 rows()를 호출 안 하고 rows_를 바로 불러와도 됨.
    cols_ = other.cols(); // 그래도 other.cols()로 할 수 있긴 하고, 비용이 비싸지 않으니 그냥 이렇게도 해봄
    contents_ = other.contents_; // 멤버함수이므로 같은 클래스 객체인 other의 other.contents_가 private이라도 쓸 수 있음
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix& other){
    return (*this) = (*this)*other;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T& scalar){
    for(std::size_t i = 0; i < rows_; ++i){
        for(std::size_t j = 0; j < cols_; ++j){
            (*this)(i,j) *= scalar;
        }
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator/=(const T& scalar){
    for(std::size_t i = 0; i < rows_; ++i){
        for(std::size_t j = 0; j < cols_; ++j){
            (*this)(i,j) /= scalar;
        }
    }
    return *this;
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& other){
    if(cols_ != other.cols())
        return false;
    if(rows_ != other.rows())
        return false;
    for(std::size_t i = 0; i < rows_; ++i){
        for(std::size_t j = 0; j < cols_; ++j){
            if((*this)(i,j) != other(i,j))
                return false;
        }
    }
    return true;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const{
    Matrix ansmat(cols_, rows_);
    for(std::size_t i = 0; i < rows_; ++i){
        for(std::size_t j = 0; j < cols_; ++j)
            ansmat(j,i) = (*this)(i,j);
    }
    return ansmat;
}

template <typename T>
Matrix<T> Matrix<T>::dagger() const {
    Matrix<T> ansmat(cols_, rows_);

    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            ansmat(j, i) = conjugate_if_needed((*this)(i, j));
        }
    }

    return ansmat;
}

// SquareMatrix member

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator+=(const SquareMatrix<T>& other){
    Matrix<T>::operator+=(other);
    return *this;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator-=(const SquareMatrix<T>& other){
    Matrix<T>::operator-=(other);
    return *this;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const SquareMatrix<T>& other){
    *this = *this * other;
    return *this;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const T& scalar){
    Matrix<T>::operator*=(scalar);
    return *this;
}

template <typename T>
SquareMatrix<T>& SquareMatrix<T>::operator/=(const T& scalar){
    Matrix<T>::operator/=(scalar);
    return *this;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::transpose() const{
    const std::size_t n = this->order();
    SquareMatrix<T> ansmat(n);

    for(std::size_t i = 0; i < n; ++i){
        for(std::size_t j = 0; j < n; ++j){
            ansmat(j, i) = (*this)(i, j);
        }
    }

    return ansmat;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::dagger() const{
    const std::size_t n = this->order();
    SquareMatrix<T> ansmat(n);

    for(std::size_t i = 0; i < n; ++i){
        for(std::size_t j = 0; j < n; ++j){
            ansmat(j, i) = conjugate_if_needed((*this)(i, j));
        }
    }

    return ansmat;
}

template <typename T>
T SquareMatrix<T>::trace() const{
    T ans{};
    for(std::size_t i=0; i<this->order(); ++i){ // this->order()이 반복마다 호출; order()이 복잡해지면 그냥 미리 값을 저장해두고 쓸 것!
        ans += (*this)(i,i);
    }
    return ans;
}

template <typename T>
SquareMatrix<T>::SquareMatrix(std::vector<std::vector<T>> vec) : Matrix<T>(vec){
    if(this->rows() != this->cols())
        throw std::invalid_argument("Rows and cols are different");
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::inverse() const{
    const std::size_t n = this->order();

    SquareMatrix<T> left(n);
    SquareMatrix<T> right = SquareMatrix<T>::identity(n);

    for(std::size_t i = 0; i < n; ++i){
        for(std::size_t j = 0; j < n; ++j){
            left(i,j) = (*this)(i,j);
        }
    }

    const double EPS = 1e-12;

    for(std::size_t col = 0; col < n; ++col){
        std::size_t pivot = col;
        double max_abs = std::abs(left(col, col));

        for(std::size_t row = col + 1; row < n; ++row){
            double now_abs = std::abs(left(row, col));

            if(now_abs > max_abs){
                max_abs = now_abs;
                pivot = row;
            }
        }

        if(max_abs < EPS){
            throw std::runtime_error("Matrix is singular.");
        }

        if(pivot != col){
            for(std::size_t j = 0; j < n; ++j){
                std::swap(left(col, j), left(pivot, j));
                std::swap(right(col, j), right(pivot, j));
            }
        }

        T div = left(col, col);

        for(std::size_t j = 0; j < n; ++j){
            left(col, j) /= div;
            right(col, j) /= div;
        }

        for(std::size_t row = 0; row < n; ++row){
            if(row == col) continue;

            T factor = left(row, col);

            for(std::size_t j = 0; j < n; ++j){
                left(row, j) -= factor * left(col, j);
                right(row, j) -= factor * right(col, j);
            }
        }
    }

    return right;
}

template <typename T>
T SquareMatrix<T>::determinant() const {
    const std::size_t n = this->order();

    SquareMatrix<T> mat(n);

    for(std::size_t i = 0; i < n; ++i){
        for(std::size_t j = 0; j < n; ++j){
            mat(i,j) = (*this)(i,j);
        }
    }

    const double EPS = 1e-12;
    int sign = 1;

    for(std::size_t col = 0; col < n; ++col){
        std::size_t pivot = col;
        double max_abs = std::abs(mat(col, col));

        for(std::size_t row = col + 1; row < n; ++row){
            double now_abs = std::abs(mat(row, col));

            if(now_abs > max_abs){
                max_abs = now_abs;
                pivot = row;
            }
        }

        if(max_abs < EPS){
            return T{};
        }

        if(pivot != col){
            for(std::size_t j = 0; j < n; ++j){
                std::swap(mat(col, j), mat(pivot, j));
            }
            sign *= -1;
        }

        for(std::size_t row = col + 1; row < n; ++row){
            T factor = mat(row, col) / mat(col, col);

            for(std::size_t j = col; j < n; ++j){
                mat(row, j) -= factor * mat(col, j);
            }
        }
    }

    T ans = T{sign};

    for(std::size_t i = 0; i < n; ++i){
        ans *= mat(i,i);
    }

    return ans;
}

template <typename T>
T SquareMatrix<T>::dominant_eigenvalue(int max_iter, double EPS) const{
    const std::size_t n = this->order();

    if(n == 0){
        throw std::invalid_argument("Matrix is empty.");
    }

    std::vector<T> v(n, T{1});
    T lambda_old{};

    for(int iter = 0; iter < max_iter; ++iter){
        std::vector<T> w(n, T{});

        for(std::size_t i = 0; i < n; ++i){
            for(std::size_t j = 0; j < n; ++j){
                w[i] += (*this)(i, j) * v[j];
            }
        }

        double norm = 0.0;
        for(std::size_t i = 0; i < n; ++i){
            double abs_val = std::abs(w[i]);
            norm += abs_val * abs_val;
        }

        norm = std::sqrt(norm);

        if(norm < EPS){
            throw std::runtime_error("Cannot find dominant eigenvalue.");
        }

        for(std::size_t i = 0; i < n; ++i){
            v[i] = w[i] / norm;
        }

        std::vector<T> Av(n, T{});
        for(std::size_t i = 0; i < n; ++i){
            for(std::size_t j = 0; j < n; ++j){
                Av[i] += (*this)(i, j) * v[j];
            }
        }

        T numerator{};
        T denominator{};

        for(std::size_t i = 0; i < n; ++i){
            numerator += conjugate_if_needed(v[i]) * Av[i];
            denominator += conjugate_if_needed(v[i]) * v[i];
        }

        T lambda = numerator / denominator;

        if(std::abs(lambda - lambda_old) < EPS){
            return lambda;
        }

        lambda_old = lambda;
    }

    return lambda_old;
}

// 전역함수

template <typename T>
Matrix<T> operator+(Matrix<T> l, const Matrix<T>& r){
    return l+=r;
}

template <typename T>
Matrix<T> operator-(Matrix<T> l, const Matrix<T>& r){
    return l-=r;
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& l, const Matrix<T>& r){
    if(l.cols() != r.rows())
        throw std::invalid_argument("Can't multiply cause of sizes.");
    Matrix<T> ansmat{l.rows(), r.cols()};
    for(std::size_t i = 0; i < l.rows(); ++i){
        for(std::size_t j = 0; j < r.cols(); ++j){
            for(std::size_t k = 0; k<l.cols(); ++k){
                ansmat(i, j)+=l(i,k)*r(k,j);
            }
        }
    }
    
    return ansmat;
}

template <typename T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& mat){
    Matrix<T> ansmat = mat;
    return ansmat*=scalar;
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& mat, const T& scalar){
    return scalar*mat;
}

template <typename T>
Matrix<T> operator/(Matrix<T> mat, const T& scalar){
    return mat/=scalar;
}


template <typename T>
SquareMatrix<T> operator+(SquareMatrix<T> l, const SquareMatrix<T>& r){
    return l += r;
}

template <typename T>
SquareMatrix<T> operator-(SquareMatrix<T> l, const SquareMatrix<T>& r){
    return l -= r;
}

template <typename T>
SquareMatrix<T> operator*(const SquareMatrix<T>& l, const SquareMatrix<T>& r){
    if(l.order() != r.order()){
        throw std::invalid_argument("Two SquareMatrix have different orders.");
    }

    const std::size_t n = l.order();
    SquareMatrix<T> ansmat(n);

    for(std::size_t i = 0; i < n; ++i){
        for(std::size_t j = 0; j < n; ++j){
            for(std::size_t k = 0; k < n; ++k){
                ansmat(i, j) += l(i, k) * r(k, j);
            }
        }
    }

    return ansmat;
}

template <typename T>
SquareMatrix<T> operator*(const T& scalar, SquareMatrix<T> mat){
    return mat *= scalar;
}

template <typename T>
SquareMatrix<T> operator*(SquareMatrix<T> mat, const T& scalar){
    return mat *= scalar;
}

template <typename T>
SquareMatrix<T> operator/(SquareMatrix<T> mat, const T& scalar){
    return mat /= scalar;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat){
    if(mat.cols() == 0 || mat.rows() == 0)
        os << 0;
    
    for(std::size_t i = 0; i < mat.rows(); ++i){
        for(std::size_t j = 0; j < mat.cols(); ++j){
            os << mat(i,j) << ' ';
        }
        os << '\n';
    }
    return os;
}
