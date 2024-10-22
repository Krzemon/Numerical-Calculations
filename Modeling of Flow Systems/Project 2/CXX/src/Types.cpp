#include "Types.hpp"
#include <algorithm>
#include <numeric>

double VarVector::norm() const
{
    return std::sqrt(std::inner_product(mData.begin(), mData.end(), mData.begin(), 0.0));
}

VarVector VarVector::operator+(const VarVector &rhs) const
{
    VarVector result(mData.size());
    std::transform(mData.begin(), mData.end(), rhs.mData.begin(), result.mData.begin(), std::plus<double>());
    return result;
}

VarVector VarVector::operator-(const VarVector &rhs) const
{
    VarVector result(mData.size());
    std::transform(mData.begin(), mData.end(), rhs.mData.begin(), result.mData.begin(), std::minus<double>());
    return result;
}

VarVector VarVector::operator*(double rhs) const
{
    VarVector result(mData.size());
    std::transform(mData.begin(), mData.end(), result.mData.begin(), [rhs](double x) { return x * rhs; });
    return result;
}

VarVector operator*(double lhs, const VarVector &rhs)
{
    VarVector result(rhs.mData.size());
    std::transform(rhs.mData.begin(), rhs.mData.end(), result.mData.begin(), [lhs](double x) { return x * lhs; });
    return result;
}

VarVector VarVector::operator/(double rhs) const
{
    VarVector result(mData.size());
    std::transform(mData.begin(), mData.end(), result.mData.begin(), [rhs](double x) { return x / rhs; });
    return result;
}

std::ostream &operator<<(std::ostream &os, const VarVector &vector)
{

     for (std::size_t i = 0; i < vector.mData.size(); i++) {
        os << vector[i] << (i == vector.mData.size() - 1 ? "" : " ");
    }

    return os;
}