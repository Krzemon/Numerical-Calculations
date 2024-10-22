#ifndef TYPES_HPP
#define TYPES_HPP

#include <cmath>
#include <vector>

#include <fstream>

class VarVector
{
private:
    std::vector<double> mData;

public:
    VarVector() = default;
    VarVector(size_t size) : mData(size) {}
    VarVector(std::vector<double> data) : mData{data} {}
    VarVector(std::initializer_list<double> data) : mData{data} {}

    std::vector<double>::iterator begin() { return mData.begin(); }
    std::vector<double>::iterator end() { return mData.end(); }

    double norm() const;

    double operator[](size_t index) const { return mData[index]; }
    friend std::ostream &operator<<(std::ostream &os, const VarVector &vector);

    VarVector operator+(const VarVector &rhs) const;
    VarVector operator-(const VarVector &rhs) const;
    VarVector operator*(double rhs) const;
    friend VarVector operator*(double lhs, const VarVector &rhs);
    VarVector operator/(double rhs) const;
};

class State
{
private:
    double time;
    VarVector position;
    VarVector velocity;

public:
    State() = default;
    State(double time, VarVector position, VarVector velocity) : time{time}, position{position}, velocity{velocity} {}

    double getTime() const { return time; }
    VarVector getPosition() const { return position; }
    VarVector getVelocity() const { return velocity; }

    friend std::ostream &operator<<(std::ostream &os, const State &state)
    {
        os << state.time << " " << state.position << " " << state.velocity;
        return os;
    }
};

#endif // TYPES_HPP