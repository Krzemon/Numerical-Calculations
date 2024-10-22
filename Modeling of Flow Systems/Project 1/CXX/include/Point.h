#ifndef POINT_H
#define POINT_H

#include "Includes.h"
#include "Obstacle.h"

class Obstacle;

class Point {
public:
    Point();
    Point(double x, double y);
    double getX() const;
    double getY() const;
    double distanceTo(const Point& other) const;
    double distanceToObstacle(const Obstacle& obstacle) const;  // Dodaj deklarację funkcji

private:
    double x, y;
};

#endif // POINT_H
