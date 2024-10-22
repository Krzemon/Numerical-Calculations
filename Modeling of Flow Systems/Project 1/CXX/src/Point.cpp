// Point.cpp
#include "Point.h"
#include "Obstacle.h"

Point::Point(double x, double y) : x(x), y(y) {}

Point::Point() : x(0.0), y(0.0) {}

double Point::getX() const {
    return x;
}

double Point::getY() const {
    return y;
}

double Point::distanceTo(const Point& other) const {
    return std::hypot(x - other.x, y - other.y);
}

double Point::distanceToObstacle(const Obstacle& obstacle) const {
    // Sprawdź odległość od przeszkody
    double minDistance = std::numeric_limits<double>::max();

    // Użyj funkcji getObstaclePoints(), aby uzyskać dostęp do obstaclePoints
    const std::vector<Point>& obstaclePoints = obstacle.getObstaclePoints();

    for (const auto& obstaclePoint : obstaclePoints) {
        double distance = distanceTo(obstaclePoint);
        if (distance < minDistance) {
            minDistance = distance;
        }
    }

    return minDistance;
}
