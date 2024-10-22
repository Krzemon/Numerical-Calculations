#ifndef OBSTACLE_H
#define OBSTACLE_H

#include "Includes.h"
#include "Point.h"
#include "Particle.h"
class Particle;
class Point;

class Obstacle {
public:
    Obstacle();
    Obstacle(const std::vector<Point>& obstacleShape, const double& mp, const double& sigma, const double& eps, const double& rcut);
    void addParticles(std::vector<Particle>& particles);
    void reactToObstacle(std::vector<Particle>& particles) const;
    bool isInsideObstacle(const Point& point) const;
    bool isNearObstacle(const Point& point, double distance) const;
    Point closestPointOnObstacle(const Point& point) const;
    const std::vector<Point>& getObstaclePoints() const;
private:
    std::vector<Point> obstaclePoints;
    double mp;
    double sigma;
    double eps;
    double rcut;
};

#endif // OBSTACLE_H