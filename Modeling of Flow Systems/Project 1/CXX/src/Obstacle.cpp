#include "Obstacle.h"
#include "LJ_Potential.h"
#include "Particle.h"

Obstacle::Obstacle() {}

Obstacle::Obstacle(const std::vector<Point>& obstacleShape, const double& mp, const double& sigma, const double& eps, const double& rcut)
    : mp(mp), sigma(sigma), eps(eps), rcut(rcut) {
    // Inicjalizuj przeszkodę na podstawie punktów cząstek
    for (const auto& obstaclePoint : obstacleShape) {
        obstaclePoints.emplace_back(obstaclePoint.getX(), obstaclePoint.getY());
    }
}

void Obstacle::addParticles(std::vector<Particle>& particles) {
    // Dodaj cząstki przeszkody
    for (const auto& obstaclePoint : obstaclePoints) {
        Particle particle(obstaclePoint.getX(), obstaclePoint.getY(), 0, mp);
        particle.xvel = particles.size();
        particle.yvel = 0;

        particles.push_back(particle);
    }
}

void Obstacle::reactToObstacle(std::vector<Particle>& particles) const {
    // React particles to the obstacle
    for (auto& particle : particles) {
        if (isInsideObstacle(Point(particle.xpos, particle.ypos))) {
            // Repulsive forces from the obstacle
            Point closestPoint = closestPointOnObstacle(Point(particle.xpos, particle.ypos));
            double dx = particle.xpos - closestPoint.getX();
            double dy = particle.ypos - closestPoint.getY();
            double dist = std::hypot(dx, dy);  // or use distanceTo function

            if (dist < particle.mass) {
                // Use Lennard-Jones potential for interaction
                LJ_Potential ljPotential = LJ_Potential::getInstance(sigma, eps, rcut);
                double forceMagnitude = ljPotential.calculateForce(dist);

                // Force components in x and y directions
                double forceX = forceMagnitude * (dx / dist);
                double forceY = forceMagnitude * (dy / dist);

                // Add force to the particle
                particle.xfor += forceX;
                particle.yfor += forceY;
            }
        }
    }
}

bool Obstacle::isInsideObstacle(const Point& point) const {
    // Point-in-polygon algorithm
    int i, j;
    bool inside = false;
    int verticesCount = obstaclePoints.size();

    for (i = 0, j = verticesCount - 1; i < verticesCount; j = i++) {
        if (((obstaclePoints[i].getY() > point.getY()) != (obstaclePoints[j].getY() > point.getY())) &&
            (point.getX() < (obstaclePoints[j].getX() - obstaclePoints[i].getX()) * (point.getY() - obstaclePoints[i].getY()) /
                             (obstaclePoints[j].getY() - obstaclePoints[i].getY()) + obstaclePoints[i].getX())) {
            inside = !inside;
        }
    }

    return inside;
}

bool Obstacle::isNearObstacle(const Point& point, double distance) const {
    // Check if the point is close to the obstacle (based on Euclidean distance)
    for (const auto& obstaclePoint : obstaclePoints) {
        double dist = std::hypot(obstaclePoint.getX() - point.getX(), obstaclePoint.getY() - point.getY());  // or use distanceTo function
        if (dist < distance) {
            return true;
        }
    }
    return false;
}

Point Obstacle::closestPointOnObstacle(const Point& point) const {
    double minDistance = std::numeric_limits<double>::max();
    Point closestPoint;

    for (const auto& obstaclePoint : obstaclePoints) {
        double distance = std::hypot(obstaclePoint.getX() - point.getX(), obstaclePoint.getY() - point.getY());  // or use distanceTo function
        if (distance < minDistance) {
            minDistance = distance;
            closestPoint = obstaclePoint;
        }
    }

    return closestPoint;
}

const std::vector<Point>& Obstacle::getObstaclePoints() const {
    return obstaclePoints;
}
