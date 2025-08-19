#ifndef RAY2ELLIPSOID_HPP
#define RAY2ELLIPSOID_HPP

#include <Eigen/Dense>
#include <cmath>

namespace s2r {
    namespace utils {

    bool RayToPlane(const Eigen::Vector3d &rayOrigin,
                      const Eigen::Vector3d &rayDirection, const Eigen::Vector4d& plane,
                      double* t) {
      double rayDotPlaneNormal = rayDirection.dot(plane.head<3>());
      if (std::abs(rayDotPlaneNormal) < 1e-6) {
        return false; // Ray is parallel to the plane
      }
      double d = plane[3] + plane.head<3>().dot(rayOrigin);
      *t = -d / rayDotPlaneNormal;
      if (*t < 0) {
        return false; // Intersection is behind the ray origin
      }
      return true; // Intersection point is in front of the ray origin
    }

    bool RayToEllipsoid(const Eigen::Vector3d &rayOrigin,
                        const Eigen::Vector3d &rayDirection, double semiXAxis,
                        double semiYAxis, double semiZAxis,
                        double* t) {
      Eigen::Vector3d scaledOrigin, scaledDirection;
      scaledOrigin << rayOrigin.x() / semiXAxis, rayOrigin.y() / semiYAxis,
          rayOrigin.z() / semiZAxis;
      scaledDirection << rayDirection.x() / semiXAxis,
          rayDirection.y() / semiYAxis, rayDirection.z() / semiZAxis;
      double a = scaledDirection.squaredNorm();
      double b = 2*scaledOrigin.dot(scaledDirection);
      double c = scaledOrigin.squaredNorm() - 1;

      double discriminant = b * b - 4 * a * c;

      if (discriminant < 0) {
        return false;
      }

      double sqrtDiscriminant = std::sqrt(discriminant);
      double t1 = (-b - sqrtDiscriminant) / (2 * a);
      double t2 = (-b + sqrtDiscriminant) / (2 * a);

      if (t2 < 0) {
        return false;
      }

      *t = (t1 < 0) ? t2 : t1;

      return true;
    }
    }
}

#endif