#ifndef NPSAT_URF_STREAMLINE_SIMPLIFICATION_H
#define NPSAT_URF_STREAMLINE_SIMPLIFICATION_H

#include "my_structures.h"
#include "streamline_reader.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

struct SimplifiedPoint{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double v = 0.0;
    double age = 0.0;
};

inline double pointSegmentDistance(const StreamlinePoint& p,
                                   const StreamlinePoint& a,
                                   const StreamlinePoint& b) {
    const double abx = b.x - a.x;
    const double aby = b.y - a.y;
    const double abz = b.z - a.z;
    const double apx = p.x - a.x;
    const double apy = p.y - a.y;
    const double apz = p.z - a.z;
    const double ab2 = abx*abx + aby*aby + abz*abz;

    double t = 0.0;
    if (ab2 > 0.0) {
        t = (apx*abx + apy*aby + apz*abz) / ab2;
        if (t < 0.0) {
            t = 0.0;
        }
        else if (t > 1.0) {
            t = 1.0;
        }
    }

    const double cx = a.x + t*abx;
    const double cy = a.y + t*aby;
    const double cz = a.z + t*abz;
    const double dx = p.x - cx;
    const double dy = p.y - cy;
    const double dz = p.z - cz;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

inline void markDouglasPeuckerPoints(const std::vector<StreamlinePoint>& points,
                                     const unsigned int first,
                                     const unsigned int last,
                                     const double tolerance,
                                     std::vector<int>& keep) {
    if (last <= first + 1) {
        return;
    }

    double maxDistance = -1.0;
    unsigned int splitIndex = first;
    for (unsigned int i = first + 1; i < last; ++i) {
        const double distance = pointSegmentDistance(points[i], points[first], points[last]);
        if (distance > maxDistance) {
            maxDistance = distance;
            splitIndex = i;
        }
    }

    if (maxDistance > tolerance) {
        keep[splitIndex] = 1;
        markDouglasPeuckerPoints(points, first, splitIndex, tolerance, keep);
        markDouglasPeuckerPoints(points, splitIndex, last, tolerance, keep);
    }
}

inline bool calculateDetailedAges(const StreamlineTrajectory& trajectory,
                                  std::vector<double>& ages) {
    if (trajectory.samples.empty()) {
        return false;
    }

    ages.assign(trajectory.samples.size(), 0.0);
    for (unsigned int i = 1; i < trajectory.samples.size(); ++i) {
        const StreamlinePoint& a = trajectory.samples[i - 1];
        const StreamlinePoint& b = trajectory.samples[i];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double dz = b.z - a.z;
        const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (b.vmag <= 0.0) {
            return false;
        }
        ages[i] = ages[i - 1] + distance / b.vmag;
    }

    return true;
}

inline void writeSimplifiedStreamline(const StreamlineTrajectory& trajectory,
                                      const URFoptions& opt,
                                      std::ofstream& simplifiedFile,
                                      std::ofstream& discardFile) {
    if (!opt.simplifyStreamline) {
        return;
    }

    if (trajectory.samples.size() < 2) {
        writeDiscardedStreamline(discardFile, trajectory, "simplification_not_enough_samples");
        return;
    }

    std::vector<double> detailedAges;
    if (!calculateDetailedAges(trajectory, detailedAges)) {
        writeDiscardedStreamline(discardFile, trajectory, "simplification_invalid_velocity");
        return;
    }

    std::vector<int> keep(trajectory.samples.size(), 0);
    keep.front() = 1;
    keep.back() = 1;
    markDouglasPeuckerPoints(trajectory.samples, 0,
                             static_cast<unsigned int>(trajectory.samples.size() - 1),
                             opt.simplifyTolerance, keep);

    std::vector<unsigned int> indices;
    for (unsigned int i = 0; i < keep.size(); ++i) {
        if (keep[i]) {
            indices.push_back(i);
        }
    }

    std::vector<SimplifiedPoint> simplified(indices.size());
    for (unsigned int i = 0; i < indices.size(); ++i) {
        const StreamlinePoint& point = trajectory.samples[indices[i]];
        simplified[i].x = point.x;
        simplified[i].y = point.y;
        simplified[i].z = point.z;
        simplified[i].age = detailedAges[indices[i]];
    }

    for (unsigned int i = 0; i + 1 < simplified.size(); ++i) {
        const double dx = simplified[i + 1].x - simplified[i].x;
        const double dy = simplified[i + 1].y - simplified[i].y;
        const double dz = simplified[i + 1].z - simplified[i].z;
        const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        const double dt = simplified[i + 1].age - simplified[i].age;
        simplified[i].v = dt > 0.0 ? distance / dt : 0.0;
    }
    simplified.back().v = simplified.size() > 1 ? simplified[simplified.size() - 2].v : 0.0;

    for (unsigned int i = 0; i < simplified.size(); ++i) {
        simplifiedFile << static_cast<unsigned long long>(trajectory.Eid) << ","
                       << static_cast<unsigned long long>(trajectory.Sid) << ","
                       << std::setprecision(10) << std::fixed
                       << simplified[i].x << ","
                       << simplified[i].y << ","
                       << simplified[i].z << ","
                       << simplified[i].v << ","
                       << simplified[i].age << std::endl;
    }
}

#endif //NPSAT_URF_STREAMLINE_SIMPLIFICATION_H
