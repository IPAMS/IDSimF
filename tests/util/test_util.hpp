#ifndef IDSIMF_CPP_TEST_UTIL_HPP
#define IDSIMF_CPP_TEST_UTIL_HPP

#include "Core_vector.hpp"
#include <string>
#include <sstream>
#include <fstream>

class VectorApproxMatcher : public Catch::MatcherBase<Core::Vector> {
    private:
        Core::Vector target_;
        double epsilon_;
        mutable bool hasMatched_ = false;

    public:
        explicit VectorApproxMatcher(const Core::Vector& target, double epsilon= 1e-6):
            target_(target),
            epsilon_(epsilon){}

        bool match(const Core::Vector& vec) const override{
            hasMatched_ = std::abs(vec.x() - target_.x()) < epsilon_ &&
                          std::abs(vec.y() - target_.y()) < epsilon_ &&
                          std::abs(vec.z() - target_.z()) < epsilon_;

            return hasMatched_;
        }

        std::string describe() const override {
            std::stringstream ss;
            if (hasMatched_){
                ss << "is approximately equal to Core::Vector("
                    << target_.x() << ", " << target_.y() << ", " << target_.z() << ") within " << epsilon_;
            }
            else {
                ss << "differs from Core::Vector("
                    << target_.x() << ", " << target_.y() << ", " << target_.z() << ") with epsilon=" << epsilon_;
            }
            return ss.str();
        }
};

inline VectorApproxMatcher ApproxEqual(const Core::Vector& target, double epsilon = 1e-6) {
    return VectorApproxMatcher(target, epsilon);
};

// Vector equality means, exact, floating point equality here, thus deactivate
// floating point comparison warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
inline bool isExactDoubleEqual(double lhs, double rhs){
    return lhs == rhs;
}
#pragma GCC diagnostic pop

inline std::string readTextFile(std::string filename){
    std::ifstream inputFilestream(filename);
    std::string result(
            (std::istreambuf_iterator<char>(inputFilestream)),
            (std::istreambuf_iterator<char>()));

    return result;
}

#endif //IDSIMF_CPP_TEST_UTIL_HPP
