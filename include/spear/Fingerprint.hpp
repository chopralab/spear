// License:

#ifndef SPEAR_FINGERPRINT_HPP
#define SPEAR_FINGERPRINT_HPP

#include <unordered_map>

namespace Spear {

class Fingerprint {
public:

    using iterator = std::unordered_map<size_t, unsigned char>::const_iterator;

    Fingerprint(size_t length) : length_(length) {}

    size_t increase(size_t i) {
        auto location = i % length_;
        auto iter = sparse_vector_.find(location);

        if (iter == sparse_vector_.end()) {
            iter = sparse_vector_.insert({i % length_, 0}).first;
        }

        ++iter->second;
        return location;
    }

    size_t value(size_t i) const {
        auto iter = sparse_vector_.find(i % length_);
        if (iter == sparse_vector_.end()) {
            return 0;
        }

        return iter->second;
    }

    size_t length() const {
        return length_;
    }

    size_t non_zero_values() const {
        return sparse_vector_.size();
    }

    iterator begin() const {
        return sparse_vector_.cbegin();
    }

    iterator end() const {
        return sparse_vector_.cend();
    }

private:

    size_t length_;
    std::unordered_map<size_t, unsigned char> sparse_vector_;

};

}

#endif
