// License:

#ifndef SPEAR_FINGERPRINT_HPP
#define SPEAR_FINGERPRINT_HPP

#include <vector>
#include <numeric>
#include <unordered_map>

namespace Spear {

class Fingerprint {
public:

    using iterator = std::unordered_map<size_t, size_t>::const_iterator;

    size_t increase(size_t i) {
        auto iter = sparse_vector_.find(i);

        if (iter == sparse_vector_.end()) {
            iter = sparse_vector_.insert({i, 0}).first;
        }

        ++iter->second;
        return i;
    }

    size_t value(size_t i) const {
        auto iter = sparse_vector_.find(i);
        if (iter == sparse_vector_.end()) {
            return 0;
        }

        return iter->second;
    }

    std::vector<bool> as_bit_vector(size_t length) const {
        std::vector<bool> bit_vector(false, length);

        for (auto&& bit : sparse_vector_) {
            bit_vector[bit.first % length] = true;
        }

        return bit_vector;
    }

    size_t size() const {
        return sparse_vector_.size();
    }

    iterator begin() const {
        return sparse_vector_.cbegin();
    }

    iterator end() const {
        return sparse_vector_.cend();
    }

private:
    std::unordered_map<size_t, size_t> sparse_vector_;

};

}

#endif
