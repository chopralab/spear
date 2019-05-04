// License:

#ifndef SPEAR_DEFAULT_PARTIALCHARGE_HPP
#define SPEAR_DEFAULT_PARTIALCHARGE_HPP

#include "spear/PartialCharge.hpp"

namespace Spear {

class Molecule;

class SPEAR_EXPORT DefaultPartialCharge final : public PartialCharge {
public:

    DefaultPartialCharge(const Molecule& /*unused*/, std::vector<double> charges) :
        PartialCharge(std::move(charges))
    { }

    const std::string& name() const override {
        return name_;
    }

    double add_atom(size_t /*new_idx*/) override {
        push_back(0.0);
        return 0.0;
    }

    virtual void remove_atom(size_t idx) override {
        erase(this->begin() + static_cast<std::ptrdiff_t>(idx));
    }

    void set(size_t idx, double val) {
        PartialCharge::set(idx, val);
    }

private:
    const std::string name_ = "default_partial_charge";
};

}

#endif
