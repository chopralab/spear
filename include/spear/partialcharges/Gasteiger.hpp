// License:

#ifndef SPEAR_GASTEIGER_PARTIALCHARGE_HPP
#define SPEAR_GASTEIGER_PARTIALCHARGE_HPP

#include "spear/PartialCharge.hpp"

namespace Spear {

class Molecule;

class SPEAR_EXPORT GasteigerCharge final : public PartialCharge {
public:

    GasteigerCharge(const Molecule& mol);

    const std::string& name() const override {
        return name_;
    }

    double add_atom(size_t /*idx*/) override {
        push_back(0.0);
        init_();
        return 0.0;
    }

    virtual void remove_atom(size_t idx) override {
        erase(this->begin() + static_cast<std::ptrdiff_t>(idx));
        init_();
    }

private:
    void init_();
    const Molecule& mol_;
    const std::string name_ = "gasteiger";
};

}

#endif
