// License:

#ifndef SPEAR_UTILITIES_HPP
#define SPEAR_UTILITIES_HPP

#include "spear/exports.hpp"

#include <vector>
#include <sstream>
#include <string>

namespace Spear {

std::string SPEAR_EXPORT base64_encode(const std::string& input);

std::string SPEAR_EXPORT base64_decode(const std::string& input);

std::string SPEAR_EXPORT encode_bitvec(const std::vector<bool>& bv);

std::vector<bool> SPEAR_EXPORT decode_bitvec(const std::string& input);

//! Packs an integer and outputs it to a stream
void append_packed_int(std::stringstream &ss, std::uint32_t num) {
    int nbytes, bix;
    unsigned int val, res;
    char tc;

    res = num;
    while (1) {
        if (res < (1 << 7)) {
            val = (res << 1);
            nbytes = 1;
            break;
        }
        res -= (1 << 7);
        if (res < (1 << 14)) {
            val = ((res << 2) | 1);
            nbytes = 2;
            break;
        }
        res -= (1 << 14);
        if (res < (1 << 21)) {
            val = ((res << 3) | 3);
            nbytes = 3;
            break;
        }
        res -= (1 << 21);
        if (res < (1 << 29)) {
            val = ((res << 3) | 7);
            nbytes = 4;
            break;
        } else {
            throw std::range_error("ERROR: Integer too big to pack\n");
        }
    }

    for (bix = 0; bix < nbytes; bix++) {
        tc = (char)(val & 255);
        ss.write(&tc, 1);
        val >>= 8;
    }
}

//! Reads an integer from a stream in packed format and returns the result.
std::uint32_t read_packed_int(std::istringstream &ss) {
    std::uint32_t val, num;
    int shift, offset;
    unsigned char tmp;
    ss.read((char*)&tmp, sizeof(tmp));
    val = tmp;
    offset = 0;
    if ((val & 1) == 0) {
        shift = 1;
    } else if ((val & 3) == 1) {
        ss.read((char *)&tmp, sizeof(tmp));
        val |= (tmp << 8);
        shift = 2;
        offset = (1 << 7);
    } else if ((val & 7) == 3) {
        ss.read((char *)&tmp, sizeof(tmp));
        val |= (tmp << 8);
        ss.read((char *)&tmp, sizeof(tmp));
        val |= (tmp << 16);
        shift = 3;
        offset = (1 << 7) + (1 << 14);
    } else {
        ss.read((char *)&tmp, sizeof(tmp));
        val |= (tmp << 8);
        ss.read((char *)&tmp, sizeof(tmp));
        val |= (tmp << 16);
        ss.read((char *)&tmp, sizeof(tmp));
        val |= (tmp << 24);
        shift = 3;
        offset = (1 << 7) + (1 << 14) + (1 << 21);
    }
    num = (val >> shift) + offset;
    return num;
}

};

#endif
