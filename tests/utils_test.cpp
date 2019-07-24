#include "spear/utils.hpp"

#include <vector>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

TEST_CASE("Encode/decode bit vectors") {
    std::vector<bool> bv(2048, false);
    bv[1] = true;

    auto encoded = base64_encode(encode_bitvec(bv));
    CHECK(encoded == "4P///wAIAAABAAAAAvkd");

    auto new_bv = decode_bitvec(base64_decode(encoded));
    CHECK(new_bv == bv);
}
