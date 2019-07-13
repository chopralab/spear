#include "spear/utils.hpp"
#include <vector>

// Stolen from RDKit
std::string Spear::base64_encode(const std::string& input) {

    static unsigned char transTable[64] = {
      'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
      'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
      'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
      'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
      '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/'
    };

    size_t resSize;
    resSize = (4 * input.size()) / 3;

    if (resSize % 4 != 0) {
        resSize += 4 - (resSize % 4);
    }

    std::string res(resSize, '\0');

    size_t i = 0;
    size_t pos = 0;

    while (i < input.size()) {
        res[pos++] = transTable[input[i] >> 2];
        if (i + 1 < input.size()) {
            res[pos++] = transTable[((input[i] & 3) << 4) | (input[i + 1] >> 4)];
            if (i + 2 < input.size()) {
                res[pos++] =
                    transTable[((input[i + 1] & 0xF) << 2) | (input[i + 2] >> 6)];
                res[pos++] = transTable[input[i + 2] & 0x3F];
            } else {
                // single padding
                res[pos++] = transTable[((input[i + 1] & 0xF) << 2)];
                res[pos++] = '=';
            }
        } else {
            // double padding
            res[pos++] = transTable[((input[i] & 3) << 4)];
            res[pos++] = '=';
            res[pos++] = '=';
        }
        i += 3;
    }
    res[resSize] = 0;
    return res;
}

std::string Spear::base64_decode(const std::string& input) {
    std::vector<unsigned char> transTable(256, 0x80);

    unsigned char c;
    for (c = 'A'; c <= 'Z'; c++) transTable[c] = c - 'A';
    for (c = 'a'; c <= 'z'; c++) transTable[c] = c - 'a' + 26;
    for (c = '0'; c <= '9'; c++) transTable[c] = c - '0' + 52;
    transTable[static_cast<size_t>('+')] = 62;
    transTable[static_cast<size_t>('/')] = 63;

    size_t outLen = 3 * input.length() / 4;
    std::string res(outLen, '\0');

    size_t pos = 0;
    size_t i = 0;
    // decode 4 bytes at a time
    unsigned char block[4];
    int nInBlock = 0;
    while (i < input.size()) {
        auto c = input[i];

        // above we set 0x80 as the junk marker in the translation table
        if (!(transTable[c] & 0x80)) {
            block[nInBlock++] = transTable[c];
            if (nInBlock == 4) {
                // finished a block
                res[pos++] = (block[0] << 2) | (block[1] >> 4);
                res[pos++] = (block[1] << 4) | (block[2] >> 2);
                res[pos++] = (block[2] << 6) | block[3];
                nInBlock = 0;
            }
        }
        i++;
    }

    // okay, now there can be 2 or 3 chars remaining to be processed
    //  (before the padding)
    if (nInBlock > 1) {
        res[pos++] = (block[0] << 2) | (block[1] >> 4);
        if (nInBlock > 2) {
            res[pos++] = (block[1] << 4) | (block[2] >> 2);
            res[pos] = (block[2] << 6);
        }
    }
    return res;
}
