#!/usr/bin/env bash

# To run it, you will need to generate a compilation database first. CMake can
# generate one for you using `cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`. It is
# better to use the same version of clang as clang-tidy when configuring with
# cmake.
#
# CMake will generate a `compile_commands.json` file, that you should copy/link
# in the root of chemfiles sources. You can then run this script from the root
# as ./scripts/tidy.sh

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/" && pwd)

if [[ ! -f "$ROOT/compile_commands.json" ]]; then
    echo "Missing compile_commands.json. See the script for how to generate it"
    exit 1
fi

HEADERS='include/spear/*.hpp'
CHECKS="clang-*,\
        modernize-*,\
        bugprone-*,\
        misc-*,\
            -misc-noexcept-move-constructor,\
        performance-*,\
        readability-*,\
            -readability-container-size-empty,\
            -readability-implicit-bool-cast,\
            -readability-else-after-return,\
	ShortStatementLines=5,\
        cppcoreguidelines-*,\
            -cppcoreguidelines-pro-type-union-access,\
            -cppcoreguidelines-pro-type-member-init,\
            -cppcoreguidelines-pro-bounds-pointer-arithmetic,\
            -cppcoreguidelines-pro-type-reinterpret-cast,\
            -cppcoreguidelines-pro-bounds-constant-array-index,\
            -cppcoreguidelines-pro-bounds-array-to-pointer-decay,\
            -cppcoreguidelines-pro-type-vararg,\
            "

pushd $ROOT > /dev/null

for file in $(find src -name '*.?pp')
do
    echo "===== checking ${file} ======"
    clang-tidy-6.0 -p . -checks="${CHECKS}" -analyze-temporary-dtors -header-filter="${HEADERS}" $file \
	-config="{CheckOptions: [{key: readability-braces-around-statements.ShortStatementLines, value: 2}]}"
done

popd > /dev/null
