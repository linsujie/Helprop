#!/usr/bin/env bash
#
# MongoDB C Driver uninstall program, generated with CMake
#
# Copyright 2009-present MongoDB, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License")
#
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -eu

__rmfile() {
    set -eu
    abs=$__prefix/$1
    printf "Remove file %s: " "$abs"
    if test -f "$abs" || test -L "$abs"
    then
        rm -- "$abs"
        echo "ok"
    else
        echo "skipped: not present"
    fi
}

__rmdir() {
    set -eu
    abs=$__prefix/$1
    printf "Remove directory %s: " "$abs"
    if test -d "$abs"
    then
        list=$(ls --almost-all "$abs")
        if test "$list" = ""
        then
            rmdir -- "$abs"
            echo "ok"
        else
            echo "skipped: not empty"
        fi
    else
        echo "skipped: not present"
    fi
}

__prefix=${DESTDIR:-}/home/jiang/helprop2/extern

__rmfile "lib/libbson-static-1.0.a"
__rmfile "lib/cmake/bson-1.0/bson_static-targets.cmake"
__rmfile "lib/cmake/bson-1.0/bson_static-targets-relwithdebinfo.cmake"
__rmfile "lib/pkgconfig/libbson-static-1.0.pc"
__rmfile "lib/libbson-1.0.so.0.0.0"
__rmfile "lib/libbson-1.0.so.0"
__rmfile "lib/libbson-1.0.so"
__rmfile "lib/cmake/bson-1.0/bson_shared-targets.cmake"
__rmfile "lib/cmake/bson-1.0/bson_shared-targets-relwithdebinfo.cmake"
__rmfile "lib/pkgconfig/libbson-1.0.pc"
__rmfile "include/libbson-1.0/bson/bson-md5.h"
__rmfile "include/libbson-1.0/bson/bson-types.h"
__rmfile "include/libbson-1.0/bson/bson-endian.h"
__rmfile "include/libbson-1.0/bson/bson-decimal128.h"
__rmfile "include/libbson-1.0/bson/bson-value.h"
__rmfile "include/libbson-1.0/bson/bson-writer.h"
__rmfile "include/libbson-1.0/bson/bson-iter.h"
__rmfile "include/libbson-1.0/bson/bson-context.h"
__rmfile "include/libbson-1.0/bson/bson-utf8.h"
__rmfile "include/libbson-1.0/bson/bson-error.h"
__rmfile "include/libbson-1.0/bson/bson-string.h"
__rmfile "include/libbson-1.0/bson/bson-cmp.h"
__rmfile "include/libbson-1.0/bson/bson-compat.h"
__rmfile "include/libbson-1.0/bson/bson-memory.h"
__rmfile "include/libbson-1.0/bson/bson.h"
__rmfile "include/libbson-1.0/bson/bson-prelude.h"
__rmfile "include/libbson-1.0/bson/bson-oid.h"
__rmfile "include/libbson-1.0/bson/bson-reader.h"
__rmfile "include/libbson-1.0/bson/bson-keys.h"
__rmfile "include/libbson-1.0/bson/bson-macros.h"
__rmfile "include/libbson-1.0/bson/bson-version-functions.h"
__rmfile "include/libbson-1.0/bson/bcon.h"
__rmfile "include/libbson-1.0/bson/bson-clock.h"
__rmfile "include/libbson-1.0/bson/bson-atomic.h"
__rmfile "include/libbson-1.0/bson/bson-json.h"
__rmfile "include/libbson-1.0/bson/bson-config.h"
__rmfile "include/libbson-1.0/bson/bson-version.h"
__rmfile "include/libbson-1.0/bson.h"
__rmfile "lib/cmake/bson-1.0/bson-1.0-config.cmake"
__rmfile "lib/cmake/bson-1.0/bson-1.0-config-version.cmake"
__rmfile "lib/cmake/bson-1.0/bson-targets.cmake"
__rmfile "lib/cmake/libbson-1.0/libbson-1.0-config.cmake"
__rmfile "lib/cmake/libbson-1.0/libbson-1.0-config-version.cmake"
__rmfile "lib/cmake/libbson-static-1.0/libbson-static-1.0-config.cmake"
__rmfile "lib/cmake/libbson-static-1.0/libbson-static-1.0-config-version.cmake"
__rmfile "share/mongo-c-driver/COPYING"
__rmfile "share/mongo-c-driver/NEWS"
__rmfile "share/mongo-c-driver/README.rst"
__rmfile "share/mongo-c-driver/THIRD_PARTY_NOTICES"
__rmfile "share/mongo-c-driver/uninstall-bson.sh"
__rmdir "share/mongo-c-driver"
__rmdir "lib/pkgconfig"
__rmdir "lib/cmake/libbson-static-1.0"
__rmdir "lib/cmake/libbson-1.0"
__rmdir "lib/cmake/bson-1.0"
__rmdir "lib/cmake"
__rmdir "include/libbson-1.0/bson"
__rmdir "include/libbson-1.0"
