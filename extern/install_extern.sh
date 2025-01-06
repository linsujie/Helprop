#!/bin/bash

current_file=`readlink -f $0`
current_dir=`dirname ${current_file}`

cd ${current_dir}
mkdir -p build/mongo-c-driver; cd build/mongo-c-driver
cmake -DCMAKE_INSTALL_PREFIX=../../ -DENABLE_MONGOC=false ../../mongo-c-driver
make
make install

cd ${current_dir}
mkdir -p build/reflect-cpp; cd build/reflect-cpp
cmake ../../reflect-cpp
cmake -DCMAKE_INSTALL_PREFIX=../../ -DREFLECTCPP_BSON=true -DREFLECTCPP_USE_VCPKG=false .
make
make install
