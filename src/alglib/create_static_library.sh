#!/bin/bash
wget http://www.alglib.net/translator/re/alglib-3.19.0.cpp.gpl.tgz
tar -xvzf alglib-3.19.0.cpp.gpl.tgz
rm -r ./include
rm -r ./lib
mkdir ./include
mkdir ./lib
mv ./alglib-cpp/src/* ./include
cd ./include
g++ -c *.cpp
ar rcs alglib.a *.o
rm -rf *.o
mv ./alglib.a ../lib/libalglib.a
echo "alglib Installed successfully!"
cd ~
