#!/usr/bin/env bash

# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/gnuconfig/config.* .

autoreconf -vfi

./configure --prefix=$PREFIX

make

make install
