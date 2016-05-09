#!/bin/sh

make clean-all

make thirdParty 

make fedora_fpic
make -f Makefile.lib fedora_fpic
