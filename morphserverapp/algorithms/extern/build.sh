#!/bin/bash
cd global-alignment
make clean
make all
cd ../qcp
make clean
make all
cd ../castellana-pevzner
make clean
make all
