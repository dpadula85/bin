#!/bin/bash

for i in `ls *.chk`;
do
    formchk -3 $i
done
