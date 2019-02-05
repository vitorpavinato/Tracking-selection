#!/bin/bash

for i in {1..1000}

do
  less selsim.$i.e* >> errorlogs.txt
done
