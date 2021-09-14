#!/bin/bash

g++ -W -Wall -O2 Fit_Gain_Simu.cc -o Fit_Gain_Simu.exe `root-config --cflags --libs` -lRooFitCore -lRooFit
