Fit_Gain_simu.exe : Fit_Gain_Simu.cc
	g++ -W -Wall -O2 $^ -o $@ `root-config --cflags --libs` -lRooFitCore -lRooFit 
Dead_time.exe : Dead_time.cc Dead_time.hh
	g++ -W -Wall -O2 $^ -o $@ `root-config --cflags --libs` -lRooFitCore -lRooFit
