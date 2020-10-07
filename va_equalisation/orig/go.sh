g++ -std=c++11 -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService va_equalisation.cpp -o va_eq
#g++ -std=c++11 -Wall `root-config --cflags --libs` langauConv.cpp corrFac.cpp -o corrFac10
#g++ -std=c++11 -Wall `root-config --cflags --libs` src/langauConv.cpp src/corrFac.cpp -o corrFac09
#g++ -std=c++11 -Wall `root-config --cflags --libs` src/langauConv.cpp src/corrFac.cpp -o corrFac12
