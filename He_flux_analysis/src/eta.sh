g++ -std=c++11 -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o eta.exe
#g++ -std=c++11 -Wall `root-config --cflags --libs` -I$DMPSWSYS/include -L$DMPSWSYS/lib -lDmpEvent -lDmpService eta.cc -o ../submit/eta.exe
