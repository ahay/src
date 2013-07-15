icc -openmp flexSIMD.cpp -I .





icc -openmp -inline-forceinline flexSIMD.cpp -I .


1013  icc flexSerial.cpp -I . -O3 -o flexSerial.exe
1.4 3.5 4.5 gflops for z,y,x loops

1014  icc flexSerial1.cpp -I . -O3 -o flexSerial1.exe
1.4 3.5 4.5 gflops



