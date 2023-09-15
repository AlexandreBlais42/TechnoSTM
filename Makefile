CXX = g++
CXXFlags = -O2 -Wall

programme: src/main.cpp src/Plateforme.cpp src/Vector3D.cpp src/Piezo.cpp
	g++ Isrc/include/ $^ -o build/programme
