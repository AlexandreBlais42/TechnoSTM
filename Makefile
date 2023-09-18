CXX = g++
CXXFlags = -O2 -Wall

programme: src/main.cpp src/Plateforme.cpp src/Vector3D.cpp  src/Serial.cpp src/Utils.cpp src/StepMotor.cpp src/GPIO.cpp
	$(CXX) -Isrc/include/ $^ -o build/programme $(CXXFlags)
