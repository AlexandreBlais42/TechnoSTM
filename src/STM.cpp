#include "STM.h"

STM::STM(const uint8_t deviceAddr, const std::string devicePath,
         std::array<uint8_t, 4> pins)
    : Aiguille(deviceAddr), Plateforme(devicePath), StepMotor(0, pins) {
  state = Initialize;
}

void STM::start(){
  int16_t voltageAiguille = 0;
  while (true){
    std::cout << "État : " << std::to_string(state) << std::endl << std::endl;
    switch (state){
      case Initialize:
        StepMotor::goToRelative(-40);
        state = Find_sample;
        break;

      case Find_sample:
        Plateforme::setPositionRelative(0, 0, 1);
        if (Plateforme::position.z >= UINT16_MAX / 2){
          state = Lower_motor;
        } else if (Aiguille::readVoltage() >= AIGUILLE_THRESHOLD_VOLTAGE){ 
          state = Mesure_height;
        }
        break;

      case Lower_motor:
        Plateforme::setPositionAbsolute(0, 0, 0);
        delay_ms(2); // Laisser le temps à la plateforme de descendre
        StepMotor::goToRelative(1);
        state = Find_sample;
        break;

      case Mesure_height:
        voltageAiguille = Aiguille::readVoltage();
        if (abs(voltageAiguille - AIGUILLE_CONSTANT_CURRENT_VOLTAGE) < 300){
          state = Save_pixel;
        } else if ( voltageAiguille < AIGUILLE_CONSTANT_CURRENT_VOLTAGE){ // Le matériel est trop loin
          Plateforme::setPositionRelative(0, 0, 1);
        } else {                                                          // Le matériel est trop proche
          Plateforme::setPositionRelative(0, 0, -1);
        }
        break;

      case Save_pixel:
        // @todo Code pour savegarder le pixel sur l'image
        break;

      case Goto_next_coordinate:
        // @todo Code pour aller à la prochaine coordonnée
        break;
    }
  }
}

void STM::acquirePixelAtConstantHeight() {}

void STM::acquirePixelAtConstantCurrent() {}

void STM::goToNextPosition() {}

void STM::fixImage() {}

void STM::exportImage(const std::string filename) {}
