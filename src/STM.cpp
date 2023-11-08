#include "STM.h"

STM::STM(const uint8_t deviceAddr, const std::string devicePath,
         std::array<uint8_t, 4> pins, const uint32_t resolutionX,
         const uint32_t resolutionY, const uint32_t scale)
    : Aiguille(deviceAddr), Plateforme(devicePath), StepMotor(0, pins),
      resolutionX(resolutionX), resolutionY(resolutionY), scale(scale) {
  state = Initialize;
  voltageAiguille = 0;
}

void STM::start() {
  int16_t voltageAiguille = 0;
  while (true) {
    /*  @todo Implémenter select() pour pouvoir recevoir des données d'un socket de manière non bloquante
     */ 
    //select(nfds, readfdset, nullptr, nullptr, timeout);

    switch (state) {
    case Initialize:
      DEBUG << "État : " << "Initialize" << std::endl;
      StepMotor::setPositionRelative(-40);
      Plateforme::setPositionAbsolute(0, 0, 0);
      state = Find_sample;
      DEBUG << "État : " << "Find_sample" << std::endl;
      break;

    case Find_sample:
      Plateforme::setPositionRelative(0, 0, 1);
      if (Plateforme::position.z >= UINT16_MAX / 2) {
        state = Lower_motor;
        DEBUG << "État : " << "Lower_motor" << std::endl;
      } else if (Aiguille::readVoltage() >= AIGUILLE_THRESHOLD_VOLTAGE) {
        state = Mesure_height;
        DEBUG << "État : " << "Mesure_height" << std::endl;
      }
      break;

    case Lower_motor:
      Plateforme::setPositionAbsolute(0, 0, 0);
      delay_ms(2); // Laisser le temps à la plateforme de descendre
      StepMotor::setPositionRelative(1);
      state = Find_sample;
      DEBUG << "État : " << "Find_sample" << std::endl;
      break;

    case Mesure_height:
      voltageAiguille = Aiguille::readVoltage();
      if (abs(voltageAiguille - AIGUILLE_CONSTANT_CURRENT_VOLTAGE) < 300) {
        state = Save_pixel;
      DEBUG << "État : " << "Save_pixel" << std::endl;
      } else if (voltageAiguille <
                 AIGUILLE_CONSTANT_CURRENT_VOLTAGE) { // Le matériel est trop
                                                      // loin
        Plateforme::setPositionRelative(0, 0, 1);
      } else { // Le matériel est trop proche
        Plateforme::setPositionRelative(0, 0, -1);
      }
      break;

    case Save_pixel:
      image.setPixel(Plateforme::position.x, Plateforme::position.y,
                     Plateforme::position.z);
      state = Goto_next_coordinate;
      DEBUG << "État : " << "Goto_next_coordinate" << std::endl;
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
