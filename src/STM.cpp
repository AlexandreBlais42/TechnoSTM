#include "STM.h"

STM::STM(const uint8_t deviceAddr, const std::string devicePath,
         std::array<uint8_t, 4> pins, const uint32_t resolutionX,
         const uint32_t resolutionY, const uint32_t scale)
    : Aiguille(deviceAddr), Plateforme(devicePath), StepMotor(0, pins),
      resolutionX(resolutionX), resolutionY(resolutionY), scale(scale) {
  image.initialize(resolutionX, resolutionY, "scan", "m");
  state = Initialize;
  voltageAiguille = 0;
}

void STM::start() {
  int16_t voltageAiguille = 0;
  bool done = false;
  while (done == false) {
    /*  @todo Implémenter select() pour pouvoir recevoir des données d'un socket
     * de manière non bloquante
     */
    // select(nfds, readfdset, nullptr, nullptr, timeout);

    switch (state) {
    case Initialize:
      DEBUG << "État : "
            << "Initialize" << std::endl;
      StepMotor::setPositionRelative(-40);
      Plateforme::setPositionAbsolute(0, 0, 0);
      state = Find_sample;
      DEBUG << "État : "
            << "Find_sample" << std::endl;
      break;

    case Find_sample:
      Plateforme::setPositionRelative(0, 0, 1);
      if (Plateforme::position.z >= UINT16_MAX / 2) {
        state = Lower_motor;
        DEBUG << "État : "
              << "Lower_motor" << std::endl;
      } else if (Aiguille::readVoltage() >= AIGUILLE_THRESHOLD_VOLTAGE) {
        state = Mesure_height;
        DEBUG << "État : "
              << "Mesure_height" << std::endl;
      }
      break;

    case Lower_motor:
      Plateforme::setPositionAbsolute(0, 0, 0);
      delay_ms(2); // Laisser le temps à la plateforme de descendre
      StepMotor::setPositionRelative(1);
      state = Find_sample;
      DEBUG << "État : "
            << "Find_sample" << std::endl;
      break;

    case Mesure_height:
      voltageAiguille = Aiguille::readVoltage();
      if (abs(voltageAiguille - AIGUILLE_CONSTANT_CURRENT_VOLTAGE) < 300) {
        state = Save_pixel;
        DEBUG << "État : "
              << "Save_pixel" << std::endl;
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
      DEBUG << "État : "
            << "Goto_next_coordinate" << std::endl;
      break;

    case Goto_next_coordinate:
      if (goToNextPosition() == false) {
        exportImage("scan.gfs");
        done = true;
      }
      state = Mesure_height;
      break;
    }
  }
}

void STM::acquirePixelAtConstantHeight() {}

void STM::acquirePixelAtConstantCurrent() {}

bool STM::goToNextPosition() {
  Vector3D<uint16_t> position = Plateforme::position;
  bool direction = position.y & 0x01 ? false : true;
  // Faux est vert la gauche et vrai vert la droite.
  // On veut scanner en faisant un serpentin
  if (position.x == resolutionX - 1 && position.y == resolutionY - 1) {
    // Coin en bas à droite atteint, fin du scan
    return false;
  }
  if (position.x == resolutionX - 1 && direction == true) {
    position.y += 1;
  } else if (position.x == 0 && direction == false) {
    position.y += 1;
  } else {
    position.x += direction ? 1 : -1;
  }
  Plateforme::setPositionAbsolute(position);
  return true;
}

void STM::exportImage(const std::string filename) {
  FILE *file = fopen(filename.c_str(), "wb");
  std::vector<uint8_t> imageData = image.getGSFFile();
  fwrite(imageData.data(), sizeof(uint8_t), imageData.size(), file);
  fclose(file);
  std::cout << filename << " written!" << std::endl;
}
