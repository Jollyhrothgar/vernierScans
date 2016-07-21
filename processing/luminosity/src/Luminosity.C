#include "Luminosity.h"
#include <iostream>

Luminosity::Luminosity() {
  std::cout << "Luminosity instantiated at " << this << std::endl;
}

Luminosity::~Luminosity() {
  std::cout << "Destroying Luminosity from " << this << std::endl;
}
