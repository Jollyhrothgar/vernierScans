#include "CrossSection.h"
#include <iostream>

CrossSection::CrossSection() {
  std::cout << "CrossSection instantiated at " << this << std::endl;
}

CrossSection::~CrossSection() {
  std::cout << "Destroying CrossSection from " << this << std::endl;
}
