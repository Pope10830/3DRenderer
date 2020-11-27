#pragma once

#include <iostream>

struct Material {
	std::string name;
	Colour colour;
  TextureMap texture;
  bool hasTexture;
	Material();
	Material(std::string n, Colour c);
  Material(std::string n, TextureMap t);
};
