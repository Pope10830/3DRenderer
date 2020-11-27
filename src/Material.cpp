#include "Material.h"
#include <utility>

Material::Material() = default;
Material::Material(std::string n, Colour c) :
		name(std::move(n)),
		colour(c),
    hasTexture(false) {}
Material::Material(std::string n, TextureMap t) :
    name(std::move(n)),
    texture(t),
    hasTexture(true) {}
