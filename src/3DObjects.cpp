#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <ModelTriangle.h>
#include "LineDrawing.cpp"
#include "Material.cpp"
#include <glm/gtx/io.hpp>
#include <math.h>
#include <RayTriangleIntersection.h>

#define WIDTH 320
#define HEIGHT 240

using namespace std;

void readMTL(vector<Material> &materials, string filename, string foldername) {
	string line;
	ifstream file;
	file.open(filename);

	string name;

	while(getline(file, line)) {
		vector<string> splitline = split(line, ' ');
		if (splitline[0] == "newmtl") {
			name = splitline[1];
		} else if (splitline[0] == "Kd") {
			int r = stof(splitline[1]) * 255;
			int g = stof(splitline[2]) * 255;
			int b = stof(splitline[3]) * 255;

			materials.push_back(Material(name, Colour(name, r, g, b)));
		} else if (splitline[0] == "map_Kd") {
			materials.back().texture = TextureMap(foldername + splitline[1]);
			materials.back().hasTexture = true;
		}
	}

	file.close();
}

void readOBJ(vector<ModelTriangle> &triangles, string filename, float scale, vector<Material> materials) {
	vector<glm::vec3> points;
	vector<glm::vec2> texturePoints;
	vector<glm::vec3> normals;
	string line;
	ifstream file;
	file.open(filename);

	Material currentMaterial = Material("Red", Colour(255, 0, 0));

	while(getline(file, line)) {
		vector<string> splitline = split(line, ' ');
		if (splitline[0] == "v") {
			glm::vec3 point(stof(splitline[1]) * scale, stof(splitline[2]) * scale, stof(splitline[3]) * scale);
			points.push_back(point);
		} else if (splitline[0] == "vt") {
			glm::vec2 point(stof(splitline[1]), stof(splitline[2]));
			texturePoints.push_back(point);
		} else if (splitline[0] == "vn") {
			glm::vec3 point(stof(splitline[1]), stof(splitline[2]), stof(splitline[3]));
			normals.push_back(point);
		} else if (splitline[0] == "f") {
			vector<string> a = split(splitline[1], '/');
			vector<string> b = split(splitline[2], '/');
			vector<string> c = split(splitline[3], '/');

			ModelTriangle triangle = ModelTriangle(points[stoi(a[0]) - 1], points[stoi(b[0]) - 1], points[stoi(c[0]) - 1], currentMaterial.colour);

			if (currentMaterial.hasTexture) {
				TexturePoint pointA = TexturePoint(texturePoints[stoi(a[1]) - 1][0] * currentMaterial.texture.width, texturePoints[stoi(a[1]) - 1][1] * currentMaterial.texture.height);
				TexturePoint pointB = TexturePoint(texturePoints[stoi(b[1]) - 1][0] * currentMaterial.texture.width, texturePoints[stoi(b[1]) - 1][1] * currentMaterial.texture.height);
				TexturePoint pointC = TexturePoint(texturePoints[stoi(c[1]) - 1][0] * currentMaterial.texture.width, texturePoints[stoi(c[1]) - 1][1] * currentMaterial.texture.height);
				triangle.texturePoints = {{pointA, pointB, pointC}};
				triangle.textureMap = currentMaterial.texture;
			}

			if (a.size() == 3) {
				triangle.vertexNormals[0] = normals[stoi(a[2]) - 1];
				triangle.vertexNormals[1] = normals[stoi(b[2]) - 1];
				triangle.vertexNormals[1] = normals[stoi(c[2]) - 1];
			} else {
				triangle.vertexNormals[0] = glm::cross(triangle.vertices[0] - triangle.vertices[1], triangle.vertices[2] - triangle.vertices[1]);
				triangle.vertexNormals[1] = glm::cross(triangle.vertices[0] - triangle.vertices[1], triangle.vertices[2] - triangle.vertices[1]);
				triangle.vertexNormals[2] = glm::cross(triangle.vertices[0] - triangle.vertices[1], triangle.vertices[2] - triangle.vertices[1]);
			}

			triangle.normal = glm::cross(triangle.vertices[0] - triangle.vertices[1], triangle.vertices[2] - triangle.vertices[1]);

			triangles.push_back(triangle);
		} else if (splitline[0] == "usemtl") {
			string materialName = splitline[1];
			for (int i = 0; i < materials.size(); i++) {
				if (materials[i].name == materialName) {
					currentMaterial = materials[i];
					break;
				}
			}
		}
	}

	file.close();
}

void projectVertices(DrawingWindow &window, vector<ModelTriangle> triangles) {
	glm::vec3 camera(0.0, 0.0, 4.0);
	float focalLength = 2.0;

	for (int i = 0; i < triangles.size(); i++) {
		for (int j = 0; j < 3; j ++) {
			glm::vec3 vertex = triangles[i].vertices[j];
			float x = vertex[0] * -150;
			float y = vertex[1] * 150;
			float z = vertex[2] - 2.0;

			float u = (focalLength * (x / z)) + (WIDTH / 2);
			float v = (focalLength * (y / z)) + (HEIGHT / 2);

			uint32_t colourInt = (255 << 24) + (255 << 16) + (255 << 8) + 255;
			window.setPixelColour(u, v, colourInt);
		}
	}
}

void renderWireframes(DrawingWindow &window, vector<ModelTriangle> triangles, glm::vec3 camera, glm::mat3 cameraOrientation) {
	float focalLength = 2.0;
	float scale = 200.0;

	for (int i = 0; i < triangles.size(); i++) {
		vector<CanvasPoint> points;

		for (int j = 0; j < 3; j ++) {
			glm::vec3 vertex = triangles[i].vertices[j];

			glm::vec3 camToVertex(vertex[0] - camera[0], vertex[1] - camera[1], vertex[2] - camera[2]);
			glm::vec3 adjustedVector = camToVertex * cameraOrientation;

			float x = (-adjustedVector[0] * scale);
			float y = (adjustedVector[1] * scale);
			float z = (adjustedVector[2]);

			float u = (focalLength * (x / z)) + (WIDTH / 2.0);
			float v = (focalLength * (y / z)) + (HEIGHT / 2.0);

			points.push_back(CanvasPoint(u, v));
		}

		CanvasTriangle triangle = CanvasTriangle(points[0], points[1], points[2]);
		drawTriangle(window, triangle, triangles[i].colour);
	}
}

void renderFilledTriangles(DrawingWindow &window, vector<ModelTriangle> triangles) {
	glm::vec3 camera(0.0, 0.0, 4.0);
	float focalLength = 2.0;

	float depth[WIDTH + 1][HEIGHT + 1];

	for (int i = 0; i < triangles.size(); i++) {
		vector<CanvasPoint> points;

		for (int j = 0; j < 3; j ++) {
			glm::vec3 vertex = triangles[i].vertices[j];
			float x = (vertex[0] * -150) + camera[0];
			float y = (vertex[1] * 150) + camera[1];
			float z = (vertex[2] - (camera[2] / 2));

			float u = (focalLength * (x / z)) + (WIDTH / 2);
			float v = (focalLength * (y / z)) + (HEIGHT / 2);

			points.push_back(CanvasPoint(u, v, z));
		}

		CanvasTriangle triangle = CanvasTriangle(points[0], points[1], points[2]);
		drawFilledTriangle(window, triangle, triangles[i].colour, depth);
	}
}

void renderTexturedTriangles(DrawingWindow &window, vector<ModelTriangle> triangles, glm::vec3 camera, glm::mat3 cameraOrientation) {
	float focalLength = 2.0;
	float scale = 200.0;

	float depth[WIDTH + 1][HEIGHT + 1];
	for (int i = 0; i < WIDTH; i++) {
		for (int j = 0; j < HEIGHT; j++) {
			depth[i][j] = 0.0;
		}
	}

	for (int i = 0; i < triangles.size(); i++) {
		vector<CanvasPoint> points;

		for (int j = 0; j < 3; j ++) {
			glm::vec3 vertex = triangles[i].vertices[j];

			glm::vec3 camToVertex(vertex[0] - camera[0], vertex[1] - camera[1], vertex[2] - camera[2]);
			glm::vec3 adjustedVector = camToVertex * cameraOrientation;

			float x = (-adjustedVector[0] * scale);
			float y = (adjustedVector[1] * scale);
			float z = (adjustedVector[2]);

			float ua = (focalLength * (x / z)) + (WIDTH / 2.0);
			float va = (focalLength * (y / z)) + (HEIGHT / 2.0);
			CanvasPoint point = CanvasPoint(ua, va, z);

			if (triangles[i].texturePoints.size() != 0) {
				point.texturePoint = triangles[i].texturePoints[j];
			}
			points.push_back(point);
		}

		CanvasTriangle triangle = CanvasTriangle(points[0], points[1], points[2]);
		if (triangle.v0().texturePoint.x != 0) {
			drawTexturedTriangle(window, triangle, triangles[i].textureMap, depth);
		} else {
			drawFilledTriangle(window, triangle, triangles[i].colour, depth);
		}
	}
}

bool pointIsInTriangle(ModelTriangle triangle, float u, float v) {
	if (u >= 0.0 && u <= 1.0) {
		if (v >= 0.0 && v <= 1.0) {
			if (u + v <= 1.0) {
				return true;
			}
		}
	}

	return false;
}

glm::vec3 calculateWorldPos(glm::vec3 solution, ModelTriangle triangle) {
	glm::vec3 leftVec = triangle.vertices[1] - triangle.vertices[0];
	glm::vec3 rightVec = triangle.vertices[2] - triangle.vertices[0];

	glm::vec3 uVec = leftVec * solution[1];
	glm::vec3 vVec = rightVec * solution[2];

	glm::vec3 point = triangle.vertices[0] + uVec + vVec;

	return point;
}

RayTriangleIntersection getClosestIntersection(glm::vec3 camera, glm::vec3 direction, vector<ModelTriangle> triangles, int triangleIndex) {
	RayTriangleIntersection closestIntersection;

	closestIntersection.distanceFromCamera = -1;

	for (int i = 0; i < triangles.size(); i++) {
		if (i != triangleIndex) {
			glm::vec3 e0 = triangles[i].vertices[1] - triangles[i].vertices[0];
			glm::vec3 e1 = triangles[i].vertices[2] - triangles[i].vertices[0];
			glm::vec3 SPVector = camera - triangles[i].vertices[0];
			glm::mat3 DEMatrix(-direction, e0, e1);
			glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

			if (pointIsInTriangle(triangles[i], possibleSolution[1], possibleSolution[2])) {
				if ((-possibleSolution[0] < closestIntersection.distanceFromCamera || closestIntersection.distanceFromCamera == -1) && -possibleSolution[0] >= 0.0) {
					closestIntersection.distanceFromCamera = -possibleSolution[0];
					closestIntersection.intersectedTriangle = triangles[i];
					closestIntersection.triangleIndex = i;
					closestIntersection.solution = possibleSolution;

					closestIntersection.intersectionPoint = calculateWorldPos(possibleSolution, triangles[i]);
				}
			}
		}
	}

	return closestIntersection;
}

float round3DP(float number) {
	return round(number * 1000.0) / 1000.0;
}

void lookAt(glm::vec3 point, glm::vec3 camera, glm::mat3 &cameraOrientation) {
	glm::vec3 forward = glm::normalize(camera - point);
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0.0, 1.0, 0.0), forward));

	if (forward[2] < 0 && round3DP(forward[0]) == 0.0) {
		right[0] *= -1;
		right[1] *= -1;
		right[2] *= -1;
	}

	glm::vec3 up = glm::normalize(glm::cross(forward, right));

	cameraOrientation = glm::mat3(right, up, forward);
}

void rotateX(glm::vec3 &camera) {
	float theta = 2;
	glm::mat3 m(
		1, 0, 0,
		0, cos(theta*M_PI/180), sin(theta*M_PI/180),
		0, -sin(theta*M_PI/180), cos(theta*M_PI/180)
	);

	camera = m * camera;
}

void rotateY(glm::vec3 &camera) {
	float theta = 2;
	glm::mat3 m(
		cos(theta*M_PI/180), 0, -sin(theta*M_PI/180),
		0, 1, 0,
		sin(theta*M_PI/180), 0, cos(theta*M_PI/180)
	);

	camera = m * camera;
}

void tilt(glm::mat3 &cameraOrientation) {
	float theta = 30;
	glm::mat3 m(
		1, 0, 0,
		0, cos(theta*M_PI/180), sin(theta*M_PI/180),
		0, -sin(theta*M_PI/180), cos(theta*M_PI/180)
	);

	cameraOrientation = cameraOrientation * m;
}

void pan(glm::mat3 &cameraOrientation) {
	float theta = 30;
	glm::mat3 m(
		cos(theta*M_PI/180), 0, -sin(theta*M_PI/180),
		0, 1, 0,
		sin(theta*M_PI/180), 0, cos(theta*M_PI/180)
	);

	cameraOrientation = cameraOrientation * m;
}

float proximityLighting(glm::vec3 point, glm::vec3 light) {
	float distanceFromLight = glm::length(light - point);

	float brightness = 1.0 / (4 * M_PI * (distanceFromLight * distanceFromLight));

	if (brightness > 1.0) brightness = 1.0;

	return brightness;
}

float angleOfIncidenceLighting(glm::vec3 point, glm::vec3 light, ModelTriangle triangle) {
	glm::vec3 vectorToLight = light - point;

	float angleOfIncidence = glm::dot(triangle.normal, vectorToLight);

	if (angleOfIncidence < 0.0) {
		glm::vec3 normal = glm::cross(triangle.vertices[1] - triangle.vertices[0], triangle.vertices[2] - triangle.vertices[0]);
		angleOfIncidence = glm::dot(normal, vectorToLight);
	}

	if (angleOfIncidence > 1.0) angleOfIncidence = 1.0;
	if (angleOfIncidence < 0.0) angleOfIncidence = 0.0;

	return pow(angleOfIncidence, 0.4);
}

float specularLighting(glm::vec3 point, glm::vec3 light, ModelTriangle triangle, glm::vec3 camera) {
	glm::vec3 incidenceVector = point - light;

	glm::vec3 normal = triangle.normal;
	float angleOfIncidence = glm::dot(normal, incidenceVector);

	if (angleOfIncidence < 0.0) {
		normal = glm::cross(triangle.vertices[1] - triangle.vertices[0], triangle.vertices[2] - triangle.vertices[0]);
		angleOfIncidence = glm::dot(normal, incidenceVector);
	}

	glm::vec3 reflectionVector = incidenceVector - (2.0f * normal * glm::dot(incidenceVector, normal));
	glm::vec3 viewVector = glm::normalize(camera - point);

	float directionDifference = glm::dot(glm::normalize(reflectionVector), viewVector);

	directionDifference = pow(directionDifference, 256);

	return directionDifference;
}

bool inShadow(glm::vec3 point, glm::vec3 light, vector<ModelTriangle> triangles, int triangleIndex) {
	glm::vec3 direction;
	direction[0] = -(light[0] - point[0]);
	direction[1] = -(light[1] - point[1]);
	direction[2] = -(light[2] - point[2]);
	direction = glm::normalize(direction);

	float distanceFromLight = glm::length(light - point);

	RayTriangleIntersection inter = getClosestIntersection(point, direction, triangles, triangleIndex);

	if (inter.triangleIndex != 0 && inter.triangleIndex != 1 && triangles[triangleIndex].colour.name != "White") {
		if (distanceFromLight > inter.distanceFromCamera && inter.distanceFromCamera != -1) {
			return true;
		}
	}

	return false;
}

float ambientLighting(float brightness) {
	float ambientLight = 0.2;

	if (brightness < ambientLight) {
		return ambientLight;
	} else {
		return brightness;
	}
}

float getBrightness(glm::vec3 point, glm::vec3 light, vector<ModelTriangle> triangles, glm::vec3 camera, int triangleIndex) {
	float brightness = 0.0;

	if (!inShadow(point, light, triangles, triangleIndex)) {
		float proxLight = proximityLighting(point, light);
		float angleOfIncidence = angleOfIncidenceLighting(point, light, triangles[triangleIndex]);
		float specLighting = specularLighting(point, light, triangles[triangleIndex], camera);

		brightness = (proxLight + angleOfIncidence + specLighting) / 3;
	}

	brightness = ambientLighting(brightness);

	if (brightness > 1.0) brightness = 1.0;
	return brightness;
}

glm::vec3 calculateDirection(float u, float v, glm::vec3 camera, glm::mat3 cameraOrientation) {
	float focalLength = 2.0;
	float scale = 200;
	float x = (u - (WIDTH / 2.0)) / -scale;
	float y = (v - (HEIGHT / 2.0)) / scale;
	glm::vec3 cameraSpacePixel(x, y, -focalLength);
	glm::vec3 worldSpacePixel = (cameraSpacePixel * (cameraOrientation)) + camera;

	glm::vec3 direction = worldSpacePixel - camera;
	direction[2] *= -1;
	direction = glm::normalize(direction);

	return direction;
}

uint32_t getTexturePixelColour(ModelTriangle triangle, glm::vec3 solution) {
	float u = solution[1];
	float v = solution[2];

	glm::vec3 leftVec = triangle.vertices[1] - triangle.vertices[0];
	glm::vec3 rightVec = triangle.vertices[2] - triangle.vertices[0];

	glm::vec3 uVec = leftVec * u;
	glm::vec3 vVec = rightVec * v;

	glm::vec3 pointVec = uVec + vVec;

	float pointTextureX = triangle.texturePoints[0].x + ((triangle.texturePoints[1].x - triangle.texturePoints[0].x) * u) + ((triangle.texturePoints[2].x - triangle.texturePoints[0].x) * v);
	float pointTextureY = triangle.texturePoints[0].y + ((triangle.texturePoints[1].y - triangle.texturePoints[0].y) * u) + ((triangle.texturePoints[2].y - triangle.texturePoints[0].y) * v);

	uint32_t colourInt = triangle.textureMap.pixels[(round(pointTextureY) * triangle.textureMap.width) + round(pointTextureX)];

	return colourInt;
}

void drawRayTrace(DrawingWindow &window, vector<ModelTriangle> triangles, glm::vec3 camera, glm::mat3 cameraOrientation, glm::vec3 light) {
	window.clearPixels();

	for (float v = 0; v < HEIGHT; v++) {
		for (float u = 0; u < WIDTH; u++) {
			std::cout << "X: " << u << ", Y: " << v << "/" << HEIGHT << std::endl;
			glm::vec3 direction = calculateDirection(u, v, camera, cameraOrientation);

			RayTriangleIntersection inter = getClosestIntersection(camera, direction, triangles, -1);
			if (inter.distanceFromCamera == -1) {
				window.setPixelColour(round(u), round(v), 0);
			} else {
				float brightness = getBrightness(inter.intersectionPoint, light, triangles, camera, inter.triangleIndex);
				uint32_t colourInt;

				if (inter.intersectedTriangle.texturePoints.size() != 0) {
					colourInt = getTexturePixelColour(inter.intersectedTriangle, inter.solution);
				} else {
					Colour colour = inter.intersectedTriangle.colour;
					colourInt = (255 << 24) + (int(colour.red * brightness) << 16) + (int(colour.green * brightness) << 8) + int(colour.blue * brightness);
				}

				window.setPixelColour(round(u), round(v), colourInt);
			}
		}
	}
}

void draw(DrawingWindow &window, vector<ModelTriangle> triangles, glm::vec3 camera, glm::mat3 cameraOrientation, glm::vec3 light, int renderMode) {
	window.clearPixels();

	if (renderMode == 0) { // Raytrace
		drawRayTrace(window, triangles, camera, cameraOrientation, light);
		std::cout << "Rendered Ray" << std::endl;
	} else if (renderMode == 1) { // Raster
		renderTexturedTriangles(window, triangles, camera, cameraOrientation);
		std::cout << "Rendered Rast" << std::endl;
	} else if (renderMode == 2) { // Wireframe
		renderWireframes(window, triangles, camera, cameraOrientation);
		std::cout << "Rendered Wire" << std::endl;
	}

}

void update(DrawingWindow &window, glm::vec3 &camera, glm::mat3 &cameraOrientation) {
	rotateY(camera);
	lookAt(glm::vec3(50.0, 45.0, 0.0), camera, cameraOrientation);
}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &camera, glm::mat3 &cameraOrientation, int &renderMode, glm::vec3 &light) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) camera[0] += 0.3;
		else if (event.key.keysym.sym == SDLK_RIGHT) camera[0] -= 0.3;
	 	else if (event.key.keysym.sym == SDLK_UP) camera[1] -= 0.3;
	 	else if (event.key.keysym.sym == SDLK_DOWN) camera[1] += 0.3;
	 	else if (event.key.keysym.sym == SDLK_RCTRL) camera[2] += 0.3;
		else if (event.key.keysym.sym == SDLK_RSHIFT) camera[2] -= 0.3;
		else if (event.key.keysym.sym == SDLK_p) rotateX(camera);
		else if (event.key.keysym.sym == SDLK_o) rotateY(camera);
		else if (event.key.keysym.sym == SDLK_l) tilt(cameraOrientation);
		else if (event.key.keysym.sym == SDLK_k) pan(cameraOrientation);
		else if (event.key.keysym.sym == SDLK_i) lookAt(glm::vec3(50.0, 45.0, 0.0), camera, cameraOrientation);
		else if (event.key.keysym.sym == SDLK_r) renderMode = 0; // Raytrace
		else if (event.key.keysym.sym == SDLK_t) renderMode = 1; // Raster
		else if (event.key.keysym.sym == SDLK_y) renderMode = 2; // Wireframe
		else if (event.key.keysym.sym == SDLK_a) light[0] -= 1;
		else if (event.key.keysym.sym == SDLK_d) light[0] += 1;
		else if (event.key.keysym.sym == SDLK_w) light[1] += 1;
		else if (event.key.keysym.sym == SDLK_s) light[1] -= 1;
		else if (event.key.keysym.sym == SDLK_q) light[2] -= 1;
		else if (event.key.keysym.sym == SDLK_e) light[2] += 1;
		std::cout << renderMode << std::endl;
	}
	else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");

}

void outputFrame(int frame, DrawingWindow window) {
	string filename = "./output/frame" + to_string(frame) + ".ppm";
	window.savePPM(filename);
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	string foldername = "./logo/";

	vector<Material> materials;
	readMTL(materials, foldername + "materials.mtl", foldername);

	vector<ModelTriangle> triangles;
	readOBJ(triangles, foldername + "logo.obj", 0.17, materials);

	glm::vec3 light(400.0, 400.0, 40.0);
	glm::vec3 camera(50.0, 45.0, 150.0);
	glm::mat3 cameraOrientation(
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
	);

	int renderMode = 0;

	int frame = 0;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, camera, cameraOrientation, renderMode, light);
		update(window, camera, cameraOrientation);
		draw(window, triangles, camera, cameraOrientation, light, renderMode);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
		outputFrame(frame, window);
		frame++;
	}
}
