#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>

#define WIDTH 320
#define HEIGHT 240

void drawLine(DrawingWindow &window, CanvasPoint start, CanvasPoint end, Colour colour) {
	float xDiff = end.x - start.x;
	float yDiff = end.y - start.y;
	float numberOfSteps = fmaxf(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;
	for (float i = 0.0; i < numberOfSteps; i++) {
		float x = start.x + (xStepSize * i);
		float y = start.y + (yStepSize * i);
		uint32_t colourInt = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		window.setPixelColour(round(x), round(y), colourInt);
	}
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	drawLine(window, triangle.v0(), triangle.v1(), colour);
	drawLine(window, triangle.v1(), triangle.v2(), colour);
	drawLine(window, triangle.v2(), triangle.v0(), colour);
}

void sortVertices(CanvasPoint (&vertices)[3], CanvasPoint v0, CanvasPoint v1, CanvasPoint v2) {
	if (v0.y <= v1.y && v0.y <= v2.y) {
		vertices[0] = v0;
		if (v1.y <= v2.y) {
			vertices[1] = v1;
			vertices[2] = v2;
		} else {
			vertices[1] = v2;
			vertices[2] = v1;
		}
	} else if (v1.y <= v0.y && v1.y <= v2.y) {
		vertices[0] = v1;
		if (v0.y <= v2.y) {
			vertices[1] = v0;
			vertices[2] = v2;
		} else {
			vertices[1] = v2;
			vertices[2] = v0;
		}
	} else {
		vertices[0] = v2;
		if (v1.y <= v0.y) {
			vertices[1] = v1;
			vertices[2] = v0;
		} else {
			vertices[1] = v0;
			vertices[2] = v1;
		}
	}
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, float (&depth)[WIDTH + 1][HEIGHT + 1]) {
	uint32_t colourInt = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);

	CanvasPoint vertices [3];
	sortVertices(vertices, triangle.v0(), triangle.v1(), triangle.v2());

	float totalDiffX = vertices[2].x - vertices[0].x;
	float totalDiffY = vertices[2].y - vertices[0].y;
	float totalDepthDiff = vertices[2].depth - vertices[0].depth;
	float diffY = vertices[1].y - vertices[0].y;
	float pointX = vertices[0].x + (totalDiffX * (diffY / totalDiffY));
	float pointDepth = vertices[0].depth + (totalDepthDiff * (diffY / totalDiffY));
	CanvasPoint point = CanvasPoint(pointX, vertices[1].y, pointDepth);

	// Top triangle
	float leftDiffX = point.x - vertices[0].x;
	float leftDiffY = point.y - vertices[0].y;
	float leftDiffDepth = point.depth - vertices[0].depth;
	float rightDiffX = vertices[1].x - vertices[0].x;
	float rightDiffY = vertices[1].y - vertices[0].y;
	float rightDiffDepth = vertices[1].depth - vertices[0].depth;

	for (float y = vertices[0].y; y < point.y; y++) {
		diffY = y - vertices[0].y;
		int minX = vertices[0].x + (leftDiffX * (diffY / leftDiffY));
		int maxX = vertices[0].x + (rightDiffX * (diffY / rightDiffY));

		float leftDepth = vertices[0].depth + (leftDiffDepth * (diffY / leftDiffY));
		float rightDepth = vertices[0].depth + (rightDiffDepth * (diffY / rightDiffY));
		float depthDiff = rightDepth - leftDepth;

		if (minX > maxX) {
			int a = minX;
			minX = maxX;
			maxX = a;
		}

		for (int x = minX; x <= maxX; x++) {
			if (x > 0 && y > 0 && x < WIDTH && y < HEIGHT) {
				float z = leftDepth + (depthDiff * (x / maxX));
				z = (1/z) * -1;

				if (z > depth[x][int(round(y))]) {
					depth[x][int(round(y))] = z;
					window.setPixelColour(round(x), round(y), colourInt);
				}
			}
		}
	}

	// Bottom Triangle
	leftDiffX = vertices[2].x - point.x;
	leftDiffY = vertices[2].y - point.y;
	leftDiffDepth = vertices[2].depth - point.depth;
	rightDiffX = vertices[2].x - vertices[1].x;
	rightDiffY = vertices[2].y - vertices[1].y;
	rightDiffDepth = vertices[2].depth - vertices[1].depth;

	for (float y = vertices[1].y; y < vertices[2].y; y++) {
		diffY = y - vertices[1].y;
		int minX = point.x + (leftDiffX * (diffY / leftDiffY));
		int maxX = vertices[1].x + (rightDiffX * (diffY / rightDiffY));

		float leftDepth = point.depth + (leftDiffDepth * (diffY / leftDiffY));
		float rightDepth = point.depth + (rightDiffDepth * (diffY / rightDiffY));
		float depthDiff = rightDepth - leftDepth;

		if (minX > maxX) {
			int a = minX;
			minX = maxX;
			maxX = a;
		}

		for (int x = minX; x <= maxX; x++) {
			if (x > 0 && y > 0 && x < WIDTH && y < HEIGHT) {
				float z = leftDepth + (depthDiff * (x / maxX));
				z = (1/z) * -1;

				if (z > depth[x][int(round(y))]) {
					depth[x][int(round(y))] = z;
					window.setPixelColour(round(x), round(y), colourInt);
				}
			}
		}
	}

}

void drawTexturedTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap textureMap, float (&depth)[WIDTH + 1][HEIGHT + 1]) {
	CanvasPoint vertices [3];
	sortVertices(vertices, triangle.v0(), triangle.v1(), triangle.v2());

	float totalDiffX = vertices[2].x - vertices[0].x;
	float totalDiffY = vertices[2].y - vertices[0].y;
	float totalDepthDiff = vertices[2].depth - vertices[0].depth;
	float diffY = vertices[1].y - vertices[0].y;
	float pointX = vertices[0].x + (totalDiffX * (diffY / totalDiffY));
	float pointDepth = vertices[0].depth + (totalDepthDiff * (diffY / totalDiffY));

	CanvasPoint point = CanvasPoint(pointX, vertices[1].y, pointDepth);

	float xAmount = ((pointX - vertices[0].x) / totalDiffX);
	float textureDiffX = vertices[2].texturePoint.x - vertices[0].texturePoint.x;
	float texturePointX = vertices[0].texturePoint.x + (textureDiffX * xAmount);

	float yAmount = (vertices[1].y - vertices[0].y) / (vertices[2].y - vertices[0].y);
	float textureDiffY = vertices[2].texturePoint.y - vertices[0].texturePoint.y;
	float texturePointY = vertices[0].texturePoint.y + (textureDiffY * yAmount);

	point.texturePoint = TexturePoint(texturePointX, texturePointY);

	CanvasPoint leftPoint;
	CanvasPoint rightPoint;

	if (point.x < vertices[1].x) {
		leftPoint = point;
		rightPoint = vertices[1];
	} else {
		leftPoint = vertices[1];
		rightPoint = point;
	}

	// Top triangle
	float leftDiffX = point.x - vertices[0].x;
	float leftDiffY = point.y - vertices[0].y;
	float leftDiffDepth = point.depth - vertices[0].depth;
	float rightDiffX = vertices[1].x - vertices[0].x;
	float rightDiffY = vertices[1].y - vertices[0].y;
	float rightDiffDepth = vertices[1].depth - vertices[0].depth;

	for (float y = vertices[0].y; y < point.y; y++) {
		diffY = y - vertices[0].y;
		float minX = vertices[0].x + (leftDiffX * (diffY / leftDiffY));
		float maxX = vertices[0].x + (rightDiffX * (diffY / rightDiffY));

		float leftDepth = vertices[0].depth + (leftDiffDepth * (diffY / leftDiffY));
		float rightDepth = vertices[0].depth + (rightDiffDepth * (diffY / rightDiffY));
		float depthDiff = rightDepth - leftDepth;

		if (minX > maxX) {
			float a = minX;
			minX = maxX;
			maxX = a;
		}

		float xPercentageLeft = (minX - vertices[0].x) / (leftPoint.x - vertices[0].x) ;
		float leftTextureDiffX = leftPoint.texturePoint.x - vertices[0].texturePoint.x;
		float leftTextureX = vertices[0].texturePoint.x + (leftTextureDiffX * xPercentageLeft);

		float yPercentageLeft = (y - vertices[0].y) / (leftPoint.y - vertices[0].y);
		float leftTextureDiffY = leftPoint.texturePoint.y - vertices[0].texturePoint.y;
		float leftTextureY = vertices[0].texturePoint.y + (leftTextureDiffY * yPercentageLeft);

		float xPercentageRight = (maxX - vertices[0].x) / (rightPoint.x - vertices[0].x);
		float rightTextureDiffX = rightPoint.texturePoint.x - vertices[0].texturePoint.x;
		float rightTextureX = vertices[0].texturePoint.x + (rightTextureDiffX * xPercentageRight);

		float yPercentageRight = (y - vertices[0].y) / (rightPoint.y - vertices[0].y);
		float rightTextureDiffY = rightPoint.texturePoint.y - vertices[0].texturePoint.y;
		float rightTextureY = vertices[0].texturePoint.y + (rightTextureDiffY * yPercentageRight);

		for (float x = minX; x <= maxX; x++) {
			if (x > 0 && y > 0 && x < WIDTH && y < HEIGHT) {
				float z = leftDepth + (depthDiff * (x / maxX));
				z = (1/z) * -1;

				// Calculate texture pixel
				float xPercentage = (x - minX) / (maxX - minX);

				float xDiff = rightTextureX - leftTextureX;
				float texturePointX = leftTextureX + (xDiff * xPercentage);

				float yDiff = rightTextureY - leftTextureY;
				float texturePointY = leftTextureY + (yDiff * xPercentage);

				if (textureMap.pixels.size() > (round(texturePointY) * textureMap.width) + round(texturePointX)) {
					uint32_t colourInt = textureMap.pixels[(round(texturePointY) * textureMap.width) + round(texturePointX)];

					if (z > depth[int(round(x))][int(round(y))]) {
						depth[int(round(x))][int(round(y))] = z;
						window.setPixelColour(round(x), round(y), colourInt);
					}
				}
			}
		}
	}

	// Bottom Triangle
	leftDiffX = vertices[2].x - point.x;
	leftDiffY = vertices[2].y - point.y;
	leftDiffDepth = vertices[2].depth - point.depth;
	rightDiffX = vertices[2].x - vertices[1].x;
	rightDiffY = vertices[2].y - vertices[1].y;
	rightDiffDepth = vertices[2].depth - vertices[1].depth;

	for (float y = vertices[1].y; y < vertices[2].y; y++) {
		diffY = y - vertices[1].y;
		float minX = point.x + (leftDiffX * (diffY / leftDiffY));
		float maxX = vertices[1].x + (rightDiffX * (diffY / rightDiffY));

		float leftDepth = point.depth + (leftDiffDepth * (diffY / leftDiffY));
		float rightDepth = point.depth + (rightDiffDepth * (diffY / rightDiffY));
		float depthDiff = rightDepth - leftDepth;

		if (minX > maxX) {
			float a = minX;
			minX = maxX;
			maxX = a;
		}

		float xPercentageLeft = (minX - leftPoint.x) / (vertices[2].x - leftPoint.x);
		float leftTextureDiffX = vertices[2].texturePoint.x - leftPoint.texturePoint.x;
		float leftTextureX = leftPoint.texturePoint.x + (leftTextureDiffX * xPercentageLeft);

		float yPercentageLeft = (y - leftPoint.y) / (vertices[2].y - leftPoint.y);
		float leftTextureDiffY = vertices[2].texturePoint.y - leftPoint.texturePoint.y;
		float leftTextureY = leftPoint.texturePoint.y + (leftTextureDiffY * yPercentageLeft);

		float xPercentageRight = (maxX - rightPoint.x) / (vertices[2].x - rightPoint.x);
		float rightTextureDiffX = vertices[2].texturePoint.x - rightPoint.texturePoint.x;
		float rightTextureX = rightPoint.texturePoint.x + (rightTextureDiffX * xPercentageRight);

		float yPercentageRight = (y - rightPoint.y) / (vertices[2].y - rightPoint.y);
		float rightTextureDiffY = vertices[2].texturePoint.y - rightPoint.texturePoint.y;
		float rightTextureY = rightPoint.texturePoint.y + (rightTextureDiffY * yPercentageRight);

		for (float x = minX; x <= maxX; x++) {
			if (x > 0 && y > 0 && x < WIDTH && y < HEIGHT) {
				float z = leftDepth + (depthDiff * (x / maxX));
				z = (1/z) * -1;

				// Calculate texture pixel
				float xPercentage = (x - minX) / (maxX - minX);

				float xDiff = rightTextureX - leftTextureX;
				float texturePointX = leftTextureX + (xDiff * xPercentage);

				float yDiff = rightTextureY - leftTextureY;
				float texturePointY = leftTextureY + (yDiff * xPercentage);

				if (textureMap.pixels.size() > (round(texturePointY) * textureMap.width) + round(texturePointX)) {
					uint32_t colourInt = textureMap.pixels[(round(texturePointY) * textureMap.width) + round(texturePointX)];

					if (z > depth[int(round(x))][int(round(y))]) {
						depth[int(round(x))][int(round(y))] = z;
						window.setPixelColour(round(x), round(y), colourInt);
					}
				}
			}
		}
	}
}
