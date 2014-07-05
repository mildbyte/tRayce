//Data about the scene, the rendering routine, object list
#ifndef SCENE_H
#define SCENE_H

#include "Ray.h"
#include "RenderablesList.h"
#include "Bitmap.h"
#include "Camera.h"
#include "PhotonMap.h"
#include "Random.h"
#include "KDNode.h"
#include "Triangle.h"
#include <list>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <limits>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


enum SamplingMode {
    STRATIFIED,
    HALTON
};

enum RenderingMode {
    RAYTRACING,
    PHOTONMAPPING,
    PATHTRACING
};

class Scene {
private:
    //Scene objects
    RenderablesList renderables_;
    std::vector<Triangle*> trianglesVector_;
    KDNode* triangles_;

    //The rendered buffer
    Bitmap* rendered_;

    int width_, height_; //Horizontal and vertical resolution

    Vector xPixel, yPixel; //X and Y directions in the image plane

    void refractRay(Intersection inter, Ray &ray, double &reflectance);
    Ray reflectRay(Intersection inter, Ray ray);

	//Returns the closest intersection of the ray from the triangles' or other renderables' list
	Intersection getClosestIntersection(Ray ray, double cullDistance);

    //Convert one pixel to world coordinates and trace it
	//dist: distance to the first intersection (INFINITY if no intersections)
    Vector tracePixel(double x, double y, double &dist);
    
    //Path tracing
    Vector pathTrace(const Ray ray, int depth, double &dist);

	Vector getColorAt(Vector point);

    int pixelsRendered_;

    double *haltonXCoords_, *haltonYCoords_;
	
	//Directions in the image plane
	Vector xPixel_, yPixel_;
    
    Vector imagePlaneX_, imagePlaneY_;
	
	std::thread *renderingThreads_;
	
	void threadDoWork(int threadId, int noThreads);
	
public:
    //Necessary to specify width and height to allocate memory
    Scene(int width, int height);

    //Add objects to the scene
    void addRenderable(Renderable* renderable);
    void addTriangle(Triangle* triangle);

    //Render the scene to a file
    void render(char* filename, BitmapPixel (*postProcess)(BitmapPixel), int noThreads);

	//Import an OBJ file into the scene
	void importObj(char* filename, Material m, Vector shift, double scale);

    //Load the photon map from a file (generates the map from scratch otherwise)
    bool loadMap(char* path);

    //Save the photon map to a file
    void saveMap(char* path);

    //Data about the image plane
    Camera camera;

    //The maximum tracing depth
    int traceDepth;

    Vector backgroundColor;

    int pathTracingSamplesPerPixel;
    int pathTracingMaxDepth;
	double pathTracingTerminationProbability;

    SamplingMode samplingMode; //Stratified or Halton?
};

#endif
