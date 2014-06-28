//Data about the scene, the rendering routine, object list
#ifndef SCENE_H
#define SCENE_H

#include "Ray.h"
#include "RenderablesList.h"
#include "Light.h"
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
    std::list<Light*> lights_;
    RenderablesList renderables_;
    
    std::vector<Triangle*> trianglesVector_;
    KDNode* triangles_;

    //The photon map
    PhotonMap* photonMap_;

    //The rendered buffer
    Bitmap* rendered_;

    //Traces a single ray and returns the resultant color
    Vector traceRay(const Ray ray, int depth, double &dist);

    int width_, height_; //Horizontal and vertical resolution

    Vector xPixel, yPixel; //X and Y directions in the image plane

    //Shading coefficient (how much a point is obscured) calculation
    double calculateShadingCoefficient(Light* light, Vector point, Vector toLight, double lightDist);

    //Phong (diffuse+specular) shading
    Vector calculatePhongColor(Intersection inter, Ray ray);

    //Transparency contribution
    Vector calculateRefraction(Intersection inter, Ray ray, int depth);
    void refractRay(Intersection inter, Ray &ray, double &reflectance);

    //Reflection contribution
    Vector calculateReflection(Intersection inter, Ray ray, int depth);
    Ray reflectRay(Intersection inter, Ray ray);

	//Returns the closest intersection of the ray from the triangles' or other renderables' list
	Intersection getClosestIntersection(Ray ray, double cullDistance);

    //Convert one pixel to world coordinates and trace it
	//dist: distance to the first intersection (INFINITY if no intersections)
    Vector tracePixel(double x, double y, double &dist);
    
    //Path tracing
    Vector pathTrace(const Ray ray, int depth, double &dist);

	Vector getColorAt(Vector point);

    //Throw photons at the wall and see what sticks
    void populatePhotonMap();

    //Gets the irradiance value at a point with a normal and two coordinages on the
    //hemisphere of samples
    Vector sampleMapAt(Vector coords, Vector normal, double x, double y);

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
    void addLight(Light* light);
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

    //Ambient lighting strength
    double ambientCoefficient;

    //The maximum tracing depth
    int traceDepth;

    //Soft shadow samples (1 means treat soft shadow lights as ordinary lights)
    int softShadowSamples;

    Vector backgroundColor;

    bool doFinalGather;//Shoot rays to nearby objects or use local irradiance data?
    int photonCount;   //How many primary photons to launch?
    int photonBounces; //Maximum number of photon bounces
    int photonGatherAmount; //How many photons to use in irradiance calculations?
    int photonGatherSamples; //How many samples for final gathering step?
    int irradiancePhotonFrequency; //every nth photon becomes an irradiance photon
    double photonGatherDotThreshold; //How close should the two normals be for an
    //irradiance photon to be used as an estimate? (1.0 is same direction, -1.0 is opposite)
    bool doIrradianceCaching;
    
    int pathTracingSamplesPerPixel;
    int pathTracingMaxDepth;

    SamplingMode samplingMode; //Stratified or Halton?
    RenderingMode renderingMode; //Raytracing, path tracing or photon mapping

    bool visualizePhotons; //if true, only shows the positions of photons.
};

#endif
