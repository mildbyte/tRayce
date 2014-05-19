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
#include <list>
#include <cstring>
#include <cstdlib>
#include <thread>

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

    //The photon map
    PhotonMap* photonMap_;

    //The rendered buffer
    Bitmap* rendered_;

    //Traces a single ray and returns the resultant color
    Vector traceRay(const Ray ray, int depth);

    int width_, height_; //Horizontal and vertical resolution

    Renderable* prevHit_; //Previously hit object, used to optimize AA (only
                          //trace several rays when encountering a new object)

    double prevDist_;      //Distance to the hit object

    //Store the previously traced row of objects and the current
    //Used to optimize AA (do AA only if the pixel above it has a different
    //object in it)
    Renderable **prevRow_, **currRow_;

    Vector xPixel, yPixel; //X and Y directions in the image plane

    //Shading coefficient (how much a point is obscured) calculation
    double calculateShadingCoefficient(Light* light, Vector point, Vector toLight, double lightDist);

    //Slightly shift the origin of the ray so that it doesn't collide with the previous object
    Ray epsilonShift(Ray ray);

    //Phong (diffuse+specular) shading
    Vector calculatePhongColor(Intersection inter, Ray ray);

    //Transparency contribution
    Vector calculateRefraction(Intersection inter, Ray ray, int depth);
    bool refractRay(Intersection inter, Ray& ray);

    //Reflection contribution
    Vector calculateReflection(Intersection inter, Ray ray, int depth);
    Ray reflectRay(Intersection inter, Ray ray);

    //Convert one pixel to world coordinates and trace it
    Vector tracePixel(double x, double y);
    
    //Path tracing
    Vector pathTrace(const Ray ray, int depth);

    //Throw photons at the wall and see what sticks
    void populatePhotonMap();

    //Gets the irradiance value at a point with a normal and two coordinages on the
    //hemisphere of samples
    Vector sampleMapAt(Vector coords, Vector normal, double x, double y);

    int pixelsRendered_;

    double *haltonXCoords_, *haltonYCoords_;
	
	//Directions in the image plane
	Vector xPixel_, yPixel_;
	
	std::thread *renderingThreads_;
	
	void threadDoWork(int threadId, int noThreads);
	
public:
    //Necessary to specify width and height to allocate memory
    Scene(int width, int height);

    //Add objects to the scene
    void addRenderable(Renderable* renderable);
    void addLight(Light* light);

    //Render the scene to a file
    void render(char* filename, BitmapPixel (*postProcess)(BitmapPixel), int noThreads);

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

    //Do antialiasing?
    bool doAA;

    //Do AA only when encountering a border?
    bool msaaOptimize;

    //The amount of antialiasing samples
    int msaaSamples;

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
    
    //If the intensity variance of the pixel drops below the cull value before exceeding the max
    //sample size, stop the sampling
    double pathTracingVarianceCull;
    int pathTracingMinBeforeCull;

    SamplingMode samplingMode; //Stratified or Halton?
    RenderingMode renderingMode; //Raytracing, path tracing or photon mapping

    bool visualizePhotons; //if true, only shows the positions of photons.
};

#endif
