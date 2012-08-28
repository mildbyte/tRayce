//Data about the scene, the rendering routine, object list
#ifndef SCENE_H
#define SCENE_H

#include "Ray.h"
#include "RenderablesList.h"
#include "Light.h"
#include "Bitmap.h"
#include "Camera.h"
#include <list>
#include <cstring>
#include <cstdlib>

class Scene {
private:
    //Scene objects
    std::list<Light*> lights_;
    RenderablesList renderables_;

    //The rendered buffer
    Bitmap* rendered_;

    //Traces a single ray and returns the resultant color
    Vector traceRay(const Ray ray, int depth);

    double width_, height_; //Horizontal and vertical resolution

    Renderable* prevHit_; //Previously hit object, used to optimize AA (only
                          //trace several rays when encountering a new object)

    double prevDist_;      //Distance to the hit object

    //Store the previously traced row of objects and the current
    //Used to optimize AA (do AA only if the pixel above it has a different
    //object in it)
    Renderable **prevRow_, **currRow_;

    //Shading coefficient (how much a point is obscured) calculation
    double calculateShadingCoefficient(Light* light, Vector point, Vector toLight, double lightDist);

    //Phong (diffuse+specular) shading
    Vector calculatePhongColor(Intersection inter, Ray ray);

    //Transparency contribution
    Vector calculateRefraction(Intersection inter, Ray ray, int depth);

    //Reflection contribution
    Vector calculateReflection(Intersection inter, Ray ray, int depth);
public:
    //Necessary to specify width and height to allocate memory
    Scene(int width, int height);

    //Add objects to the scene
    void addRenderable(Renderable* renderable);
    void addLight(Light* light);

    //Render the scene to a file
    void render(char* filename, BitmapPixel (*postProcess)(BitmapPixel));

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
};

#endif
