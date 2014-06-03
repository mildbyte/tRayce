#include <iostream>

#include "Scene.h"
#include "Sphere.h"
#include "Box.h"
#include "Plane.h"
#include "Triangle.h"
#include <ctime>

using namespace std;

//Generates a random scene with a bunch of spheres
void randomScene(Scene* scene) {
    for (int i = 0; i < 30; i++) {
        Sphere *testSphere = new Sphere(
            (double)(rand() % 30 / 15.0f),
            Vector( (rand() % 240 - 120)/10.0, (rand() % 160 - 80) / 10.0, (rand() % 100) / 10.0 )
        );

        testSphere->material.color = Vector((rand()% 256 / 256.0), (rand()% 256 / 256.0), (rand()% 256 / 256.0));
        testSphere->material.isReflective = false;
        testSphere->material.reflectivity = 0.1;
        testSphere->material.diffuse = 1.0;
        scene->addRenderable(testSphere);
    }
}

//A sample postprocess routine. This one makes the image grayscale.
BitmapPixel postProcess(BitmapPixel pix) {
	double lum = pix.color.getX() * 0.30 + pix.color.getY() * 0.59 + pix.color.getZ() * 0.11;
    pix.color.setX (lum);
    pix.color.setY (lum);
    pix.color.setZ (lum);        
    return pix;
}

BitmapPixel amplify(BitmapPixel pix) {
	BitmapPixel result;
	result.color.set(pix.color.getX() * 100, pix.color.getY() * 100, pix.color.getZ() * 100);
	return result;
}

void addQuad(Scene &s, Vector origin, Vector s1, Vector s2, Material m) {
    Triangle *t1 = new Triangle(origin, origin + s1 + s2, origin + s1);
    t1->material = m;
    
    Triangle *t2 = new Triangle(origin, origin + s2, origin + s1 + s2);
    t2 -> material = m;
    
    s.addRenderable(t1);
    s.addRenderable(t2);
}


int main()
{
    srand((unsigned)time(0));
    Scene scene(640, 480);
    //scene.camera.width = 32;
    scene.camera.width = 16;
    scene.camera.height = 12;
    scene.camera.position.setX(1.5);
    scene.backgroundColor.set(0, 0, 0);
    scene.msaaSamples = 2;
    scene.msaaOptimize = false;
    scene.softShadowSamples = 1;
    scene.traceDepth = 1;
    scene.camera.position.setZ(-20);
    scene.camera.planeDistance = 15;
    scene.camera.lensRadius = 0;
    scene.camera.focalDistance = 25;

    scene.renderingMode = PATHTRACING;
    scene.pathTracingSamplesPerPixel = 32; //spp squared is actually cast
    scene.pathTracingMaxDepth = 7; // Too few samples and rays that go through
    // a sphere, bounce off a wall, through the sphere again and to the light
    // will terminate too early
    
    scene.doFinalGather = true;
    scene.photonCount = 65536;
    scene.photonBounces = 5;
    scene.photonGatherAmount = 32;
    scene.photonGatherSamples = 16;
    scene.irradiancePhotonFrequency = 4;
    scene.photonGatherDotThreshold = 0.9;
    scene.doIrradianceCaching = true;
    scene.samplingMode = STRATIFIED;
    
    // Problem with this scene: when we do a final gather in the corner
    // shaded by the big sphere, most rays hit the sphere. We then look up
    // an irradiance photon with a similar normal, but the closest one we find
    // is one on the emitting sphere. Hence we do light the corner that's supposed
    // to be shaded. Temporary solution: more photons so that at least some get
    // to the corner. Better solution?

    Sphere *redSphere = new Sphere(5, Vector(-4, 5, 13));
    Sphere *greenSphere = new Sphere(3, Vector(8, 7, 7));
    /*
    Sphere *blueLight = new Sphere(3, Vector(-1.5, -10.5, 6.5));
    Sphere *redLight = new Sphere(3, Vector(4.5, -10.5, 6.5));
    Sphere *greenLight = new Sphere(3, Vector(1.5, -10.5, 10.5));
    */
    
    Sphere *oneLight = new Sphere(100, Vector(1.5, -109.8, 8.5));

    Plane *bottomPlane = new Plane(Vector(0, 10, 0), Vector(0, -1, 0));
    Plane *upPlane = new Plane(Vector(0, 0, 18), Vector(0, 0, -1));
    Plane *leftPlane = new Plane(Vector(-11, 0, 0), Vector(1, 0, 0));
    Plane *rightPlane = new Plane(Vector(14, 0, 0), Vector(-1, 0, 0));
    Plane *topPlane = new Plane(Vector(0, -10, 0), Vector(0, 1, 0));
    Plane *backPlane = new Plane(Vector(0, 0, -10), Vector(0, 0, 1));
    
    oneLight->material.color = Vector(.95, .95, .95);
    oneLight->material.emittance = Vector(30, 30, 30);

    /*
    blueLight->material.color = Vector(.35,.35,.95);
    blueLight->material.emittance.set(30,20,20);
    redLight->material.color = Vector(.95,.35,.35);
    redLight->material.emittance.set(30,20,20);
    greenLight->material.color = Vector(.35,.95,.35);
    greenLight->material.emittance.set(30,30,20);
    */
    
    greenSphere->material.color.set(0.35, 0.95, 0.35);
    greenSphere->material.refrIndex = 1.42;
    greenSphere->material.transparency = 0.9;
    greenSphere->material.isTransparent = true;
    greenSphere->material.reflectivity = 0.1;
    greenSphere->material.isReflective = true;
    
    redSphere->material.color = Vector(.95,.35,.35);
    redSphere->material.refrIndex = 1.42;
    redSphere->material.transparency = 0.75;
    redSphere->material.isTransparent = true;
    redSphere->material.reflectivity = 0.25;
    redSphere->material.isReflective = true;

    bottomPlane->material.color = Vector(.95, .95, .95);
    upPlane->material.color = Vector(.95, .95, .95);
    leftPlane->material.color = Vector(.95, .95, .95);
    rightPlane->material.color = Vector(.95, .95, .95);
    topPlane->material.color = Vector(.95, .95, .95);
    backPlane->material.color = Vector(.95, .95, .95);
    
    Light *pointLight = new PointLight();
    pointLight->position = Vector(1.5, 16, 6);
    pointLight->brightness = 10;
    scene.addLight(pointLight);
    
        
    Material m;
    m.isTransparent = true;
    m.transparency = 0.95;
    m.reflectivity = 0.05;
    m.isReflective = true;
    m.refrIndex = 1.42;
    m.color = Vector(0, 1, 0);
    
    Triangle *t1 = new Triangle(Vector(4, 10, 13), Vector(8, -3, 9), Vector(12, 10, 13));
    Triangle *t2 = new Triangle(Vector(12, 10, 13), Vector(8, -3, 9), Vector(8, 10, 5));
    Triangle *t3 = new Triangle(Vector(8, 10, 5), Vector(4, 10, 13), Vector(8, -3, 9));
    
    t1->material = m;
    t2->material = m;
    t3->material = m;
    
    //addQuad(scene, Vector(5, 10, 17), Vector(0, -7, -7), Vector(-7, 0, 0), m);
    /*
    scene.addRenderable(redLight);
    scene.addRenderable(greenLight);
    scene.addRenderable(blueLight);
    */
    
    scene.addRenderable(oneLight);
    
    scene.addRenderable(t1);
    scene.addRenderable(t2);
    scene.addRenderable(t3);
    
    scene.addRenderable(bottomPlane);
    scene.addRenderable(upPlane);
    scene.addRenderable(leftPlane);
    scene.addRenderable(rightPlane);
    scene.addRenderable(topPlane);
    scene.addRenderable(backPlane);
    scene.ambientCoefficient = 5;
    //Check if the precalculated map exists and is valid
    //bool mapExists = scene.loadMap("map.dat");

    scene.render("test.bmp", NULL, 8);

    //Save the calculated map for future use
    //if (!mapExists) scene.saveMap("map.dat");

    return 0;
}
