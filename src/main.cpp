#include <iostream>

#include "Scene.h"
#include "Sphere.h"
#include "Box.h"
#include "Plane.h"
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

int main()
{
    srand((unsigned)time(0));
    Scene scene(640, 480);
    //scene.camera.width = 32;
    scene.camera.width = 16;
    scene.camera.height = 12;
    scene.backgroundColor.set(0, 0, 0);
    //scene.doAA = true;
    scene.msaaSamples = 2;
    scene.msaaOptimize = false;
    scene.softShadowSamples = 1;
    scene.traceDepth = 1;
    scene.camera.position.setZ(-10);
    scene.camera.planeDistance = 10;

    scene.renderingMode = PHOTONMAPPING;
    scene.pathTracingSamplesPerPixel = 10;
    scene.pathTracingMaxDepth = 5;
    
    scene.doFinalGather = true;
    //scene.visualizePhotons = true;
    scene.photonCount = 65536;
    scene.photonBounces = 5;
    scene.photonGatherAmount = 32;
    scene.photonGatherSamples = 16;
    scene.irradiancePhotonFrequency = 4;
    scene.photonGatherDotThreshold = 0.9;
    scene.samplingMode = STRATIFIED;
    
/*    Sphere *sphere1 = new Sphere(1000, Vector(0, 1000, 0));
    sphere1->material.color = Vector(0.5, 0.8, 0.5);
    Sphere *sphere2 = new Sphere(1000, Vector(0, 0, 1020));
    sphere2->material.color = Vector(0.8, 0.5, 0.5);
    Sphere *sphere3 = new Sphere(5, Vector(0, -7, 10));
    sphere3->material.color = Vector(1, 1, 1);
    sphere3->material.emittance = Vector(5000, 5000, 5000);
    
    scene.addRenderable(sphere1);
    scene.addRenderable(sphere2);
    scene.addRenderable(sphere3);
  */  
    
    // Problem with this scene: when we do a final gather in the corner
    // shaded by the big sphere, most rays hit the sphere. We then look up
    // an irradiance photon with a similar normal, but the closest one we find
    // is one on the emitting sphere. Hence we do light the corner that's supposed
    // to be shaded. Temporary solution: more photons so that at least some get
    // to the corner. Better solution?

    Sphere *redSphere = new Sphere(5, Vector(-1, 5, 13));
    Sphere *blueSphere = new Sphere(3, Vector(10, 0, 15));

//    Box *redSphere = new Box(Vector(-6, -2, 9), Vector(6, 6, 6));
//    Box *greenSphere = new Box(Vector(1, -2, 8), Vector(6, 6, 6));
//    Box *blueSphere = new Box(Vector(-6.5, 4, 7.5), Vector(6, 6, 6));
//    Box *yellowSphere = new Box(Vector(0.5, 4, 7), Vector(6, 6, 6));

    Plane *bottomPlane = new Plane(Vector(0, 10, 0), Vector(0, -1, 0));
    Plane *upPlane = new Plane(Vector(0, 0, 18), Vector(0, 0, -1));
    Plane *leftPlane = new Plane(Vector(-11, 0, 0), Vector(1, 0, 0));
    Plane *rightPlane = new Plane(Vector(14, 0, 0), Vector(-1, 0, 0));
    Plane *topPlane = new Plane(Vector(0, -10, 0), Vector(0, 1, 0));
    Plane *backPlane = new Plane(Vector(0, 0, -5), Vector(0, 0, 1));

    blueSphere->material.color = Vector(.35,.35,.75);
    blueSphere->material.emittance.set(2000,2000,4000);
    
    redSphere->material.color = Vector(.75,.35,.35);
    redSphere->material.emittance = Vector(0, 0, 0);

    bottomPlane->material.color = Vector(.95, .95, .95);
    upPlane->material.color = Vector(.25, .75, .25);
    leftPlane->material.color = Vector(.75, .25, .25);
    rightPlane->material.color = Vector(.25, .25, .75);
    topPlane->material.color = Vector(.75, .75, .75);
    backPlane->material.color = Vector(.75, .75, .75);


    scene.addRenderable(redSphere);
    scene.addRenderable(blueSphere);
    scene.addRenderable(bottomPlane);
    scene.addRenderable(upPlane);
    scene.addRenderable(leftPlane);
    scene.addRenderable(rightPlane);
    scene.addRenderable(topPlane);
    scene.addRenderable(backPlane);
    scene.ambientCoefficient = 0;
/*
    AreaLight* topLight = new AreaLight();
    topLight->position.set(1.5, -8, 8);
    topLight->dir1.set(1, 0, 0);
    topLight->dir2.set(0, 0, 1);
    topLight->size1 = 3;
    topLight->size2 = 3;
    topLight->brightness = 20000;
    scene.addLight(topLight);
*/    
    //Check if the precalculated map exists and is valid
    bool mapExists = scene.loadMap("map.dat");

    scene.render("test.bmp", NULL, 8);

    //Save the calculated map for future use
    if (!mapExists) scene.saveMap("map.dat");

    return 0;
}
