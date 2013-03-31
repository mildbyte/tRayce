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

int main()
{
    srand((unsigned)time(0));
    Scene scene(160, 120);
//    scene.camera.height = 15;
    scene.backgroundColor.set(0, 0, 0);
    scene.doAA = false;
    scene.msaaSamples = 4;
    scene.msaaOptimize = false;
    scene.softShadowSamples = 1;
    scene.traceDepth = 1;
    scene.camera.position.setY(0);
    scene.camera.position.setZ(-20);
    scene.camera.planeDistance = 20;

    scene.photonMapping = true;
    scene.photonCount = 100000; 
    scene.photonBounces = 3;
    scene.photonGatherAmount = 10;
    scene.photonGatherSamples = 256;

    
    Sphere *redSphere = new Sphere(5, Vector(-1, 5, 13));
    Sphere *greenSphere = new Sphere(3, Vector(-7, 7, 8));
    Sphere *blueSphere = new Sphere(3, Vector(10, 7, 15));

    scene.camera.direction = Vector(-1, 5, 13) - scene.camera.position;
    scene.camera.direction.normalize();


//    Box *redSphere = new Box(Vector(-6, -2, 9), Vector(6, 6, 6));
//    Box *greenSphere = new Box(Vector(1, -2, 8), Vector(6, 6, 6));
//    Box *blueSphere = new Box(Vector(-6.5, 4, 7.5), Vector(6, 6, 6));
//    Box *yellowSphere = new Box(Vector(0.5, 4, 7), Vector(6, 6, 6));

    Plane *bottomPlane = new Plane(Vector(0, 10, 0), Vector(0, -1, 0));
    Plane *upPlane = new Plane(Vector(0, 0, 16), Vector(0, 0, -1));
    Plane *leftPlane = new Plane(Vector(-11, 0, 0), Vector(1, 0, 0));
    Plane *rightPlane = new Plane(Vector(14, 0, 0), Vector(-1, 0, 0));

    redSphere->material.color = Vector(1.0, 0.2, 0.2);
    redSphere->material.reflectivity = 0.4;
    redSphere->material.isReflective = false;
    redSphere->material.specular = 0.3;

    greenSphere->material.color = Vector(0.2, 1.0, 0.2);
    greenSphere->material.reflectivity = 0.2;
    greenSphere->material.isReflective = false;
    greenSphere->material.specular = 0.3;

    blueSphere->material.color = Vector(0.2, 0.2, 1.0);
    blueSphere->material.reflectivity = 0.4;
    blueSphere->material.isReflective = false;
    blueSphere->material.specular = 0.3;

/*    yellowSphere->material.color = Vector(1.0, 1.0, 0.2);
    yellowSphere->material.reflectivity = 0.3;
    yellowSphere->material.isReflective = false;
    yellowSphere->material.specular = 0.3;
*/
    bottomPlane->material.color = Vector(1, 1, 1);
    bottomPlane->material.isReflective = true;
    bottomPlane->material.reflectivity = 1;

    upPlane->material.color = Vector(1, 1, 1);
    upPlane->material.isReflective = false;
    upPlane->material.specular = 0;
    upPlane->material.reflectivity = 0.2;

    leftPlane->material.color = Vector(1, 0, 0);
    rightPlane->material.color = Vector(0, 1, 0);

    scene.addRenderable(redSphere);
    scene.addRenderable(greenSphere);
    scene.addRenderable(blueSphere);
//    scene.addRenderable(yellowSphere);
    scene.addRenderable(bottomPlane);
    scene.addRenderable(upPlane);
    scene.addRenderable(leftPlane);
    scene.addRenderable(rightPlane);
    scene.ambientCoefficient = 0.05;

    AreaLight* topLight = new AreaLight();
    topLight->position.set(0, -9, 0);
    topLight->dir1.set(1, 0, 0);
    topLight->dir2.set(0, 0, 1);
    topLight->size1 = 3;
    topLight->size2 = 3;
    topLight->brightness = 0.3;
    scene.addLight(topLight);

    scene.render("test.bmp", NULL);

    return 0;
}
