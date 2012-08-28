#include <iostream>

#include "Scene.h"
#include "Sphere.h"
#include "Box.h"
#include "Plane.h"
#include <ctime>

using namespace std;

/*bool collides(Sphere sphere1, Sphere sphere2) {
    Vector distvect = sphere1.position_ - sphere2.position;
    double distance = distvect.normalize();
    return distvect >= (sphere1.radius_ + sphere2.radius_);
}*/

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

BitmapPixel noPostProcess (BitmapPixel pix) { return pix; }

int main()
{
    srand((unsigned)time(0));
    Scene scene(1280, 800);
    scene.camera.height = 15;
    scene.backgroundColor.set(0, 0, 0);
    scene.doAA = false;
    scene.msaaSamples = 4;
    scene.msaaOptimize = true;
    //scene.camera.height = 15;
    scene.softShadowSamples = 1;
    scene.traceDepth = 4;
    scene.camera.position.setY(2);
    scene.camera.position.setZ(-20);
    scene.camera.planeDistance = 20;

    
    Sphere *redSphere = new Sphere(5, Vector(-1, 5, 13));
    Sphere *greenSphere = new Sphere(3, Vector(-7, 7, 8));
    Sphere *blueSphere = new Sphere(3, Vector(10, 7, 15));


//    Box *redSphere = new Box(Vector(-6, -2, 9), Vector(6, 6, 6));
//    Box *greenSphere = new Box(Vector(1, -2, 8), Vector(6, 6, 6));
//    Box *blueSphere = new Box(Vector(-6.5, 4, 7.5), Vector(6, 6, 6));
//    Box *yellowSphere = new Box(Vector(0.5, 4, 7), Vector(6, 6, 6));

    Plane *bottomPlane = new Plane(Vector(0, 10, 0), Vector(0, -1, 0));
    Plane *upPlane = new Plane(Vector(0, 0, 16), Vector(0, 0, -1));

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
    bottomPlane->material.color = Vector(0.3, 0.3, 0.3);
    bottomPlane->material.isReflective = true;
    bottomPlane->material.reflectivity = 1;

    upPlane->material.color = Vector(0.3, 0.3, 0.3);
    upPlane->material.isReflective = false;
    upPlane->material.specular = 0;
    upPlane->material.reflectivity = 0.2;

    scene.addRenderable(redSphere);
    scene.addRenderable(greenSphere);
    scene.addRenderable(blueSphere);
//    scene.addRenderable(yellowSphere);
    scene.addRenderable(bottomPlane);
    scene.addRenderable(upPlane);
    scene.ambientCoefficient = 0.05;

    AreaLight* topLight = new AreaLight();
    topLight->position.set(10, -9, 0);
    topLight->dir1.set(1, 0, 0);
    topLight->dir2.set(0, 0, 1);
    topLight->size1 = 3;
    topLight->size2 = 3;
    topLight->brightness = 0.3;
    scene.addLight(topLight);

    AreaLight* topLight2 = new AreaLight();
    topLight2->position.set(-1.5, -10, 10);
    topLight2->dir1.set(1, 0, 0);
    topLight2->dir2.set(0, 0, -1);
    topLight2->size1 = 3;
    topLight2->size2 = 3;
    topLight2->brightness = 0.3;
    scene.addLight(topLight2);

    AreaLight* topLight3 = new AreaLight();
    topLight3->position.set(-3, 0, 0);
    topLight3->dir1.set(0, -1, 0);
    topLight3->dir2.set(1, 0, 0);
    topLight3->size1 = 6;
    topLight3->size2 = 6;
    topLight3->brightness = 0.3;
    scene.addLight(topLight3);

    scene.render("test.bmp", &noPostProcess);

    return 0;
}
