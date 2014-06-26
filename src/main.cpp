#include <iostream>

#include "Scene.h"
#include "Sphere.h"
#include "AABB.h"
#include "Plane.h"
#include "Triangle.h"
#include <ctime>
#include "Random.h"

using namespace std;

//Generates a random scene with a bunch of spheres
void randomScene(Scene* scene) {
    for (int i = 0; i < 30; i++) {
        double r = drand() * 2.0 + 0.3;
        Sphere *testSphere = new Sphere(
            r,
            Vector(drand() * (25.0-2*r) - (11.0-r), drand() * (20.0-2*r) - (10.0-r), drand() * (18.0-2*r) + r)
        );

        testSphere->material.color = Vector(drand(), drand(), drand());
        
        if (drand() > 0.8) {
            testSphere->material.reflectivity = drand();
        }
        
        if (drand() > 0) {
            testSphere->material.transparency = drand() * (1.0 - testSphere->material.reflectivity);
            testSphere->material.refrIndex = 1.42;
        }
        
        
        /*if (!testSphere->material.isTransparent && !testSphere->material.isReflective) {
            if (drand() > 0.5) {
                testSphere->material.emittance = Vector(drand(), drand(), drand()) * 10.0;
                
            }
        }*/
        
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
    Scene scene(640, 360);
    //scene.camera.width = 32;
    scene.camera.width = 16;
    scene.camera.height = 9;

    scene.backgroundColor.set(0, 0, 0);
    scene.softShadowSamples = 1;
    scene.traceDepth = 5;
	scene.camera.position = Vector(0, -6, -15);
	scene.camera.direction = Vector(0, 0, 1);
	scene.camera.direction.normalize();
    scene.camera.planeDistance = 15;
    scene.camera.lensRadius = 0;
    scene.camera.focalDistance = 25;
	
    scene.renderingMode = PHOTONMAPPING;
    scene.pathTracingSamplesPerPixel = 16; //spp squared is actually cast
    scene.pathTracingMaxDepth = 5; // Too few samples and rays that go through
    // a sphere, bounce off a wall, through the sphere again and to the light
    // will terminate too early
    
    scene.doFinalGather = true;
	scene.visualizePhotons = false;
	scene.photonCount = 100000;
    scene.photonBounces = 7;
    scene.photonGatherAmount = 128;
    scene.photonGatherSamples = 4;
    scene.irradiancePhotonFrequency = 8;
    scene.photonGatherDotThreshold = 0.9;
    scene.doIrradianceCaching = true;
    scene.samplingMode = STRATIFIED;
    
    // Problem with this scene: when we do a final gather in the corner
    // shaded by the big sphere, most rays hit the sphere. We then look up
    // an irradiance photon with a similar normal, but the closest one we find
    // is one on the emitting sphere. Hence we do light the corner that's supposed
    // to be shaded. Temporary solution: more photons so that at least some get
    // to the corner. Better solution?

    Sphere *oneLight = new Sphere(7, Vector(0, 15, 10));

    Plane *bottomPlane = new Plane(Vector(0, -10, 0), Vector(0, 1, 0));
    Plane *upPlane = new Plane(Vector(0, 0, 18), Vector(0, 0, -1));
    Plane *leftPlane = new Plane(Vector(-10.0, 0, 0), Vector(1, 0, 0));
    Plane *rightPlane = new Plane(Vector(10.0, 0, 0), Vector(-1, 0, 0));
    Plane *topPlane = new Plane(Vector(0, 10, 0), Vector(0, -1, 0));
    Plane *backPlane = new Plane(Vector(0, 0, -10), Vector(0, 0, 1));
    
    Light *testLight = new PointLight();
    testLight->position = Vector(0, 12, 5);
    testLight->brightness = 5.0;
  
    oneLight->material.color = Vector(.95, .95, .95);
    oneLight->material.emittance = Vector(30, 30, 30);

    bottomPlane->material.color = Vector(.95, .95, .95);
    upPlane->material.color = Vector(.95, .95, .95);
    leftPlane->material.color = Vector(.95, .95, .95);
    rightPlane->material.color = Vector(.95, .95, .95);
    topPlane->material.color = Vector(.95, .95, .95);
    //topPlane->material.emittance = Vector(5, 5, 5);
    backPlane->material.color = Vector(.95, .95, .95);
    
    Material m;
    m.transparency = 1.0;
	//m.reflectivity = 0.5;
    //m.diffuse = 1.0;
    m.refrIndex = 1.42;
    m.color = Vector(0.01, 0.5, 0.929);
	//m.color = Vector(1, 1, 1);
/*    
    Triangle *t1 = new Triangle(Vector(4, 10, 13), Vector(8, -3, 9), Vector(12, 10, 13));
    Triangle *t2 = new Triangle(Vector(12, 10, 13), Vector(8, -3, 9), Vector(8, 10, 5));
    Triangle *t3 = new Triangle(Vector(8, 10, 5), Vector(4, 10, 13), Vector(8, -3, 9));
    
    t1->material = m;
    t2->material = m;
    t3->material = m;
    */
    scene.addRenderable(oneLight);
    
    scene.addLight(testLight);
    
    //scene.addRenderable(t1);
    //scene.addRenderable(t2);
    //scene.addRenderable(t3);

    scene.addRenderable(bottomPlane);
    scene.addRenderable(upPlane);
    //scene.addRenderable(leftPlane);
    //scene.addRenderable(rightPlane);
    //scene.addRenderable(topPlane);
    //scene.addRenderable(backPlane);
    
    //seed_drand(39332);
    //randomScene(&scene);
    
    printf("Loading the object file...\n");
	//scene.importObj("sphere.obj", m, Vector(0, -10, 5), 4.0);

	Sphere *s = new Sphere(4, Vector(0, -6, 5));
	s->material = m;
	scene.addRenderable(s);
    
    scene.render("test.bmp", NULL, 8);

    return 0;
}
