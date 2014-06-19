#include <iostream>

#include "Scene.h"
#include "Sphere.h"
#include "AABB.h"
#include "Plane.h"
#include "Triangle.h"
#include <ctime>
#include "Random.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

//Adds a model inside an OBJ file into the scene with a given
//material. Does not reuse vertices, so each triangle takes
//3xvertexsize + overhead space. Only supports triangles.
void importObj(Scene* scene, char* filename, Material m,
    Vector shift, double scale) {
    vector<Vector> vectors;
    vector<Vector> normalVectors;
    
    ifstream stream;
    stream.open(filename);
    
    int faces = 0;
    
    while (!stream.eof()) {
        string type;
        stream >> type;
        
        if (type == "v") {
            double x, y, z;
            stream >> x >> y >> z;
            vectors.push_back(Vector(x, y, z) * scale + shift);
        } else if (type == "vn") {
            double x, y, z;
            stream >> x >> y >> z;
            normalVectors.push_back(Vector(x, y, z));
        } else if (type == "f") {
            faces++;
            //Face format: f vertex vertex vertex
            //vertex: v/t/n
            //where v is the vertex id, t is the texture UV coordinate id
            //n is the normal direction at the vertex
            //Any of the second or the third can be skipped:
            //v, v/t, v//n are all valid
            
            int vIds[3];
            int normalvIds[3] = {-1, -1, -1};
            
            for (int i = 0; i < 3; i++) {
                string vertices;
                stream >> vertices;
                
                //Split on '/'
                istringstream ss(vertices);
                string tok;
                getline(ss, tok, '/'); //First number: vertex id
                vIds[i] = atoi(tok.c_str());
                
                if (!ss.eof()) getline(ss, tok, '/'); //Second number: texture UV coordinates (discard)
                
                if (!ss.eof()) {
                    getline(ss, tok, '/'); //Last number: normal vertex id
                    normalvIds[i] = atoi(tok.c_str());
                }
            }
            
            Triangle *t;
            
            //Vertex ids in OBJ are 1-based (so 0 never appears in a face description)
            if (normalvIds[0] != -1 && normalvIds[1] != -1 && normalvIds[2] != -1) {
                t = new Triangle(vectors[vIds[0]-1], vectors[vIds[1]-1], vectors[vIds[2]-1],
                                 normalVectors[normalvIds[0]-1], normalVectors[normalvIds[1]-1], normalVectors[normalvIds[2]-1]);
            } else {
                t = new Triangle(vectors[vIds[0]-1], vectors[vIds[1]-1], vectors[vIds[2]-1]);
            }
            
            t->material = m;
            scene->addTriangle(t);
            //scene->addRenderable(t);
        } else {
            stream.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }
    }
    
    cout << "Loaded " << vectors.size() << " vertices, " << normalVectors.size() <<
        " normal vectors and " << faces << " faces." << endl;
    
    stream.close();
}

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
    scene.camera.width = 8;
    scene.camera.height = 4.5;

    scene.backgroundColor.set(0, 0, 0);
    //scene.doAA = true;
    //scene.msaaSamples = 4;
    scene.msaaOptimize = false;
    scene.softShadowSamples = 1;
    scene.traceDepth = 5;
	scene.camera.position = Vector(0, 3, -10);
	scene.camera.direction = Vector(0, -6.2, 10);
	scene.camera.direction.normalize();
    scene.camera.planeDistance = 15;
    scene.camera.lensRadius = 0;
    scene.camera.focalDistance = 25;
	
    scene.renderingMode = PATHTRACING;
    scene.pathTracingSamplesPerPixel = 16; //spp squared is actually cast
    scene.pathTracingMaxDepth = 10; // Too few samples and rays that go through
    // a sphere, bounce off a wall, through the sphere again and to the light
    // will terminate too early
    
    scene.doFinalGather = true;
	scene.visualizePhotons = true;
	scene.photonCount = 65536;
    scene.photonBounces = 5;
    scene.photonGatherAmount = 32;
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

    Sphere *oneLight = new Sphere(10, Vector(0, 10, -10));

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
    //m.isTransparent = true;
    //m.transparency = 1.0;
	m.reflectivity = 0.5;
    //m.diffuse = 1.0;
    //m.refrIndex = 1.42;
    m.color = Vector(0.929, 0.5, 0.01);
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
    importObj(&scene, "torus.obj", m, Vector(0, -10, 10), 4.0);
    
    scene.render("test.bmp", NULL, 8);

    return 0;
}
