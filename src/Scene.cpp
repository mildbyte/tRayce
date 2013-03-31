#include "Scene.h"

void Scene::addRenderable(Renderable* renderable) {
    renderables_.addRenderable(renderable);
}

void Scene::addLight(Light* light) {
    lights_.push_back(light);
}

Scene::Scene(int width, int height) {
    rendered_ = new Bitmap(width, height);
    width_ = width;
    height_ = height;

    //Set up the default camera position
    camera.position.set(0, 0, -10);
    camera.direction.set(0, 0, 1);
    camera.width = 24;
    camera.height = 18;
    camera.planeDistance = 10;

    traceDepth = 5;
    backgroundColor.set(0, 0, 0);

    doAA = false;
    msaaSamples = 1;
    msaaOptimize = true;

    softShadowSamples = 1;

    //Set up the AA optimization structures
    prevHit_ = NULL;

    prevRow_ = new Renderable*[width];
    memset(prevRow_, NULL, width * sizeof (Renderable*));

    currRow_ = new Renderable*[width];
    memset(currRow_, NULL, width * sizeof (Renderable*));

    photonMapping = false;
    photonCount = 1000;
    photonBounces = 2;

    photonExposure = 50;
    photonGatherRadius = 0.7;
}

double Scene::calculateShadingCoefficient(Light* light, Vector point, Vector toLight, double lightDist) {
    //Calculates the shading coefficient (how much a point is obscured)
    //Different calculations for different light types
    switch(light->type) {
        case LT_POINT: {
            //Cast a ray from the point to the light source
            Ray pointToLight;
            pointToLight.direction = toLight;
            pointToLight.origin = point + toLight * 0.0001;


            //If the ray intersects something closer than the light, the point
            //is fully shaded, otherwise, it's fully illuminated.
            if (renderables_.intersectsCloser(pointToLight, lightDist)) return 0;
            return 1;

        } break;
        case LT_AREA: {
            if (softShadowSamples == 1) {
                //Treat as a point light if only one sample
                //The following has been copied from the LT_POINT case.

                //Cast a ray from the point to the light source
                Ray pointToLight;
                pointToLight.direction = toLight;
                pointToLight.origin = point + toLight * 0.0001;


                //If the ray intersects something closer than the light, the point
                //is fully shaded, otherwise, it's fully illuminated.
                if (renderables_.intersectsCloser(pointToLight, lightDist)) return 0;
                return 1;
            }

            double contribution = 1.0 / (softShadowSamples*softShadowSamples);
            //Each ray contributes 1/ssSamples^2 of the shading coefficient

            double totalShade = 0;
            AreaLight *currLight = ((AreaLight*)(light));

            //dx, dy: direction of one square on the grid
            Vector dx = currLight->dir1 * (currLight->size1 / softShadowSamples);
            Vector dy = currLight->dir2 * (currLight->size2 / softShadowSamples);

            //Position of the left top corner of the light
            Vector basePos = currLight->position - currLight->dir1 *
                currLight->size1 * 0.5 - currLight->dir2 * currLight->size2 * 0.5;

            Vector x (0, 0, 0);
            for (int gridX = 0; gridX < softShadowSamples; gridX++) {
                Vector y (0, 0, 0);
                for (int gridY = 0; gridY < softShadowSamples; gridY++) {
                    //Current global position on the grid
                    Vector gridPos = basePos + x + y;
                    //Pick a random position in the current square on the grid
                    gridPos += dx * (rand() / RAND_MAX) + dy * (rand() / RAND_MAX);
                    //Calculate the position relative to the point in question
                    gridPos -= point;
                    double distance = gridPos.normalize();
                    Ray pointToLight;

                    pointToLight.direction = gridPos;
                    pointToLight.origin = point + gridPos * 0.0001;

                    if (!renderables_.intersectsCloser(pointToLight, distance)){
                        //If the ray reaches the light safely, increase the
                        //shading value
                        totalShade += contribution;
                    }
                    //Next vertical position
                    y += dy;
                }
                //Next horizontal position
                x += dx;
            }

            return totalShade;

        } break;
    }

    return 0;
}

Vector Scene::calculatePhongColor(Intersection inter, Ray ray) {
    //Implements Phong shading
    //Adds point color at the intersection using Lambert law
    //(the light intensity depends on the angle between the light ray
    //and the normal to a point) and point color due to specular highlight

    Vector totalColor(0, 0, 0); //Total calculated color
    double shade = 0; //Shading coefficient for the current light

    //Iterate through every light in the scene
    for (std::list<Light*>::iterator it = lights_.begin();
        it != lights_.end(); it++) {

        Vector toLight = ((Light*)(*it))->position - inter.coords;
        double distance = toLight.normalize();

        //See if the point is affected by light
        shade = calculateShadingCoefficient(*it, inter.coords, toLight, distance)
            * inter.normal.dot(toLight);

        //If the hit happened from the inside of the object, adjust the dot
        //product (so that we e.g. can light a sphere from the inside
        if (inter.fromTheInside) shade = -shade;
        //Negative shading values (if the material is facing away
        //from the light) are clamped to zero
        if (shade < 0) shade = 0;

        //The reflected light vector
        Vector reflLight = toLight - inter.normal
            * 2.0f * toLight.dot(inter.normal);

        //The more the reflected vector coincides with the view->object ray,
        //the stronger is the specular coefficient.
        double specCoef = ray.direction.dot(reflLight);
        if (specCoef > 0) {
            specCoef = powf(specCoef,
                            inter.object->material.specularExp) * shade;

            //Apply the specular highlight to the color
            totalColor += inter.object->material.color
                * specCoef * inter.object->material.specular;
        }

        //Apply the diffuse lighting
        totalColor += inter.object->material.color *
        ((shade * ((PointLight*)(*it))->brightness * inter.object->material.diffuse) *
         (1 - ambientCoefficient) + ambientCoefficient);
    }

    return totalColor;
}

Vector Scene::calculateRefraction(Intersection inter, Ray ray, int level) {
    //Refracts the ray and sends it further

    //Refraction index
    double n = 1 / inter.object->material.refrIndex;
    double c1 = inter.normal.dot(ray.direction);

    //If the ray hit the primitive from the inside, reverse the refraction index.
    if (inter.fromTheInside) n = 1/n;

    //If negative, the ray hit the primitive from the outside
    //However, we still need to make the cosine positive
    if (c1 < 0) c1 = -c1;

    double c2 = 1 - n * n * (1 - c1 * c1);

    //Bail out if the cosine of the angle squared is > 1 (total internal
    //reflection) or negative (cannot calculate the square root)

    if ((c2 > 1) || (c2 < 0)) return Vector(0, 0, 0);

    c2 = sqrtf(c2);

    //Calculate the direction of the refracted ray
    Vector refrDir = ray.direction * n + inter.normal * (n * c1 - c2);
    refrDir.normalize();

    //Form the refracted ray and trace it
    ray.origin = inter.coords + refrDir * 0.001;
    ray.direction = refrDir;
    return traceRay(ray, level - 1) * inter.object->material.transparency;
}

Vector Scene::calculateReflection(Intersection inter, Ray ray, int level) {
    //Reflects the ray and traces it further

    //Reflect the incident ray in the surface of the primitive
    Vector reflDir = ray.direction - inter.normal
                * 2.0f * ray.direction.dot(inter.normal);

    reflDir.normalize();

    //The origin of the ray is the point of the original collision
    Ray reflRay;
    reflRay.origin = inter.coords + reflDir * 0.001;
    reflRay.direction = reflDir;

    //Send the reflected ray further (and decrease the tracing level)
    return traceRay(reflRay, level - 1) * inter.object->material.reflectivity;
}

//Generates a random vector in the hemisphere of the normal
Vector sampleLambertianBRDF(Vector normal) {
    Vector result(0, 0, 0);
    do {
        result.setX(rand() * 2 - 1);
        result.setY(rand() * 2 - 1);
        result.setZ(rand() * 2 - 1);
        result.normalize();
    } while (normal.dot(result) < 0);
    return result;
}

Vector Scene::traceRay(const Ray ray, int level) {
    //Traces a single ray through the scene; returns its color.
    //This is where the magic happens.

    //If the maximum tracing depth has been reached, return
    if (level <= 0) return backgroundColor;

    //If the ray hits something, return the background color.
    Intersection inter = renderables_.getFirstIntersection(ray);

    if (!inter.happened) {
        prevHit_ = NULL;
        return backgroundColor;
    }

    //Calculate the resultant color of the pixel
    inter.normal.normalize();
    Vector resultColor(0, 0, 0);

    if (photonMapping) {
        //Gathering the photons replaces classic raytracing
        resultColor = photonMap_->gatherPhotons(inter.coords,
            inter.normal, photonExposure, photonGatherRadius, 10);
        resultColor += calculatePhongColor(inter, ray);
/*
        for (int i = 0; i < 5; i++) {
            //Sampling ray from the surface of the entity
            Ray samplingRay;

            do {
                samplingRay.direction.setX(rand() * 2 - 1);
                samplingRay.direction.setY(rand() * 2 - 1);
                samplingRay.direction.setZ(rand() * 2 - 1);
                samplingRay.direction.normalize();
            } while (inter.normal.dot(samplingRay.direction) < 0);


            Intersection sampleInter = renderables_.getFirstIntersection(samplingRay);
            if (!sampleInter.happened) continue;
            sampleInter.normal.normalize();

            resultColor += photonMap_->gatherPhotons(sampleInter.coords,
                sampleInter.normal, photonExposure, photonGatherRadius, 10);
        }*/
//        resultColor *= 0.05;
    } else {
        //Phong (diffuse+specular) pass
        resultColor += calculatePhongColor(inter, ray);

        //Add the transmitted ray
        if (inter.object->material.isTransparent)
            resultColor += calculateRefraction(inter, ray, level);

        //Add the reflected ray
        if (inter.object->material.isReflective)
            resultColor += calculateReflection(inter, ray, level);
    }

    //Store the distance (used for postprocessing)
    prevDist_ = inter.distance;
    prevHit_ = inter.object;
    return resultColor;
}

//Converts screen to image coordinates and traces them.
Vector Scene::tracePixel(double x, double y) {
    Ray ray;

    //Uses conic projection, casts the rays from the same point.
    ray.origin = camera.position;

    Vector xWorld = xPixel * (x - 0.5 * width_ + 0.5);
    Vector yWorld = yPixel * (y - 0.5 * height_ + 0.5);

    //First move to the centre of the image plane, then in the image plane
    //to the needed point.
    ray.direction = (camera.direction * camera.planeDistance)
                  + xWorld + yWorld;
    ray.direction.normalize();

    return traceRay(ray, traceDepth);
}

inline double min(double a, double b) {return a < b ? a : b;}

Vector combineColors(Vector c1, Vector c2) {
    return Vector(min(c1.getX(), c2.getX()), min(c1.getY(), c2.getY()),
                  min(c1.getZ(), c2.getZ()));
}

void Scene::populatePhotonMap() {
    //Construct the map
    //There will be one photon per bounce
    photonMap_ = new PhotonMap(photonCount * photonBounces);

    //The number of photons emitted per light depends on the light's intensity
    double totalIntensity = 0;
    
    int photonsHit = 0;
    int allHits = 0;

    for (std::list<Light*>::iterator it = lights_.begin();
        it != lights_.end(); it++) totalIntensity += ((Light*)(*it))->brightness;


    for (std::list<Light*>::iterator it = lights_.begin();
        it != lights_.end(); it++) {
        int photonsToCast = (int)((((Light*)(*it))->brightness/totalIntensity) 
                          * (double)photonCount);
        printf("Casting %d photons...\n", photonsToCast);

        for (int i = 0; i < photonsToCast; i++) {
            //Make the light->scene ray
            Ray photonRay;
            photonRay.origin = ((Light*)(*it))->position;
            
            //Generate a random direction for our ray
            photonRay.direction.setX(rand() * 2 - 1);
            photonRay.direction.setY(rand() * 2 - 1);
            photonRay.direction.setZ(rand() * 2 - 1);
            photonRay.direction.normalize();
            
            Vector photonEnergy(1, 1, 1);

            int currBounces = 1;

            Intersection inter = renderables_.getFirstIntersection(photonRay);
            if (inter.happened) photonsHit++;

            while (inter.happened && currBounces < photonBounces) {
                //Record the photon               
                photonEnergy = combineColors(photonEnergy, 
                                             inter.object->material.color);
                photonMap_->addPhoton(inter.coords, photonRay.direction, 
                                      photonEnergy);

                //New point to cast the ray from
                photonRay.origin = inter.coords;
                
                inter.normal.normalize(); //not normalized by the intersection
                                          //to save CPU cycles

                //Reflect the ray
                photonRay.direction = sampleLambertianBRDF(inter.normal);

                //Send the ray onwards
                inter = renderables_.getFirstIntersection(photonRay);
                if (inter.happened) allHits++;
                currBounces++;
            }
        }
    }
    
    photonMap_->makeTree();

    printf("Total primary hits: %d\n", photonsHit);
    printf("Total hits: %d\n", photonsHit + allHits);
}

void Scene::render(char* filename, BitmapPixel (*postProcess)(BitmapPixel)) {
    //Renders the scene to a file

    //Populate the photon map (should really be done once per scene, not every render,
    //is here for testing).
    if (photonMapping) {
        printf("Populating the photon map...\n");
        populatePhotonMap();
    }

    printf("Rendering...\n");

    //imagePlaneX and imagePlaneY are the unit x and y vectors in the image plane.
    Vector imagePlaneX(
        camera.direction.getZ(),
        0,
        -camera.direction.getX());
    
    Vector imagePlaneY(
         -camera.direction.getY()*camera.direction.getX(),
         camera.direction.getZ()*camera.direction.getZ() 
            + camera.direction.getX() * camera.direction.getX(),
         -camera.direction.getZ()*camera.direction.getY());

    imagePlaneX.normalize();
    imagePlaneY.normalize();

    //xPixel and yPixel: one pixel in the image plane
    xPixel = imagePlaneX * (camera.width / (double)width_);
    yPixel = imagePlaneY * (camera.height / (double)height_);

    //The resulting color of the pixel
    Vector resultColor;

    //Renderable that has been hit previously
    Renderable* prevHit = NULL;
    
    //Each MSAA sample's contribution to the final pixel
    double contribution = 1.0 / (msaaSamples*msaaSamples);

    for (int realy = 0; realy < height_; realy++) {
        for (int realx = 0; realx < width_; realx++) {
            //Trace a preemptive ray through the pixel on the image plane
            resultColor = tracePixel(realx, realy);

            currRow_[realx] = prevHit_;

            //Runs if either no optimizations enabled or a new object is hit
            //either in the pixel on the left or above the current pixel
            if (doAA && (!msaaOptimize || (prevHit_ != prevHit)
                         || (prevHit_ != prevRow_[realx]))){
                //Antialiased image
                //Divides every pixel into an msaaSamples x msaaSamples grid
                //Traces a ray through the centre of each square
                prevHit = prevHit_;

                resultColor.set(0, 0, 0);

                //Traverse the grid
                for (int i = 1; i <= msaaSamples; i++) {
                    for (int j = 1; j <= msaaSamples; j++) {
                        //Trace the ray (contributes only a part of the
                        //resultant color)
                        resultColor += tracePixel(
                            realx + 1.0/(msaaSamples + 1.0) * (i + 0.5) - 0.5,
                            realy + 1.0/(msaaSamples + 1.0) * (j + 0.5) - 0.5)
                            * contribution;
                    }
                }
                //Mark the antialiased pixels (for debugging purposes)
                //resultColor.set(1, 1, 1);
            }

            //Set the result color (either AA'd or w/o AA)
            rendered_->setPixel(realx, realy, resultColor, prevDist_);

            //Copy the collision data in the current row
            memcpy(prevRow_, currRow_, sizeof(Renderable*) * width_);
        } 
    } 
    if (postProcess) rendered_->foreach(postProcess);

    //Save the resulting bitmap to file
    rendered_->saveToFile(filename);
}
