#include "Scene.h"

//Multiplies all elements of the two colors together
Vector combineColors(Vector c1, Vector c2) {
    return Vector(c1.getX() * c2.getX(), c1.getY() * c2.getY(),
                  c1.getZ() * c2.getZ());
}

//Returns the indexth number of base 'base' Halton sequence; adapted from
//http://en.wikipedia.org/wiki/Halton_sequence
double halton(int index, int base) {
    double result = 0.0;
    double f = 1.0 / (double)base;
    int i = index;
    while (i) {
        result += f * (i % base);
        i /= base;
        f /= base;
    }
    return result;
}

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

    photonMap_ = NULL;

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

    //Some default settings
    doFinalGather = false;
    visualizePhotons = false;
    photonCount = 1000;
    photonBounces = 2;
    
    photonGatherAmount = 10;
    photonGatherSamples = 10;

    photonGatherDotThreshold = 0.9;
    irradiancePhotonFrequency = 4;
    
    pathTracingMaxDepth = 5;
    pathTracingSamplesPerPixel = 100;

    samplingMode = STRATIFIED;
    renderingMode = RAYTRACING;
}

double Scene::calculateShadingCoefficient(Light* light, Vector point, Vector toLight, double lightDist) {
    //Calculates the shading coefficient (how much a point is obscured)
    //Different calculations for different light types
    switch(light->type) {
        case LT_POINT: {
            //Cast a ray from the point to the light source
            Ray pointToLight;
            pointToLight.direction = toLight;
            pointToLight.origin = point;

             pointToLight = epsilonShift(pointToLight);

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
                pointToLight.origin = point;

                pointToLight = epsilonShift(pointToLight);

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
                    pointToLight.origin = point;

                    pointToLight = epsilonShift(pointToLight);

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
        if (ray.direction.dot(inter.normal) > 0) shade = -shade;
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
             
    //Add the object's emittance
    totalColor += inter.object->material.emittance;

    return totalColor;
}

Ray Scene::epsilonShift(Ray ray) {
    Ray result = ray;
    result.origin += result.direction * 0.0001;

    return result;
}

bool Scene::refractRay(Intersection inter, Ray& ray) {
    //Refraction index
    double n = 1 / inter.object->material.refrIndex;
    double c1 = inter.normal.dot(ray.direction);

    //Make c1 positive. If it was positive, the ray hit the sphere from the inside
    //and so the normal must be inverted.
    if (c1 < 0) {
        c1 = -c1;
    } else {
        inter.normal = -inter.normal;
        n = 1/n;
    }

    double c2 = 1 - n * n * (1 - c1 * c1);

    //Bail out if the cosine of the angle squared is > 1 (total internal
    //reflection) or negative (cannot calculate the square root)

    if ((c2 > 1) || (c2 < 0)) return false;

    c2 = sqrtf(c2);

    //Calculate the direction of the refracted ray
    Vector refrDir = ray.direction * n + inter.normal * (n * c1 - c2);
    refrDir.normalize();

    ray.origin = inter.coords;
    ray.direction = refrDir;
    
    return true;
}

Vector Scene::calculateRefraction(Intersection inter, Ray ray, int level) {
    //Refracts the ray and sends it further
    //Form the refracted ray and trace it
    if (refractRay(inter, ray)) {
        return traceRay(epsilonShift(ray), level - 1) * inter.object->material.transparency;
    } else {
        return backgroundColor;
    }
}

Ray Scene::reflectRay(Intersection inter, Ray ray) {
    //Reflect the incident ray in the surface of the primitive
    Vector reflDir = ray.direction - inter.normal
                * 2.0f * ray.direction.dot(inter.normal);

    reflDir.normalize();

    //The origin of the ray is the point of the original collision
    Ray result;

    result.origin = inter.coords;
    result.direction = reflDir;

    return result;
}

Vector Scene::calculateReflection(Intersection inter, Ray ray, int level) {
    //Send the reflected ray further (and decrease the tracing level)
    return traceRay(epsilonShift(reflectRay(inter, ray)), level - 1) 
           * inter.object->material.reflectivity;
}

/*
Vector sampleHemisphere(Vector normal, double a, double b) {
    Vector result;
    do {
        result = Vector(2*drand()-1, 2*drand()-1, 2*drand()-1);
        if (result.normalize() > 1) continue;
        if (result.dot(normal) < 0) continue;
        return result;
    } while (true);
}*/

Vector sampleHemisphere2(Vector normal) {
    while (true) {
        Vector direction(2*drand()-1, 2*drand()-1, 2*drand()-1);
        
        double dot = normal.dot(direction);
        if (dot >= 0.0 && dot <= 1.0) {
            direction.normalize();
            return direction;
        }
    }
}

Vector sampleHemisphere(Vector normal, double Xi1, double Xi2) {
    //Cosine-weighted hemisphere sampling.
    //Adapted from http://pathtracing.wordpress.com/2011/03/03/cosine-weighted-hemisphere/
    double theta = acos(sqrt(1.0-Xi1));
    double phi = 2.0 * 3.1415926535897932384626433832795 * Xi2;

    double xs = sin(theta) * cos(phi);
    double ys = cos(theta);
    double zs = sin(theta) * sin(phi);

    Vector y = normal;
    Vector h = y;

    if (abs(h.getX())<=abs(h.getY()) && abs(h.getX())<=abs(h.getZ()))
    h.setX(1.0);
    else if (abs(h.getY())<=abs(h.getX()) && abs(h.getY())<=abs(h.getZ()))
    h.setY(1.0);
    else
    h.setZ(1.0);

    Vector x = (h.cross(y));
    Vector z = (x.cross(y));
    x.normalize();
    z.normalize();
    
    Vector direction = x * xs+ y * ys + z * zs;
    direction.normalize();
    return direction;
}

Vector Scene::sampleMapAt(Vector coords, Vector normal, double x, double y) {
    //Sampling ray from the surface of the entity
    Ray samplingRay;

    samplingRay.origin = coords;
    samplingRay.direction = sampleHemisphere(normal, x, y);
    
    double dot = normal.dot(samplingRay.direction);
    
    //Shift the origin slightly to avoid collisions with itself
    samplingRay = epsilonShift(samplingRay);

    Intersection sampleInter = renderables_.getFirstIntersection(samplingRay, 0);
    if (!sampleInter.happened) return backgroundColor;
    sampleInter.normal.normalize();

    //Read the precomputed radiance value
    return photonMap_->acceleratedIrradiance(sampleInter.coords, 
                                             sampleInter.normal, 
                                             photonGatherDotThreshold) * dot;
}

Vector Scene::pathTrace(const Ray ray, int depth) {
    //If the maximum tracing depth has been reached, return
    if (depth <= 0) return backgroundColor;
    
    Intersection inter;
    
    //If we are casting a ray from the eye into the scene, cull collisions behind the image plane
    if (depth == pathTracingMaxDepth) {
        inter = renderables_.getFirstIntersection(ray, camera.planeDistance);
    } else {
        inter = renderables_.getFirstIntersection(ray, 0);
    }

    if (!inter.happened) {
        return backgroundColor;
    }

    inter.normal.normalize();
    
    if (inter.object->material.isTransparent) {

        Vector refracted;
        Ray nextRay = ray;
        
        if (!refractRay(inter, nextRay)) {
            nextRay = reflectRay(inter, nextRay);
        }

        return inter.object->material.emittance + pathTrace(epsilonShift(nextRay), depth - 1);
    } else {    
        Ray nextRay;
        nextRay.origin = inter.coords;
        nextRay.direction = sampleHemisphere2(inter.normal);
        
        nextRay = epsilonShift(nextRay);
        
        Vector brdf = inter.object->material.color * (2.0 * nextRay.direction.dot(inter.normal));
        Vector reflected = pathTrace(nextRay, depth - 1);
        
        return inter.object->material.emittance + combineColors(brdf, reflected);
    }
}

Vector Scene::traceRay(const Ray ray, int level) {
    //Traces a single ray through the scene; returns its color.
    //This is where the magic happens.
    
    //If the maximum tracing depth has been reached, return
    if (level <= 0) return backgroundColor;

    //If the ray hits something, return the background color.
    Intersection inter = renderables_.getFirstIntersection(ray, camera.planeDistance);

    if (!inter.happened) {
        prevHit_ = NULL;
        return backgroundColor;
    }

    inter.normal.normalize();
    
    Vector resultColor(0, 0, 0);
    
    switch (renderingMode) {
        case RAYTRACING:
             //Phong (diffuse+specular) pass
            resultColor += calculatePhongColor(inter, ray);

            //Add the transmitted ray
            if (inter.object->material.isTransparent)
                resultColor += calculateRefraction(inter, ray, level);

            //Add the reflected ray
            if (inter.object->material.isReflective)
                resultColor += calculateReflection(inter, ray, level);
        break;
        case PHOTONMAPPING:
            //Gathering the photons replaces classic raytracing
            if (visualizePhotons) {
                resultColor = photonMap_->visualizePhoton(inter.coords, 0.01);
            } else if (!doFinalGather) {
            resultColor = photonMap_->acceleratedIrradiance(inter.coords,
                            inter.normal, photonGatherDotThreshold);
            } else {
                switch (samplingMode) {
                    case STRATIFIED: {
                        double squareSide = 1.0 / photonGatherSamples;
                        for (int y = 0; y < photonGatherSamples; y++) {
                            double ybase = squareSide * y;
                            for (int x = 0; x < photonGatherSamples; x++) {
                                //Jittered position in the grid
                                double xbase = squareSide * x;
                        
                                double xPos = xbase + sqrt(drand()) * squareSide;
                                double yPos = ybase + drand() * squareSide;

                                resultColor += sampleMapAt(inter.coords, inter.normal, xPos, yPos);
                            }
                        }
                        resultColor /= (double)(photonGatherSamples * photonGatherSamples);
                        break;
                    }
                    case HALTON: {
                        for (int sampleIndex = 0; sampleIndex < photonGatherSamples; sampleIndex++) {
                            resultColor += sampleMapAt(inter.coords, inter.normal,
                                                       haltonXCoords_[sampleIndex],
                                                       haltonYCoords_[sampleIndex]);
                        }

                        resultColor /= (double)photonGatherSamples;
                        break;
                    }
                }
                //Combine with the color of the object
                resultColor = combineColors(inter.object->material.color, resultColor);
                
                //Add the raytraced color (direct illumination)
                             //Phong (diffuse+specular) pass
                resultColor += calculatePhongColor(inter, ray);

                //Add the transmitted ray
                if (inter.object->material.isTransparent)
                    resultColor += calculateRefraction(inter, ray, level);

                //Add the reflected ray
                if (inter.object->material.isReflective)
                    resultColor += calculateReflection(inter, ray, level);
            }
        break;       
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

    Vector xWorld = xPixel_ * (x - 0.5 * (double)width_ + 0.5);
    Vector yWorld = yPixel_ * (y - 0.5 * (double)height_ + 0.5);

    //First move to the centre of the image plane, then in the image plane
    //to the needed point.
    ray.direction = (camera.direction * camera.planeDistance)
                  + xWorld + yWorld;
    ray.direction.normalize();
    
    if (renderingMode == PATHTRACING) {
        Vector resultColor(0, 0, 0);
        for (int i = 0; i < pathTracingSamplesPerPixel; i++) {
            resultColor += pathTrace(ray, pathTracingMaxDepth);
        }
        resultColor /= (double)pathTracingSamplesPerPixel;
        return resultColor;
    } else {
        return traceRay(ray, traceDepth);
    }
}

void Scene::populatePhotonMap() {
    //Construct the map
    //There will be one photon per bounce at most
    photonMap_ = new PhotonMap(photonCount * photonBounces);

    //The number of photons emitted per light depends on the light's intensity
    double totalIntensity = 0;
    
    int allHits = 0;

    //for (std::list<Light*>::iterator it = lights_.begin();
    //    it != lights_.end(); it++) totalIntensity += ((Light*)(*it))->brightness;
    
    for (std::list<Renderable*>::iterator it = renderables_.begin();
        it != renderables_.end(); it++) totalIntensity += ((Renderable*)(*it))->material.emittance.modulus();


//    for (std::list<Light*>::iterator it = lights_.begin();
//        it != lights_.end(); it++) {
//        int photonsToCast = (int)((((Light*)(*it))->brightness/totalIntensity) 
//                          * (double)photonCount);

    for (std::list<Renderable*>::iterator it = renderables_.begin();
        it != renderables_.end(); it++) {
        
        int photonsToCast = (int)((((Renderable*)(*it))->material.emittance.modulus()/totalIntensity)
                          * (double)photonCount);
        
        printf("Casting %d photons...\n", photonsToCast);

        for (int i = 0; i < photonsToCast; i++) {
            //Make the light->scene ray
            Ray photonRay;
            //photonRay.origin = ((Light*)(*it))->position;
            
            photonRay.origin = ((Renderable*)(*it))->sampleSurface();
            photonRay.direction = sampleHemisphere(((Renderable*)(*it))->getNormalAt(photonRay.origin), drand(), drand());
            photonRay = epsilonShift(photonRay);

            //double brightness = ((Light*)(*it))->brightness;
            
            Vector photonEnergy(((Renderable*)(*it))->material.emittance);

            int currBounces = 0;

            Intersection inter = renderables_.getFirstIntersection(photonRay, 0);

            while (inter.happened && currBounces < photonBounces) {
                currBounces++;
                
                double randVar = drand();

                Material objMat = inter.object->material;

                double avgDiffuse = (objMat.color.getX() + objMat.color.getY()+
                                    objMat.color.getZ()) / 3.0;

                if (randVar > avgDiffuse) break; //The photon has been absorbed
                allHits++;
                
                //Record the photon							 
                inter.normal.normalize(); //not normalized by the intersection
                                          //to save CPU cycles
                
                photonEnergy *= (1.0 / avgDiffuse);

                photonEnergy = combineColors(photonEnergy, 
                                             inter.object->material.color);
				
                photonMap_->addPhoton(inter.coords, photonRay.direction, 
                                      photonEnergy, inter.normal);


                //Diffuse the ray
                photonRay.direction = sampleHemisphere(inter.normal, drand(), drand());
                
                //New point to cast the ray from (+ avoid collision with itself)
                photonRay.origin = inter.coords;
                photonRay = epsilonShift(photonRay);

                //Send the ray onwards
                inter = renderables_.getFirstIntersection(photonRay, 0);
            }
        }
    }
    
    printf("Total hits: %d\n", allHits);
    
    //The energy of the light is spread evenly amongst all photons
    photonMap_->scalePhotonPower(1.0/allHits);

    photonMap_->makeTree();

    printf("Precalculating irradiance...\n");
    photonMap_->precalculateIrradiance(irradiancePhotonFrequency, photonGatherAmount);
}

bool Scene::loadMap(char* path) {
    //Does the map file exist?
    ifstream mapFile(path);
    if (!mapFile) return false; 
    mapFile.close();

    photonMap_ = PhotonMap::makeFromFile(path);
    if (photonMap_ == NULL) return false;
    return true;
}

void Scene::saveMap(char* path) {
    photonMap_->saveToFile(path);
}

void Scene::threadDoWork(int threadId) {
	//A rendering task executed in a separate thread. Renders the pixels assigned to it
	//and stores them in the bitmap.
	
    //The resulting color of the pixel
    Vector resultColor;

    //Renderable that has been hit previously
    //Renderable* prevHit = NULL;
    
    //Each MSAA sample's contribution to the final pixel
    double contribution = 1.0 / (msaaSamples*msaaSamples);
	
    for (int pixel = threadStartPixels_[threadId]; pixel < threadEndPixels_[threadId]; pixel++) {
		int realx = pixel % width_;
		int realy = pixel / width_;

		//Trace a preemptive ray through the pixel on the image plane
		resultColor = tracePixel(realx, realy);

		currRow_[realx] = prevHit_;
		
		//Uses global state, MSAA optimizations disabled for now
		if (doAA) {
		//Runs if either no optimizations enabled or a new object is hit
		//either in the pixel on the left or above the current pixel
		//if (doAA && (!msaaOptimize || (prevHit_ != prevHit)
		//			 || (prevHit_ != prevRow_[realx]))){
			//Antialiased image
			//Divides every pixel into an msaaSamples x msaaSamples grid
			//Traces a ray through the centre of each square
			//prevHit = prevHit_;

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
		//TODO: prevDist_ is global state, remove
		rendered_->setPixel(realx, realy, resultColor, prevDist_);
		
		//Changes global state: bad for threading. MSAA optimizations disabled for now
		////Copy the collision data in the current row
		//memcpy(prevRow_, currRow_, sizeof(Renderable*) * width_);

		pixelsRendered_++;
		if (pixelsRendered_ % 1000 == 0) {
			printf("Rendered %d pixel(s) out of %d (%f\%)\n", pixelsRendered_, 
				(width_ * height_),
				pixelsRendered_ / (double)(width_ * height_) * 100);
		}
    } 
}

void Scene::render(char* filename, BitmapPixel (*postProcess)(BitmapPixel), int noThreads) {
    //Renders the scene to a file

    //Do we have the photon map yet?
    if (renderingMode == PHOTONMAPPING && photonMap_ == NULL) {
        printf("Populating the photon map...\n");
        populatePhotonMap();
    }

    //Populate the Halton sequence cache if we are using Halton sampling
    if (samplingMode == HALTON) {
        haltonXCoords_ = new double[photonGatherSamples];
        haltonYCoords_ = new double[photonGatherSamples];

        for (int i = 0; i < photonGatherSamples; i++) {
            haltonXCoords_[i] = halton(i, 2);
            haltonYCoords_[i] = halton(i, 3);
        }
    }

    printf("Rendering...\n");

    pixelsRendered_ = 0;

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
    xPixel_ = imagePlaneX * (camera.width / (double)width_);
    yPixel_ = imagePlaneY * (camera.height / (double)height_);
	
	//Set up the information for threads
	renderingThreads_ = new std::thread[noThreads];
	threadStartPixels_ = new int[noThreads];
	threadEndPixels_ = new int[noThreads];
	
	//Pieces of picture each thread must render
	int noPixels = height_ * width_;
	if (noThreads == 1) {
		threadStartPixels_[0] = 0;
		threadEndPixels_[0] = noPixels;
	} else {
		for (int i = 0; i < noThreads - 1; i++) {
			threadStartPixels_[i] = noPixels / noThreads * i;
			threadEndPixels_[i] = noPixels / noThreads * (i + 1);
		}
		threadStartPixels_[noThreads-1] = threadEndPixels_[noThreads-2];
		threadEndPixels_[noThreads-1] = noPixels;
	}
	
	for (int i = 1; i < noThreads; i++) {
		renderingThreads_[i] = std::thread(&Scene::threadDoWork, this, i);
	}
	
	threadDoWork(0);
	
	for (int i = 1; i < noThreads; i++) {
		renderingThreads_[i].join();
	}
    
    //Dispose of the Halton cache
    if (samplingMode == HALTON) {
        delete (haltonXCoords_);
        delete (haltonYCoords_);
    }

    printf("Rendering complete, postprocessing...\n");
    if (postProcess) rendered_->foreach(postProcess);

    printf("Saving...\n");
    //Save the resulting bitmap to file
    rendered_->saveToFile(filename);

    printf("Done. kd-tree visits per pixel: %f\n", (double)(photonMap_->kdTreeVisited_) / (double)(width_ * height_));
}
