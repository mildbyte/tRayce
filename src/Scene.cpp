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

void Scene::addTriangle(Triangle* triangle) {
    trianglesVector_.push_back(triangle);
}

void Scene::addRenderable(Renderable* renderable) {
    renderables_.addRenderable(renderable);
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
    camera.lensRadius = 0.0;
    camera.focalDistance = 10.0;

    backgroundColor.set(0, 0, 0);
    
    pathTracingMaxDepth = 5;
    pathTracingSamplesPerPixel = 10;
	pathTracingTerminationProbability = 0.5;

    samplingMode = STRATIFIED;
}

//Adds a model inside an OBJ file into the scene with a given
//material. Does not reuse vertices, so each triangle takes
//3xvertexsize + overhead space. Only supports triangles.
void Scene::importObj(char* filename, Material m, Vector shift, double scale) {
	vector<Vector> vectors;
	vector<Vector> normalVectors;
	vector<Vector> textureVectors;

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
		}
		else if (type == "vn") {
			double x, y, z;
			stream >> x >> y >> z;
			normalVectors.push_back(Vector(x, y, z));
		}
		else if (type == "vt") {
			double u, v, w = 0;

			string vertices;
			stream >> vertices;

			//u, v are mandatory, whereas w is optional
			istringstream ss(vertices);
			ss >> u >> v;

			if (!ss.eof()) ss >> w;

			textureVectors.push_back(Vector(u, v, w));
		}
		else if (type == "f") {
			faces++;
			//Face format: f vertex vertex vertex
			//vertex: v/t/n
			//where v is the vertex id, t is the texture UV coordinate id
			//n is the normal direction at the vertex
			//Any of the second or the third can be skipped:
			//v, v/t, v//n are all valid

			int vIds[3];
			int normalvIds[3] = { -1, -1, -1 };
			int textureIds[3] = { -1, -1, -1 };

			for (int i = 0; i < 3; i++) {
				string vertices;
				stream >> vertices;

				//Split on '/'
				istringstream ss(vertices);
				string tok;
				getline(ss, tok, '/'); //First number: vertex id
				vIds[i] = atoi(tok.c_str());

				if (!ss.eof()) {
					getline(ss, tok, '/'); //Second number: texture UV coordinates
					if (tok != "") textureIds[i] = atoi(tok.c_str());
				}

				if (!ss.eof()) {
					getline(ss, tok, '/'); //Last number: normal vertex id
					normalvIds[i] = atoi(tok.c_str());
				}
			}

			Triangle *t;

			//Vertex ids in OBJ are 1-based (so 0 never appears in a face description)
			if (normalvIds[0] != -1 && normalvIds[1] != -1 && normalvIds[2] != -1) {
				t = new Triangle(vectors[vIds[0] - 1], vectors[vIds[1] - 1], vectors[vIds[2] - 1],
					normalVectors[normalvIds[0] - 1], normalVectors[normalvIds[1] - 1], normalVectors[normalvIds[2] - 1]);
			}
			else {
				t = new Triangle(vectors[vIds[0] - 1], vectors[vIds[1] - 1], vectors[vIds[2] - 1]);
			}

			if (textureIds[0] != -1 && textureIds[1] != -1 && textureIds[2] != -1) {
				t->setTextureCoords(textureVectors[textureIds[0]], textureVectors[textureIds[1]], textureVectors[textureIds[2]]);
			}

			t->material = m;
			addTriangle(t);
			//addRenderable(t);
		}
		else {
			stream.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
	}

	cout << "Loaded " << vectors.size() << " vertices, " << normalVectors.size() <<
		" normal vectors and " << faces << " faces." << endl;

	stream.close();
}

void Scene::refractRay(Intersection inter, Ray &ray, double &reflectance) {
	//Refracts a ray, using Fresnel equations to calculate the fraction that should also be reflected.

    double c1 = inter.normal.dot(ray.direction);
	double n_in, n_out;

    //Make c1 positive. If it was positive, the ray hit the sphere from the inside
    //and so the normal must be inverted.
    if (c1 < 0) {
        c1 = -c1;
		n_in = 1;
		n_out = inter.object->material.refrIndex;
    } else {
        inter.normal = -inter.normal;
		n_in = inter.object->material.refrIndex;
		n_out = 1;
    }

	double refrRatio = n_in / n_out;

    double c2 = 1 - pow(n_in / n_out, 2) * (1 - c1 * c1);

    //Bail out if the cosine of the angle squared is > 1 (total internal
    //reflection) or negative (cannot calculate the square root)

    if ((c2 > 1) || (c2 < 0)) {
		reflectance = 1;
		return;
    }

    c2 = sqrt(c2);

	double r_s = pow((n_in * c1 - n_out * c2) / (n_in * c1 + n_out * c2), 2);
	double r_p = pow((n_out * c1 - n_in * c2) / (n_out * c1 + n_in * c2), 2);
	reflectance = (r_s + r_p) * 0.5;

    //Calculate the direction of the refracted ray
    Vector refrDir = ray.direction * refrRatio + inter.normal * (refrRatio * c1 - c2);
    refrDir.normalize();

    ray.origin = inter.coords;
    ray.direction = refrDir;
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

Vector sampleHemisphere(Vector normal) {
    while (true) {
        Vector direction(2*drand()-1, 2*drand()-1, 2*drand()-1);
        
        double dot = normal.dot(direction);
        if (dot >= 0.0 && dot <= 1.0) {
            direction.normalize();
            return direction;
        }
    }
}

Vector sampleHemisphereCosine(Vector normal, double Xi1, double Xi2) {
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

Vector Scene::getColorAt(Vector point) {
	double thresh = 0.05;
	if (abs(point[0] - round(point[0])) < thresh || abs(point[2] - round(point[2])) < thresh) return Vector(0, 0, 0);
	return Vector(1, 1, 1);
}

Vector Scene::pathTrace(const Ray ray, int depth, double &dist) {
    //If the maximum tracing depth has been reached, return
	//(do Russian Roulette for path termination)
	if (depth <= 0) {
		double decision = drand();
		if (decision < pathTracingTerminationProbability) return Vector(0, 0, 0);
	}

    //If we are casting a ray from the eye into the scene, cull collisions behind the image plane
    Intersection inter = getClosestIntersection(ray, depth == pathTracingMaxDepth ? camera.planeDistance : 0);
    
    if (!inter.happened) {
		dist = INFINITY;
        return backgroundColor;
    }

    inter.normal.normalize();
	dist = inter.distance;

    double nonDiffuse = inter.object->material.reflectivity + inter.object->material.transparency;

    double decision = drand();
    
    // The material can be part diffuse and part either perfectly specular or transparent.
	// For specular, we reflect the ray and trace it; for transparent, use Fresnel equations
	// to compute the weights of the reflected and the transmitted components.
	double dummy;
	Vector result;

	if (decision > nonDiffuse) {
		Ray nextRay;
		nextRay.origin = inter.coords;
		nextRay.direction = sampleHemisphereCosine(inter.normal, drand(), drand());

		Vector color = inter.object->material.color;
		//See if we can texture the object

		if (inter.object->material.isTextured) {
			double u;
			double v;

			if (inter.object->getUVAt(inter.coords, &u, &v)) {
				color = inter.object->material.texture->textureSample(u, v);
			}
		}

		Vector brdf = color * nextRay.direction.dot(inter.normal);

		Vector reflected = pathTrace(nextRay, depth - 1, dummy);

		result = inter.object->material.emittance + combineColors(brdf, reflected);
	} else if (inter.object->material.transparency > 0) {
		// Transparent material, use Fresnel's equations to compute the reflectance
		double reflectance;

		Ray nextRay = ray;
		refractRay(inter, nextRay, reflectance);

		decision = drand();

		if (decision > reflectance) {
			result = inter.object->material.emittance
				+ combineColors(inter.object->material.color, pathTrace(nextRay, depth - 1, dummy));
		} else {
			if (inter.normal.dot(ray.direction) > 0) inter.normal = -inter.normal;
			Ray nextRay = reflectRay(inter, ray);
			result = inter.object->material.emittance
				+ combineColors(inter.object->material.color, pathTrace(nextRay, depth - 1, dummy));
		}
	} else {
		// Specular reflection
		Ray nextRay = reflectRay(inter, ray);
		result = inter.object->material.emittance
			+ combineColors(inter.object->material.color, pathTrace(nextRay, depth - 1, dummy));
    }

	//The estimate will be divided by the termination probabiliy to compensate for the RR
	if (depth <= 0) return result * (1.0 / (1.0 - pathTracingTerminationProbability));
	else return result;
}

Intersection Scene::getClosestIntersection(Ray ray, double cullDistance) {
	Intersection inter = renderables_.getFirstIntersection(ray, cullDistance);

	Intersection inter2 = triangles_->getFirstIntersection(ray, cullDistance);
	if (!inter.happened || (inter2.happened && inter.distance > inter2.distance)) inter = inter2;

	return inter;
}

//Converts screen to image coordinates and traces them.
Vector Scene::tracePixel(double x, double y, double &dist) {    
	dist = 0;

    int totalSamples = pathTracingSamplesPerPixel * pathTracingSamplesPerPixel;
    
    Vector result(0, 0, 0);
	for (int i = 0; i < totalSamples; i++) {
		Ray ray;

		//Exploit multiple SPP to perform AA
		double xShift = drand() - 0.5;
		double yShift = drand() - 0.5;

		Vector xWorld = xPixel_ * (x - 0.5 * (double)width_ + xShift);
		Vector yWorld = yPixel_ * (y - 0.5 * (double)height_ + yShift);
		ray.direction = (camera.direction * camera.planeDistance)
			+ xWorld + yWorld;
		ray.direction.normalize();

		Vector pixelPoint = camera.position + ray.direction * camera.focalDistance;

		//Similar to the raytracing case, but the ray origin
		//is perturbed
		double r = drand() * camera.lensRadius;
		double theta = drand() * 2 * PI;

		ray.origin = camera.position + imagePlaneX_ * r * cos(theta)
			+ imagePlaneY_ * r * sin(theta);

		//Alter the direction so that the ray goes through the same pixel in the
		//focal plane
		ray.direction = pixelPoint - camera.position;
		ray.direction.normalize();

		double d;
		result += pathTrace(ray, pathTracingMaxDepth, d);
		dist += d;
	}

	result /= (double)totalSamples;
	dist /= (double)totalSamples;

	return result;
}

void Scene::threadDoWork(int threadId, int noThreads) {
	//A rendering task executed in a separate thread. Renders the pixels assigned to it
	//and stores them in the bitmap.
	
    //The resulting color of the pixel
    Vector resultColor;

    //Renderable that has been hit previously
    //Renderable* prevHit = NULL;
    
    int onePercent = width_ * height_ / 1000;
    
    for (int pixel = threadId; pixel < width_ * height_; pixel += noThreads) {    
		int realx = pixel % width_;
		int realy = pixel / width_;

		//Trace a ray through the pixel on the image plane
		double dist;
		resultColor = tracePixel(realx, realy, dist);
		rendered_->setPixel(realx, realy, resultColor, dist);

		pixelsRendered_++;
		if (pixelsRendered_ % onePercent == 0) {
			printf("Rendered %d pixel(s) out of %d (%f%%)\r", pixelsRendered_, 
				(width_ * height_),
				pixelsRendered_ / (double)(width_ * height_) * 100);
            fflush(stdout);
		}
    } 
}

void Scene::render(char* filename, BitmapPixel (*postProcess)(BitmapPixel), int noThreads) {
    //Renders the scene to a file

    //Build a kd-tree for the triangles
    printf("Building a kd-tree for the triangles...\n");
    triangles_ = KDNode::build(trianglesVector_, 0);
    
    printf("Built, %d triangles\n", triangles_->getItems().size());
    printf("Rendering...\n");

    pixelsRendered_ = 0;

    //imagePlaneX and imagePlaneY are the unit x and y vectors in the image plane.
    imagePlaneX_ = Vector(
        camera.direction.getZ(),
        0,
        -camera.direction.getX());
    
	// Negate imagePlaneY so that it increases downwards (since it does in the bitmap)
    imagePlaneY_ = -Vector(
         -camera.direction.getY()*camera.direction.getX(),
         camera.direction.getZ()*camera.direction.getZ() 
            + camera.direction.getX() * camera.direction.getX(),
         -camera.direction.getZ()*camera.direction.getY());

    imagePlaneX_.normalize();
    imagePlaneY_.normalize();

    //xPixel and yPixel: one pixel in the image plane
    xPixel_ = imagePlaneX_ * (camera.width / (double)width_);
    yPixel_ = imagePlaneY_ * (camera.height / (double)height_);
	
	//Set up the information for threads
	renderingThreads_ = new std::thread[noThreads];
	
	//Start up the threads
	for (int i = 1; i < noThreads; i++) {
		renderingThreads_[i] = std::thread(&Scene::threadDoWork, this, i, noThreads);
	}
	
	threadDoWork(0, noThreads);
	
	for (int i = 1; i < noThreads; i++) {
		renderingThreads_[i].join();
	}

    printf("\nRendering complete, postprocessing...\n");
    if (postProcess) rendered_->foreach(postProcess);

    printf("Saving...\n");
    //Save the resulting bitmap to file
    rendered_->saveToFile(filename);

    printf("Done.\n");
}
