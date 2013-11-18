#include "PhotonMap.h"
#include <limits>
#define DINFINITY std::numeric_limits<double>::infinity()

//for printf debugging
#include <cstdio>

#include <cstdlib>

//Determines which axis we are currently using to compare photons
char currAxis;
//Has to be global because it is accessed from a comparator routine.

inline double sqr(double a) {return a*a;}

double getVectorComponent(Vector a, char axis) {
    switch(axis) {
        case 0: return a.getX(); break;
        case 1: return a.getY(); break;
        case 2: return a.getZ(); break;
    }
}

PhotonMap::PhotonMap(int size) {
    printf("Allocating the space for the map (%d bytes per photon)...\n", sizeof(Photon));
    photons_ = new Photon[size];
    currPtr_ = 0;

    kdTreeVisited_ = 0;
}

//Load the map from file
PhotonMap* PhotonMap::makeFromFile(char* path) {
    ifstream file(path, ios::in | ios::binary);

    int mapLength;
    int mapPhoFreq;

    file.read((char*)&mapLength, sizeof(int));
    file.read((char*)&mapPhoFreq, sizeof(int));

    //Find the size of the kd-tree
    int mapkdTreeSize = 1;
    while (mapkdTreeSize < mapLength + 1) mapkdTreeSize = mapkdTreeSize << 1;

    //Check the file size
    file.seekg(0, ios_base::end);
    int requiredSize = (2 + mapkdTreeSize) * sizeof(int) + sizeof(Photon) * mapLength;

	int length = file.tellg();
    if (length != requiredSize) {
        printf("Error: invalid photon map size! The photon map will be regenerated\n");
        return NULL;
    }

    //Allocate the space for the new photon map and populate it
    PhotonMap* map = new PhotonMap();
    map->currPtr_ = mapLength;
    map->irradiancePhotonFrequency_ = mapPhoFreq;
    map->kdTreeSize_ = mapkdTreeSize;

    //Get back to the original position
    file.seekg(2 * sizeof(int), ios_base::beg);
    
    //Allocate the space for the data structures
    map->photons_ = new Photon[mapLength];
    map->kdTree_ = new int[mapkdTreeSize];
    
    file.read((char*)(map->photons_), sizeof(Photon)*mapLength); 
    file.read((char*)(map->kdTree_), sizeof(int)*mapkdTreeSize);

    file.close();

    return map;
}

//Dumps the photon map into a file for further loading
void PhotonMap::saveToFile(char* path) {
    ofstream file(path, ios::out | ios::binary);

    //Format: number of photons in the map, then the irradiance photon frequency,
    //then number of photons 128-byte records, then (nearest power of 2 higher than the number
    //of photons) kd-tree entries (integers).
    file.write((char*)&currPtr_, sizeof(int));
    file.write((char*)&irradiancePhotonFrequency_, sizeof(int));

    file.write((char*)photons_, sizeof(Photon)*currPtr_);
    file.write((char*)kdTree_, sizeof(int)*kdTreeSize_);
    
    file.close();
}

PhotonMap::~PhotonMap() {
    delete(photons_);
    delete(kdTree_);
}

//The comparator used to determine which photon comes earlier on some axis
int photonComparator (const void* p1, const void* p2) {
    Photon *photon1 = (Photon*) p1;
    Photon *photon2 = (Photon*) p2;
    return(getVectorComponent(photon1->position, currAxis) -
           getVectorComponent(photon2->position, currAxis)) > 0 ? 1 : -1;
}

//Returns the axis with the largest min-max photon position spread
char PhotonMap::getGreatestSpreadAxis(int left, int right) {
    double greatestSpread = 0;
    char greatestSpreadAxis = 0;

    double min = DINFINITY;
    double max = 0;

    for (int i = 0; i < 3; i++) {
        min = DINFINITY;
        max = 0;

        for (int j = left; j < right; j++) {
            double position = getVectorComponent(photons_[j].position, i);
            if (position < min) min = position;
            if (position > max) max = position;
        }

        if (max - min > greatestSpread) {
            greatestSpread = max - min;
            greatestSpreadAxis = i;
        }
    }

    return greatestSpreadAxis;
}

//Populates the kdTree_ array with integer pointers to photons that form a heap
//(2n element and 2n+1st element are the children of the nth element)
void PhotonMap::balanceTree(int left, int right, int depth, int position) {
    //This case is weird.
    if (right == left) return;
    
    //One element remaining, record it
    if (right - left == 1) {
        kdTree_[position] = left;
        return;
    }
    
    //Make sure we get the most out of our subdivision
    currAxis = getGreatestSpreadAxis(left, right);
    
    //Find the median of the current sector (there is a faster way, but building the
    //tree is quite fast as it is)
    qsort(&(photons_[left]), right - left, sizeof(Photon), &photonComparator);
    int median = (left + right) / 2;
    photons_[median].axis = currAxis;

    //Record the node of the tree and balance its children
    kdTree_[position] = median;
   
    balanceTree(left, median, depth + 1, position * 2);
    balanceTree(median + 1, right, depth + 1, position * 2 + 1);
}

//Prints a photon for debugging purposes
void dumpPhoton(Photon p) {
 printf("photon at %f, %f, %f from %f, %f, %f with %f, %f, %f, axis ",
        p.position.getX(), p.position.getY(), p.position.getZ(),
        p.direction.getX(), p.direction.getY(), p.direction.getZ(),
        p.energy.getX(), p.energy.getY(), p.energy.getZ());
    if (p.axis == 0) printf("X\n");
    if (p.axis == 1) printf("Y\n");
    if (p.axis == 2) printf("Z\n");
}

//Outputs all photons for debugging purposes
void PhotonMap::dumpList() {
    for (int i = 0; i < currPtr_; i++) {
        printf("%d ", i);
        dumpPhoton(photons_[i]);
    }
}

//Outputs the photons in the way they are ordered in the kd-tree
void PhotonMap::dumpTree() {
    for (int i = 1; i < kdTreeSize_; i++) {
        printf("%d ", i);
        if (kdTree_[i] == -1) printf("NULL\n");
        else dumpPhoton(photons_[kdTree_[i]]);
    }
}


void PhotonMap::makeTree() {
    printf("Making the kd-tree for %d photons...\n", currPtr_);
    
    //Allocate the required array for the tree (next power of two greater than
    //the number of photons)
    kdTreeSize_ = 1;
    while (kdTreeSize_ < currPtr_ + 1) kdTreeSize_ = kdTreeSize_ << 1;
    printf("Size %d\n", kdTreeSize_);
    kdTree_ = new int[kdTreeSize_];

    //-1 denotes that a node is empty (only happens in the last half of the tree, as
    //it is balanced
    for (int i = 0; i < kdTreeSize_; i++) kdTree_[i] = -1;

    balanceTree(0, currPtr_, 0, 1);
}

void PhotonMap::scalePhotonPower(double factor) {
    for (int i = 0; i < currPtr_; i++) {
        photons_[i].energy *= factor;
    }
}

//Avoids the square root computation
inline double sqDist(Vector a, Vector b) {
    Vector diff = b - a;
    return diff.dot(diff);
}

//Populates the neighbours heap with amount photons closest to the point
void PhotonMap::findNearestNeighbours(Vector point, int treePos, int amount, priority_queue<Neighbour> &neighbours) {
    kdTreeVisited_++;

    //Are we at an empty node?
    if (treePos >= kdTreeSize_) return;
    if (kdTree_[treePos] == -1) return;

    //Try and add ourselves to the heap
	int newPh = kdTree_[treePos];
	double newPhDist = sqDist(point, photons_[kdTree_[treePos]].position);
    if (neighbours.size() < amount) {
        Neighbour newNeighbour;
        newNeighbour.id = newPh;
        newNeighbour.distance = newPhDist;

        neighbours.push(newNeighbour);
    } else {
		//What's the furthest we have?
		double maxDist = neighbours.top().distance;

		//If we can do better, remove that neighbour and add a new one
		if (maxDist > newPhDist) {
			Neighbour newNeighbour;
			newNeighbour.id = newPh;
			newNeighbour.distance = newPhDist;

			neighbours.pop();
			neighbours.push(newNeighbour);
		}
	}

    //The coordinate of the median
    double subdivideLocation = getVectorComponent(photons_[kdTree_[treePos]].position,
                                                  photons_[kdTree_[treePos]].axis);

    //The signed distance to the median
    double distToMedian = 
        getVectorComponent(point, photons_[kdTree_[treePos]].axis) - subdivideLocation;
    
    //Is the point to the left of the median?
    if (distToMedian < 0) {
        findNearestNeighbours(point, 2*treePos, amount, neighbours);
    } else {
        findNearestNeighbours(point, 2*treePos+1, amount, neighbours);
    }
    
    double sqDistToMedian = sqr(distToMedian);
    
    //Can there be a neighbour in the other branch?
    if (sqDistToMedian < neighbours.top().distance) {
        if (distToMedian < 0) {
            findNearestNeighbours(point, 2*treePos+1, amount, neighbours);
        } else {
            findNearestNeighbours(point, 2*treePos, amount, neighbours);
        }
    }
}

//Performs the correct call to findNearestNeighbours
priority_queue<Neighbour> PhotonMap::nearestNeighboursWrapper(Vector point, int amount) {
    priority_queue<Neighbour> result = priority_queue<Neighbour>();

    findNearestNeighbours(point, 1, amount, result);
	
	return result;
}

void PhotonMap::addPhoton(Vector position, Vector direction, Vector energy, Vector normal) {
    Photon p;
    p.position = position; p.direction = direction; p.energy = energy; p.normal = normal;
    photons_[currPtr_++] = p;
}

//Filters the photons based on their distance from the sought point
double simpsonKernel(double sqx) {
    return 3.0 / PI * sqr(1 - sqx);
}

void PhotonMap::precalculateIrradiance(int frequency, int noPhotons) {
    irradiancePhotonFrequency_ = frequency;

    //Precompute the radiance for every frequencyth photon
    for (int i = 0; i < currPtr_; i+= frequency) {
        if (i % 1000 == 0) {
            printf ("\r%d out of %d photons processed (%f\%)\n", 
                i, currPtr_, (double)i/currPtr_*100);
        }
        photons_[i].irradiance = irradianceEstimate(
            photons_[i].position, photons_[i].normal, noPhotons);
    }
}

//Mostly copied from findNearestNeighbours, but only looks for one photon and certain
//conditions
void PhotonMap::findIrradiancePhoton(Vector point, Vector normal,
                                     double threshold, int treePos,
									 double &bestDist, int &bestId) {
     
    if (treePos >= kdTreeSize_) return;
    if (kdTree_[treePos] == -1) return;

    //Set the photon as best if it's indeed an irradiance photon, its normal is close to
    //our needed normal and it's closer to the sought point than the previous best
    if (kdTree_[treePos] % irradiancePhotonFrequency_ == 0) {
        if (normal.dot(photons_[kdTree_[treePos]].normal) >= threshold) {
            double sqdistance = sqDist(point, photons_[kdTree_[treePos]].position);
            if (sqdistance < bestDist) {
				bestDist = sqdistance;
                bestId = kdTree_[treePos];
            }
        }
    }
    
    double subdivideLocation = getVectorComponent(photons_[kdTree_[treePos]].position,
                                                  photons_[kdTree_[treePos]].axis);

    double distToMedian = 
        getVectorComponent(point, photons_[kdTree_[treePos]].axis) - subdivideLocation;
    
    if (distToMedian < 0) {
        findIrradiancePhoton(point, normal, threshold, 2*treePos, bestDist, bestId);
    } else {
        findIrradiancePhoton(point, normal, threshold, 2*treePos + 1, bestDist, bestId);
    }
    
    double sqDistToMedian = sqr(distToMedian);
    
    if (sqDistToMedian < bestDist) {
        if (distToMedian < 0) {
            findIrradiancePhoton(point, normal, threshold, 2*treePos + 1, bestDist, bestId);
        } else {
            findIrradiancePhoton(point, normal, threshold, 2*treePos, bestDist, bestId);
        }
    }
}

//Looks up the closest irradiance photon (wrapper for findIrradiancePhoton)
Vector PhotonMap::acceleratedIrradiance(Vector point, Vector normal, double threshold) {
    double irradiancePhotonDist = DINFINITY;
    int irradiancePhotonId = -1;

    findIrradiancePhoton(point, normal, threshold, 1, irradiancePhotonDist, irradiancePhotonId);
    
    //This shouldn't happen.
    if (irradiancePhotonId == -1) return Vector(0, 0, 0);

    return photons_[irradiancePhotonId].irradiance;
}

//Gathers the photons in a given radius to determine the illumination of an entity at a certain point
//and a certain normal.
Vector PhotonMap::irradianceEstimate(Vector point, Vector normal, int noPhotons) {
    Vector result(0, 0, 0);
    priority_queue<Neighbour> neighbours = nearestNeighboursWrapper(point, noPhotons);

    double sqRadius = neighbours.top().distance;
    
    double factor = 1.0 / (sqRadius);
    while (!neighbours.empty()) { 
        Neighbour neighbour = neighbours.top();
        neighbours.pop();
		
		double dot = normal.dot(-photons_[neighbour.id].direction);
		if (dot < 0.0) continue;
        double weight = simpsonKernel(neighbour.distance / sqRadius) * dot;
		
        result += photons_[neighbour.id].energy * weight;

    }
    return result * factor;
}

//Returns a pixel on the gradient from black to white depending on the distance of the closest
//photon. The smaller the weight, the faster the falloff.
//Radiance photons are highlighted in red
Vector PhotonMap::visualizePhoton(Vector point, double weight) {
    priority_queue<Neighbour> neighbours = nearestNeighboursWrapper(point, 1);
    double distance = neighbours.top().distance;
    double factor = weight / (weight + distance);
	
    if (neighbours.top().id % irradiancePhotonFrequency_ == 0) {
        return Vector(factor, 0, factor);
    } else {
        return Vector(factor, factor, factor);
    }
}
