#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include "Vector.h"

#include <queue>
#include <iostream>
#include <fstream>

using namespace std;

//Used for photon mapping, stores information about photons on the scene
struct Photon {
    Vector position;
    Vector direction;
    Vector energy;
    char axis; //Used to store in a kd tree, axis where the subdivision happened
    
    //Precomputed irradiance estimates and the normal of the surface for final gather
    //acceleration
    Vector normal;
    Vector irradiance;
};

//Is put into the neighbours priority queue during the kNN search
struct Neighbour {
    int id;         //Position in the photons_ array
    double distance;//Distance to the point we are looking for
    bool operator< (const Neighbour& n) const {
        return (distance < n.distance);
    }
};

class PhotonMap {
private:
    Photon* photons_;
    int currPtr_;
    int kdTreeSize_;

    int* kdTree_;

    int irradiancePhotonFrequency_;

    int countPhotonsAt(Vector point, double distsq);
    double getDistance(Vector point, int noPhotons);

    void balanceTree(int left, int right, int depth, int location);
    double getComponent(Vector a, char axis);
    char getGreatestSpreadAxis(int left, int right);
    
    void addNearestNeighbour(int newPh, double newPhDist);
    void findNearestNeighbours(Vector point, int treePos);
    void nearestNeighboursWrapper(Vector point, int amount);
    Vector irradianceEstimate(Vector point, Vector normal, int noPhotons);
    
    void findIrradiancePhoton(Vector point, Vector normal, double threshold, int treePos);

    void dumpTree();
    void dumpList();

    int neighboursNeeded_;
    priority_queue<Neighbour> neighbours_;

    int irradiancePhotonId_;
    double irradiancePhotonDist_;

public:
    PhotonMap(int size);
    PhotonMap(char* path);
    ~PhotonMap();
    
    unsigned long long kdTreeVisited_;

    void addPhoton(Vector position, Vector direction, Vector energy, Vector normal);
    void makeTree();
    Vector acceleratedIrradiance(Vector point, Vector normal, double threshold);

    void scalePhotonPower(double factor);

    void precalculateIrradiance(int frequency, int noPhotons);
    Vector visualizePhoton(Vector point, double weight);

    void saveToFile(char* path);
};

#endif
