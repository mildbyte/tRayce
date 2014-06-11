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
    
    PhotonMap() {};

    int countPhotonsAt(Vector point, double distsq);
    double getDistance(Vector point, int noPhotons);

    void balanceTree(int left, int right, int depth, int location);
    double getComponent(Vector a, char axis);
    char getGreatestSpreadAxis(int left, int right);

    void findNearestNeighbours(Vector point, int treePos, int amount, priority_queue<Neighbour> &neighbours);
    priority_queue<Neighbour> nearestNeighboursWrapper(Vector point, int amount);
    
    void findIrradiancePhoton(Vector point, Vector normal, double threshold, int treePos, double &bestDist, int &bestId);

    void dumpTree();
    void dumpList();

public:
    PhotonMap(int size);
    static PhotonMap* makeFromFile(char* path);

    ~PhotonMap();
    
    unsigned long long kdTreeVisited_;

    void addPhoton(Vector position, Vector direction, Vector energy, Vector normal);
    void makeTree();
    Vector acceleratedIrradiance(Vector point, Vector normal, double threshold);
    Vector irradianceEstimate(Vector point, Vector normal, int noPhotons);

    void scalePhotonPower(double factor);

    void precalculateIrradiance(int frequency, int noPhotons);
    Vector visualizePhoton(Vector point, double weight);

    void saveToFile(char* path);
};

#endif
