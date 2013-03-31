#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include "Vector.h"
//Used for photon mapping, stores information about photons on the scene

struct Photon {
    Vector position;
    Vector direction;
    Vector energy;
    char axis; //Used to store in a kd tree, axis where the subdivision happened
};

class PhotonMap {
private:
    Photon* photons_;
    int currPtr_;

    Photon* kdTree_;

    int countPhotonsAt(Vector point, double distsq);
    double getDistance(Vector point, int noPhotons);

    void balanceTree(int left, int right, int depth, int location);
    double getComponent(Vector a, char axis);
    
    void replaceMaxDist(Photon newPh, double newPhDist);
    void findNearestNeighbours(Vector point, int treePos);
    void nearestNeighboursWrapper(Vector point, int amount);

    int neighboursNeeded_;
    int foundNeighbours_;
    Photon* neighbours_;
    double* neighbourDists_;

public:
    PhotonMap(int size);
    ~PhotonMap();
    
    void addPhoton(Vector position, Vector direction, Vector energy);
    Vector gatherPhotons(Vector point, Vector normal, double radius, double exposure, int noPhotons);
    void makeTree();
};

#endif
