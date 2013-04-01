#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include "Vector.h"

#include <queue>
//Used for photon mapping, stores information about photons on the scene

struct Photon {
    Vector position;
    Vector direction;
    Vector energy;
    char axis; //Used to store in a kd tree, axis where the subdivision happened
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
    int* indices_;


    int countPhotonsAt(Vector point, double distsq);
    double getDistance(Vector point, int noPhotons);

    void balanceTree(int left, int right, int depth, int location);
    double getComponent(Vector a, char axis);
    char getGreatestSpreadAxis(int left, int right);
    
    void addNearestNeighbour(int newPh, double newPhDist);
    void findNearestNeighbours(Vector point, int treePos);
    void nearestNeighboursWrapper(Vector point, int amount);

    void dumpTree();
    void dumpList();
//    void dumpNeighbours();

    int neighboursNeeded_;
    
    std::priority_queue<Neighbour> neighbours_;

public:
    PhotonMap(int size);
    ~PhotonMap();
    
    int kdTreeVisited_;

    void addPhoton(Vector position, Vector direction, Vector energy);
    Vector gatherPhotons(Vector point, Vector normal, int noPhotons);
    void makeTree();
};

#endif
