#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include "Vector.h"

//Used for photon mapping, stores information about photons on the scene

struct Photon {
    Vector position;
    Vector direction;
    Vector energy;
};

class PhotonMap {
private:
    Photon* photons_;
    int currPtr_;

    int countPhotonsAt(Vector point, double distsq);
    double getDistance(Vector point, int noPhotons);
public:
    PhotonMap(int size);
    ~PhotonMap();
    
    void addPhoton(Vector position, Vector direction, Vector energy);
    Vector gatherPhotons(Vector point, Vector normal, double radius, double exposure, int noPhotons);
};

#endif
