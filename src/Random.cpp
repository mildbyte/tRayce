#include "Random.h"

std::mt19937 twister;
std::uniform_real_distribution<double> dist(0, 1);

//Generates a random number between 0.0 and 1.0
double drand() {
    return dist(twister);
}

void seed_drand(uint64_t seed) { twister.seed(seed); }