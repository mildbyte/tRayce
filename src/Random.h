#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <cstdint>

//Generates a random number between 0.0 and 1.0
double drand();

//Seeds the RNG with a value
void seed_drand(uint64_t seed);

#endif
