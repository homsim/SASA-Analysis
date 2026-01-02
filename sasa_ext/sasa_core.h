#ifndef SASA_CORE_H
#define SASA_CORE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constants
#define VMD_SEED 38572111  // Match VMD's fixed seed for reproducibility
#define PI 3.14159265358979323846

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    int *neighbors;
    int count;
    int capacity;
} NeighborList;

typedef struct {
    Vec3 *points;
    int count;
    int capacity;
} PointList;

// Core SASA computation
int compute_sasa(
    const float *coords,      // atom coordinates [3*n_atoms]
    const float *radii,       // atom radii [n_atoms]
    int n_atoms,
    float probe_radius,       // srad parameter
    int n_samples,           // samples parameter
    unsigned int seed,       // for reproducibility
    float *total_sasa,       // output: total SASA area
    PointList *surface_points // output: surface point coordinates
);

// Neighbor list construction
int build_neighbor_lists(
    const float *coords,
    const float *radii,
    int n_atoms,
    float probe_radius,
    NeighborList **neighbor_lists
);

// Random sphere point generation
void generate_sphere_points(
    Vec3 *sphere_points,
    int n_samples,
    unsigned int seed
);

// Point burial test
int is_point_buried(
    Vec3 point,
    const float *coords,
    const float *radii,
    const int *neighbors,
    int n_neighbors,
    float probe_radius
);

// Memory management
void free_neighbor_lists(NeighborList *lists, int n_atoms);
void free_point_list(PointList *list);

#endif