#include "sasa_core.h"

// Simple linear congruential generator for reproducible random numbers
typedef struct {
    unsigned int state;
} RNG;

static void rng_init(RNG *rng, unsigned int seed) {
    rng->state = seed;
}

static float rng_uniform(RNG *rng) {
    // LCG constants from Numerical Recipes
    rng->state = (1664525U * rng->state + 1013904223U);
    return (float)(rng->state) / (float)(0xFFFFFFFFU);
}

void generate_sphere_points(Vec3 *sphere_points, int n_samples, unsigned int seed) {
    RNG rng;
    rng_init(&rng, seed);

    for (int i = 0; i < n_samples; i++) {
        float u1 = rng_uniform(&rng);
        float u2 = rng_uniform(&rng);

        // Uniform distribution on sphere using spherical coordinates
        float z = 2.0f * u1 - 1.0f;
        float phi = 2.0f * PI * u2;
        float r = sqrtf(1.0f - z * z);

        sphere_points[i].x = r * cosf(phi);
        sphere_points[i].y = r * sinf(phi);
        sphere_points[i].z = z;
    }
}

int build_neighbor_lists(
    const float *coords,
    const float *radii,
    int n_atoms,
    float probe_radius,
    NeighborList **neighbor_lists
) {
    *neighbor_lists = (NeighborList*)malloc(n_atoms * sizeof(NeighborList));
    if (!*neighbor_lists) return -1;

    // Find maximum radius for cutoff calculation
    float max_radius = 0.0f;
    for (int i = 0; i < n_atoms; i++) {
        if (radii[i] > max_radius) max_radius = radii[i];
    }
    float cutoff = 2.0f * (max_radius + probe_radius);

    for (int i = 0; i < n_atoms; i++) {
        NeighborList *list = &(*neighbor_lists)[i];
        list->capacity = 32;  // Initial capacity
        list->count = 0;
        list->neighbors = (int*)malloc(list->capacity * sizeof(int));
        if (!list->neighbors) return -1;

        Vec3 atom_i = {coords[3*i], coords[3*i+1], coords[3*i+2]};

        for (int j = 0; j < n_atoms; j++) {
            if (i == j) continue;

            Vec3 atom_j = {coords[3*j], coords[3*j+1], coords[3*j+2]};
            float dx = atom_i.x - atom_j.x;
            float dy = atom_i.y - atom_j.y;
            float dz = atom_i.z - atom_j.z;
            float dist = sqrtf(dx*dx + dy*dy + dz*dz);

            if (dist < cutoff) {
                // Resize if needed
                if (list->count >= list->capacity) {
                    list->capacity *= 2;
                    list->neighbors = (int*)realloc(list->neighbors,
                                                   list->capacity * sizeof(int));
                    if (!list->neighbors) return -1;
                }
                list->neighbors[list->count++] = j;
            }
        }
    }

    return 0;
}

int is_point_buried(
    Vec3 point,
    const float *coords,
    const float *radii,
    const int *neighbors,
    int n_neighbors,
    float probe_radius
) {
    for (int i = 0; i < n_neighbors; i++) {
        int neighbor_idx = neighbors[i];
        Vec3 neighbor = {
            coords[3*neighbor_idx],
            coords[3*neighbor_idx+1],
            coords[3*neighbor_idx+2]
        };

        float dx = point.x - neighbor.x;
        float dy = point.y - neighbor.y;
        float dz = point.z - neighbor.z;
        float dist = sqrtf(dx*dx + dy*dy + dz*dz);

        float neighbor_cutoff = radii[neighbor_idx] + probe_radius;
        if (dist < neighbor_cutoff) {
            return 1;  // Point is buried
        }
    }
    return 0;  // Point is exposed
}

int compute_sasa(
    const float *coords,
    const float *radii,
    int n_atoms,
    float probe_radius,
    int n_samples,
    unsigned int seed,
    float *total_sasa,
    PointList *surface_points
) {
    // Initialize surface points list
    surface_points->capacity = n_atoms * n_samples;  // Upper bound
    surface_points->count = 0;
    surface_points->points = (Vec3*)malloc(surface_points->capacity * sizeof(Vec3));
    if (!surface_points->points) return -1;

    // Build neighbor lists for efficiency
    NeighborList *neighbor_lists;
    if (build_neighbor_lists(coords, radii, n_atoms, probe_radius, &neighbor_lists) != 0) {
        return -1;
    }

    // Generate unit sphere points
    Vec3 *sphere_points = (Vec3*)malloc(n_samples * sizeof(Vec3));
    if (!sphere_points) {
        free_neighbor_lists(neighbor_lists, n_atoms);
        return -1;
    }
    generate_sphere_points(sphere_points, n_samples, seed);

    *total_sasa = 0.0f;

    for (int atom_idx = 0; atom_idx < n_atoms; atom_idx++) {
        Vec3 atom_center = {coords[3*atom_idx], coords[3*atom_idx+1], coords[3*atom_idx+2]};
        float atom_radius = radii[atom_idx] + probe_radius;
        int exposed_points = 0;

        for (int point_idx = 0; point_idx < n_samples; point_idx++) {
            // Scale unit sphere point by expanded radius
            Vec3 surface_point = {
                atom_center.x + atom_radius * sphere_points[point_idx].x,
                atom_center.y + atom_radius * sphere_points[point_idx].y,
                atom_center.z + atom_radius * sphere_points[point_idx].z
            };

            // Test if point is buried by neighbors
            if (!is_point_buried(surface_point, coords, radii,
                               neighbor_lists[atom_idx].neighbors,
                               neighbor_lists[atom_idx].count, probe_radius)) {
                exposed_points++;

                // Store surface point
                if (surface_points->count < surface_points->capacity) {
                    surface_points->points[surface_points->count++] = surface_point;
                }
            }
        }

        // Calculate area contribution for this atom
        float atom_area = (4.0f * PI * atom_radius * atom_radius * exposed_points) / n_samples;
        *total_sasa += atom_area;
    }

    free(sphere_points);
    free_neighbor_lists(neighbor_lists, n_atoms);

    return 0;
}

void free_neighbor_lists(NeighborList *lists, int n_atoms) {
    if (lists) {
        for (int i = 0; i < n_atoms; i++) {
            if (lists[i].neighbors) {
                free(lists[i].neighbors);
            }
        }
        free(lists);
    }
}

void free_point_list(PointList *list) {
    if (list && list->points) {
        free(list->points);
        list->points = NULL;
        list->count = 0;
        list->capacity = 0;
    }
}