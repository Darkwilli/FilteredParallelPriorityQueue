#pragma once
#include <string>

const std::string MESH_PATH = "D:/models/lucy.obj";
const std::string OUTPUT_DIRECTORY = "D:/PPQOutput";
const unsigned int DESIRED_NUMBER_OF_VERTICES = 14000;

constexpr int slotSize = 8;
const int NUM_THREADS = 32;

constexpr int countIterations = 100;

//#define MULTI_QUEUE
//#define SINGLE_THREADED
#define NEXT_BEST_ERROR // If the collapse is invalid (topology check) the next best vertex to collapse to is searched instead and the PQ updated instead of disallowing HECs from this vertex 
#define TRIANGLE_INDICES // Uses Triangle indices instead of saving every triangle in every vertex (faster)
#define SINGLE_THREADED_LOADING // makes the loading single threaded MT seems to work but was not used for the benchmarks (because loading is not included in the total times)
#define ERROR_UPDATE_OPTIMISED // Checks if the collapsed vertex or the vertex to collapse to is the previous best vertex. Only in that case all errors have to be recalculated to find the new best error otherwise only one  error has to be calculated
#define NORMAL_VALIDATION_CHECK // Checks if the updated face normals deviate too much from the old ones
#define NORMAL_THRESHOLD 0.5 // Dot Product of new Face normal and old Face normal has to be over this for every face 
#define BENCHMARK 1
#define BENCHMARK_DEBUG 1 // have NUM_THREADS as first benchmark and decimate to DESIRED_NUMBER_OF_VERTICES
//#define NO_HAUSDORFF
//#define EXPORT // Export the obj after decimation
#define NO_DEBUG_CHECKS // Disables some debug checks
#define VECTOR_REDUCTION // Use a vector as a map for the final reduction step (is faster but uses more memory)


#ifdef MULTI_QUEUE
#define PRIORITY_STRUCTURE MultiQueue
#else
#ifdef SINGLE_THREADED
#define PRIORITY_STRUCTURE PriorityQueueSingleThreaded
#else
#define PRIORITY_STRUCTURE FPPQMesh
#endif
#endif