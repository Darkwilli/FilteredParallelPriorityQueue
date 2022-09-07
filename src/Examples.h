#pragma once

#include <vector>
#include <iostream>

#include "FPPQ.h"
#include "FPPQMesh.h"
#include "ParallelMesh.h"

void FPPQexample() {
    int maxCapacity = 9999999;
    FPPQ<float> fppq(maxCapacity); // Create a priority Queue with float errors and a maximum capacity of 9999999 elements

    std::vector<std::pair<float, int>> initialValues; // set up some initial elements for the queue
    for(int i = 0; i < 1000000; i++) {
        initialValues.push_back({ 1.f + sin(float(i) / 100.f),i }); // Setup each element with an inital error
    }
    int countThreads = 4;
    fppq.fill(initialValues, countThreads); // Initialise the queue with 4 Threads (the parallel std::sorting function inside can actually use more than that)

#pragma omp parallel for num_threads(countThreads)
    for(int i = 0; i < 1000000; i++) {
        fppq.update(i, cos(float(i) / 100.f)); // Update the element with a new value
    }

#pragma omp parallel for num_threads(countThreads)
    for(int i = 0; i < 1000000; i++) {
        int elementId = fppq.pop(); // pop an element
        // Do something with the element

        fppq.insert(i + 1000000, sin(float(i + elementId) / 100.f)); // Insert a new element
        int elementId2 = fppq.pop(); // pop an other element
        // Do something with the other element
    }
}

void meshDecimationExample() {
    auto& pMesh = ParallelMesh::getInstance(); // get the mesh
    if(!pMesh.loadFromObj(std::filesystem::path(MESH_PATH))) { // Load an obj
        std::cout << "Could not load mesh: "
            << MESH_PATH << std::endl;
        return;
    }


    int countThreads = 4;
    float decimation = .01; // Decimation to 1% of original vertices
    pMesh.computeQuadricErrorMatrices(countThreads); // Initialise the error quadrics with 4 threads (the parallel std::sorting function inside can actually use more than that)
    int finalVertexCount = pMesh.reduceVerticesTo(decimation, countThreads); // Reduce the mesh with 4 threads to 1%
    pMesh.reduceVertices(finalVertexCount, countThreads); // Perform the deletion of the Collapsed vertices and repair the mesh structure (vertex ids etc.)
    pMesh.deletePriorityStructure(); // Delete the Priority Queue
}