# FPPQ

This is an implementation of the "Relaxed Parallel Priority Queue with Filter Levels" from our VMV 2022 paper: https://diglib.eg.org/handle/10.2312/vmv20221202. 

There is a general implementation of the Priority Queue in FPPQ.h, which contains the changeable parameters for the filter levels.
Additionally, there is an implementation used for mesh decimation in FPPQMesh.h and the framework for mesh decimation in most of the other files. There are implementations for a Multi Queue and  a single threaded Priority Queue that can be used for mesh decimation as well.

The main is contained in MeshSimplification.cpp and starts a benchmark for mesh decimation with parameters defined in defines.h.

In Examples.h are examples of how to use the normal filtered priority queue and the one for mesh decimation.

Important notes:
- The quality metrics in MeshMetrics.h do not work out of the box. 
    - We used a library for the hausdorff distance, which is not included. 
    - For The Mesh Similarity one needs a single threaded decimated version of the mesh (generated with this implementation to ensure same decimation steps) for comparison which has to be called: <br> 
&emsp; \<ObjName\>\<DecimationPercent\>%.obj 
- The Priority Queue pop() for mesh decimation will not terminate if there are no valid options (because of topology checks) to collapse.
