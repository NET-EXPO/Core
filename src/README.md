### Package components

#### networkexposure (.c, .h)
Component to compute Network Exposure Model

#### affiliationexposure (.c, .h)
Component to compute Affiliation Exposure Model

#### main.c
main.c demonstrates how to wield networkexposure.c and affiliationexposure.c and also perform an evaluation of its performance to handle large network files. The coding below uses either multi-threading of exposure computation to handle multiple files or just a single thread for one network or a handful of network files. See the main function for further details on executing either the multi-threaded or single threaded execution.