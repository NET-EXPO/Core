### Package components

#### networkexposure (.c, .h)
Component to compute Network Exposure Model

#### affiliationexposure (.c, .h)
Component to compute Affiliation Exposure Model

#### main.c
main.c demonstrates how to wield networkexposure.c and affiliationexposure.c and also perform an evaluation of its performance to handle large network files. The coding below uses either multi-threading of exposure computation to handle multiple files or just a single thread for one network or a handful of network files. 

See the "main" function for further details on executing either the multi-threaded or single threaded execution.

Example for single thread execution of network exposure model

```C
//Define the path for output of the performance results
char fileName [] = "path_to_performance_output.csv";
    
//Define the path location of the network files
char target_directory [] = "path_to_folder_network_file_and_attribute_file";

runBatchTest(target_directory, fileName, NETWORK, NULL, true);

return(0);
```
