### Package components

#### networkexposure (.c, .h)
Component to compute Network Exposure Model

#### affiliationexposure (.c, .h)
Component to compute Affiliation Exposure Model

#### main.c
main.c demonstrates how to wield networkexposure.c and affiliationexposure.c and also perform an evaluation of its performance to handle large network files. The coding below uses either multi-threading of exposure computation to handle multiple files or just a single thread for one network or a handful of network files. 

See the "main" function and the examples below for further details on executing either the multi-threaded or single threaded execution.


Example for single thread execution of network exposure model

```C
//Define the path for output of the performance results
char fileName [] = "path_to_performance_output.csv";
    
//Define the path location of the network files
char target_directory [] = "path_to_folder_network_file_and_attribute_file";

//N.B. "NETWORK" flags for single threaded computation 
runBatchTest(target_directory, fileName, NETWORK, NULL, true);

return(0);
```
Example for single thread execution of affiliation exposure model

```C
//Define the path for output of the performance results
char fileName [] = "path_to_performance_output.csv";
    
//Define the path location of the network files
char target_directory [] = "path_to_folder_network_file_and_attribute_file";

//N.B. "AFFILIATION" flags for single threaded computation. 
//The character parameter is used to indicate label used to identify nodes that are in the first mode. 
runBatchTest(target_directory, fileName, AFFILIATION, "1",true);

return(0);
```

Example for multi-thread execution of network exposure model

```C
if (pthread_mutex_init(&lock, NULL) != 0)
{
	printf("\n thread initializing failed\n");
	return 1;
}

//Define the path for output of the performance results
char fileName [] = "path_to_performance_output.csv";
    
//Define the path location of the network files
char target_directory [] = "path_to_folder_network_file_and_attribute_file";

//N.B. "NETWORK_THREAD" flags for multi-threaded computation  
runBatchTest(target_directory, fileName, NETWORK_THREAD, NULL, true);

pthread_mutex_destroy(&lock);
pthread_exit(NULL);
```

Example for multi-thread execution of affiliation exposure model

```C
if (pthread_mutex_init(&lock, NULL) != 0)
{
	printf("\n thread initializing failed\n");
	return 1;
}

//Define the path for output of the performance results
char fileName [] = "path_to_performance_output.csv";
    
//Define the path location of the network files
char target_directory [] = "path_to_folder_network_file_and_attribute_file";

//N.B. "AFFILIATION_THREAD" flags for multi-threaded computation  
//The character parameter is used to indicate label used to identify nodes that are in the first mode. 
runBatchTest(target_directory, fileName, AFFILIATION_THREAD, "1",true);

pthread_mutex_destroy(&lock);
pthread_exit(NULL);
```