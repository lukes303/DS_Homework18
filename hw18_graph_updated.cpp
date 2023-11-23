//
// ========================================
// HW18: Implement Graph, Graph Traverse
//       and Dijkstra's algorithm. 
// ========================================
// 
// In this assignment, we will implement an 
// undirected graph, two traverse algorithms 
// and the Dijkstra's algorithm. 
//  
// This template is based on matrix-based 
// graph implementation, but you are free 
// to work on list-based graph implementation. 
// 
// For simplicity, the template assumes no 
// satellite data is used and focuses on 
// dealing with integers (node numbers) only. 
// 
// Finally, feel free to modify the given template 
// including the classes, functions, etc. Just make 
// sure the outputs of functions tested in "main" 
// stay the same (for auto-grading). 
// 

#include <iostream>
#include <stack>
#include <queue>
#include <climits>
using namespace std;

class Result {
public:
    int length; 
    int weight;
    int* path;
    Result(int l, int w, int* p);
};

Result::Result(int l, int w, int* p){
    length = l;
    weight = w;
    path = p;
}

// Please omplete the graph class. 
// Remember this is an undirected graph 
// of "size" nodes labeled from 0 to size-1. 
class Graph {

private:

    // this variable holds the matrix 
    int** m;

    // this variable holds the number of nodes in graph 
    int size;

public:

    // this function returns the degree of node i 
    int Degree(int i);


    // this function adds an edge between node i and 
    // node j, and assigns a weight "w" to the edge 
    // 
    // recall: if we do not want a weighted graph, 
    // simply set w = 1 for all edges 
    //
    // you may also check boundary but in the testing 
    // samples we assume all inputs are within boundary 
    void Add(int i, int j, int w);


    // this function returns 1 if node i and 
    // node j are connected and returns 0 otherwise.
    // 
    // note: you can also let it return the weight;  
    //       weight is zero if no edge exists 
    int IsEdge(int i, int j);


    // this function returns 1 if there is a path 
    // from node i to node j and returns 0 otherwise.  
    int IsPath(int i, int j);

    // this function performs depth-first search 
    // starting at node i. break any tie based on 
    // node number (smaller node goes first) e.g., 
    // if you can pick node 2 or node 3, pick 2. 
    // at last, this function should return an array 
    // of size n holding the traverse sequence of nodes. 
    // (Assume input graph is always connected.)
    int* DFT(int i);


    // this function performs breadth-first search 
    // starting at node i. when exploring neighbors 
    // of a set of nodes, explore them based on the 
    // order of nodes in the queue. 
    // 
    // this means once you pop a node from the queue, 
    // add its neighbors to the queue. (here, break 
    // ties based on neighbor node numbers -- smaller 
    // node gets added to the queue first)
    // 
    // at last, this function should return an array 
    // of size n holding the traverse sequence of nodes. 
    int* BFT(int i);

    // 
    // The following performs the Dijkstra's algorithm
    // to find the shorest path from node i to node j.  
    // 
    // It returns address of an object of the 
    // Result class, which contains three 
    // public variables (see definition at top): 
    // (i) int length: length of the shorest path 
    // (ii) int weight: total weight of the shortest path
    // (iii) int *path: array of nodes on the path 
    // Example: 
    // If the shortest path is 2 -> 3 -> 0, and 
    // weight on (2,3) is 5 and weight on (3,0) is 1, 
    // then path[0] = 2, path[1] = 3, path[2] = 0
    // and length = 3 and weight = 6. 
    // 
    Result* Dijkstra(int i, int j);


    // this is the constructor with input arguments 
    // remember to allocate space for matrix "m" and 
    // initialize all entries to zero 
    // 
    // Initialization is important because in "main" 
    // we only add edges to the matrix and assume its 
    // rest entries are zeros. 
    Graph(int n);
};

// This function returns the degree of node i 
int Graph::Degree(int i){
    
    int degree = 0;
    int k,j;

    for(j = 0; j < size; j++){
        if(m[i][j] != 0) degree++; 
    }

    return degree;
}

// Add
void Graph::Add(int i, int j, int w){
    m[i][j] = w;
    m[j][i] = w;
}

// Is edge
int Graph::IsEdge(int i, int j){
    if(m[i][j] != 0) return 1;
    else return 0;
}

// Is path  
int Graph::IsPath(int i, int j){
    
    //Use BFT to determine path
    int* BFTSequence = BFT(i);

    for(int k = 0; k < size; k++){
        if(BFTSequence[k] == j) return 1;
    }

    return 0;
}

// DFT
int* Graph::DFT(int i){

    // Create an array to store the traversal sequence
    int* traversalSequence = new int[size];

    // Create a boolean array to keep track of visited vertices
    bool* visited = new bool[size];

    // Initialize all vertices as not visited
    for(int k = 0; k < size; k++) {
        visited[k] = false;
    }

    // Create a stack for iterative traversal
    std::stack<int> stack;

    // Push the i onto the stack
    stack.push(i);

    // Counter to keep track of the number of visited vertices
    int count = 0;

    // Iterative Depth-First Traversal using a stack
    while(!stack.empty()) {
        
        // Get the current vertex from the stack
        int current = stack.top();
        stack.pop();

        // Process the current vertex if not visited
        if(!visited[current]) {
            // Mark the current vertex as visted
            visited[current] = true;
            // Store the current vertex in the traversal sequence
            traversalSequence[count++] = current;

            // Push adjacent vertecies onto the stack from largest to lowest
            for(int j = size - 1; j >= 0; j--){
                if(m[current][j] != 0 && !visited[j]){
                    stack.push(j);
                }
            }
        }
    }

    delete[] visited;

    return traversalSequence;
}


//BFT
int* Graph::BFT(int i){
    
    // Create an array to store the traversal sequence
    int* traversalSequence = new int[size];

    // Create a boolean array to keep track of visited vertices
    bool* visited = new bool[size];


    // Initialize all vertices as not visited
    for(int k = 0; k < size; k++) {
        visited[k] = false;
    }

    // Make queue
    std::queue<int> q;  
    q.push(i);
    // Flag starting node i as visited 
    visited[i] = true;

    // Counter to keep track of the number of visited vertices
    int count = 0;

    while(!q.empty()){
        // Get the front of the queue
        int current = q.front(); 
        // Pop
        q.pop();
        //Store current int traversalSequence
        traversalSequence[count++] = current;

        //Explore neighbors of the current node
        for(int j = 0; j < size; j++){
            //Check if there is an edge and the neighbor is not visited
            if(m[current][j] != 0 && !visited[j]){
                // Add the neighbor to the queue
                q.push(j);
                // Mark the neighbor as visited
                visited[j] = true; 
            }
        }
    }

    delete[] visited;

    return traversalSequence;
}

//Dijkstra
Result* Graph::Dijkstra(int i, int j){

    // Array to store the shortest distances
    int* dist = new int[size];

    // Array to store the shortest path
    int* path = new int[size];

    // Array to keep track of visited nodes
    bool* visited = new bool[size];

    for(int k = 0; k < size; k++){
        dist[k] = INT_MAX; //Initialize distances to infinity
        path[k] = -1; // Initialize path to -1 (indicating no path)
        visited[k] = false; // Initialize all nodes as not visited
    }

    // Distance from source to self is 0
    dist[i] = 0;

    //Queue to store vertices
    queue<int> q;

    q.push(i);

    while(!q.empty()){
        // Find the vertex with the minimum distance in the current queue
        int u = q.front();

        for(int k = 0; k < size; k++){
            if(!visited[k] && dist[k] < dist[u]){
                u = i;
            }
        }

        // Remove the vertex with the min distance from the queue
        q.pop();

        // Mark the current node as visited
        visited[u] = true;

        //Update the distance and path for each neighbor of u
        for (int v = 0; v < size; ++v) {
            if (!visited[v] && m[u][v] && dist[u] != INT_MAX && dist[u] + m[u][v] < dist[v]) {
                dist[v] = dist[u] + m[u][v];
                path[v] = u;
                q.push(v);
            }
        }
    }

    // Build the path array
    int pathLength = 0;
    int current = j;
    while (current != -1) {
        current = path[current];
        ++pathLength;
    }

    int* shortestPath = new int[pathLength];
    current = j;
    for (int i = pathLength - 1; i >= 0; --i) {
        shortestPath[i] = current;
        current = path[current];
    }

    Result* result = new Result(pathLength, dist[j], shortestPath);

    delete[] dist;
    delete[] path;
    delete[] visited;

    return result;
}

//Constructor
Graph::Graph(int n){
    size = n;

    //allocate space for array of int pointers
    m = new int*[size];

    //for each pointer allocate an array
    for(int i = 0; i < size; i++){
        m[i] = new int[size];
    }

    //initialize all values to 0
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size ; j++){
            m[i][j] = 0;
        }
    }
}


int main()
{

    int mode, size, i, j, w;

    int a, b; // node numbers used for testing 

    cin >> mode >> size >> a >> b;

    Graph x(size);

    // each time input a pair 
    // of indices and a weight 
    // remember to separate them 
    // when inputing from keyboard 
    // also, we assume inputs are 
    // within boundary 
    while (cin >> i >> j >> w) {
        x.Add(i, j, w);
    }

    // Mode 0: test IsEdge()
    if (mode == 0) {
        cout << x.IsEdge(a, b);
    }
    // Mode 1: test IsPath()
    else if (mode == 1) {
        cout << x.IsPath(a, b);
    }
    // Mode 2: test Degree()
    else if (mode == 2) {
        cout << x.Degree(a);
    }
    // Mode 3: test DFT()
    else if (mode == 3) {
        int* s = new int[size];
        s = x.DFT(a);
        for (int i = 0; i < size; i++) {
            cout << s[i];
        }
    }
    // Mode 4: test BFT()
    else if (mode == 4) {
        int* s = new int[size];
        s = x.BFT(a);
        for (int i = 0; i < size; i++) {
            cout << s[i];
        }
    }
    // Mode 5: test Dijkstra()
    else if (mode == 5) {
        Result *z = x.Dijkstra(a, b);
        cout << z->length << endl;
        cout << z->weight << endl;
        for (int i = 0; i < z->length; i++) {
            cout << z->path[i];
        }
    }

    return 0;

}
