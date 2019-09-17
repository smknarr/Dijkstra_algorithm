#include <stdio.h> 
#include <stdlib.h> 
#include <limits.h> 
#include <math.h>
// A structure to represent a node in adjacency list 
struct ListNode 
{ 
	int dest; 
	struct ListNode* next; 
}; 
// A structure to represent an adjacency liat 
struct List 
{ 
	struct ListNode *head; // pointer to head node of list 
}; 
// A structure to represent a graph. A graph is an array of adjacency lists. 
// Size of array will be V (number of vertices in graph) 
struct Graph 
{ 
	int V; 
	struct List* array; 
}; 
typedef struct Graph* Graph;
// A utility function to create a new adjacency list node 
struct ListNode* newListNode(int dest) 
{ 
	struct ListNode* newNode = 
			(struct ListNode*) malloc(sizeof(struct ListNode)); 
	newNode->dest = dest;
	newNode->next = NULL; 
	return newNode; 
} 




// A utility function that creates a graph of V vertices 
Graph createGraph(int V) 
{ 
	int i;
	Graph graph = (Graph) malloc(sizeof(struct Graph)); 
	graph->V = V; 
	// Create an array of adjacency lists. Size of array will be V 
	graph->array = (struct List*) malloc(V * sizeof(struct List)); 
	// Initialize each adjacency list as empty by making head as NULL 
	for (i = 0; i < V; ++i) 
		graph->array[i].head = NULL; 
	return graph; 
} 
// Adds an edge to an undirected graph 
void addEdge(Graph graph, int src, int dest) 
{ 
	// Add an edge from src to dest. A new node is added to the adjacency list of src. The node is added at the begining 
	struct ListNode* newNode = (struct ListNode*) malloc(sizeof(struct ListNode)); 
	newNode->dest = dest;
	newNode->next = NULL;
	newNode->next = graph->array[src].head; 
	graph->array[src].head = newNode; 
	// Since graph is undirected, add an edge from dest to src also 
	newNode = (struct ListNode*) malloc(sizeof(struct ListNode)); 
	newNode->dest = src;
	newNode->next = NULL; 
	newNode->next = graph->array[dest].head; 
	graph->array[dest].head = newNode; 
} 
// Structure to represent a min heap node 
typedef struct  
{ 
	int v; 
	int dist; 
}MinHeapNode; 
// Structure to represent a min heap 
struct MinHeap 
{ 
	int size;	 // Number of heap nodes present currently 
	int capacity; // Capacity of min heap 
	int *pos;	 // This is needed for decreaseKey() 
	MinHeapNode **array; 
}; 



// A utility function to create a new Min Heap Node 
MinHeapNode* newMinHeapNode(int v, int dist) 
{ 
	MinHeapNode* minHeapNode = (MinHeapNode*)malloc(sizeof(MinHeapNode)); 
	minHeapNode->v = v; 
	minHeapNode->dist = dist; 
	return minHeapNode; 
} 
// A utility function to create a Min Heap 
struct MinHeap* createMinHeap(int capacity) 
{ 
	struct MinHeap* minHeap = (struct MinHeap*)malloc(sizeof(struct MinHeap)); 
	minHeap->pos = (int *)malloc(capacity * sizeof(int)); 
	minHeap->size = 0; 
	minHeap->capacity = capacity; 
	minHeap->array = (MinHeapNode*)malloc(capacity * sizeof(MinHeapNode)); 
	return minHeap; 
} 
// A utility function to swap two nodes of min heap. Needed for min heapify 
void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b) 
{ 
	MinHeapNode* t = *a; 
	*a = *b; 
	*b = t; 
} 
// A standard function to heapify at given idx. This function also updates position of nodes when they are swapped. 
// Position is needed for decreaseKey() 
void minHeapify(struct MinHeap* minHeap, int idx) 
{ 
	int smallest, left, right; 
	smallest = idx; 
	left = 2 * idx + 1; 
	right = 2 * idx + 2; 

	if (left < minHeap->size && minHeap->array[left]->dist < minHeap->array[smallest]->dist ) 
		smallest = left; 

	if (right < minHeap->size && minHeap->array[right]->dist <= minHeap->array[smallest]->dist ) 
		smallest = right; 

	if (smallest != idx) 
	{ 
		// The nodes to be swapped in min heap 
		MinHeapNode *smallestNode = minHeap->array[smallest]; 
		MinHeapNode *idxNode = minHeap->array[idx]; 
		
                        // Swap positions 
		minHeap->pos[smallestNode->v] = idx; 
		minHeap->pos[idxNode->v] = smallest; 

		// Swap nodes 
		swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]); 

		minHeapify(minHeap, smallest); 
	} 
} 
// A utility function to check if the given minHeap is ampty or not 
int isEmpty(struct MinHeap* minHeap) 
{ 
	return minHeap->size == 0; 
} 
// Standard function to extract minimum node from heap 
MinHeapNode* extractMin(struct MinHeap* minHeap) 
{ 
	if (isEmpty(minHeap)) 
		return NULL; 
	// Store the root node 
	MinHeapNode* root = minHeap->array[0]; 
	// Replace root node with last node 
	MinHeapNode* lastNode = minHeap->array[minHeap->size - 1]; 
	minHeap->array[0] = lastNode; 
	// Update position of last node 
	minHeap->pos[root->v] = minHeap->size-1; 
	minHeap->pos[lastNode->v] = 0; 
	// Reduce heap size and heapify root 
	--minHeap->size; 
	minHeapify(minHeap, 0); 
	return root; 
} 
// Function to decreasy dist value of a given vertex v. This function uses pos[] of min heap to get the current index of node in min heap 
void decreaseKey(struct MinHeap* minHeap, int v, int dist) 
{ 
	// Get the index of v in heap array 
	int i = minHeap->pos[v]; 
	// Get the node and update its dist value 
	minHeap->array[i]->dist = dist; 
	// Travel up while the complete tree is not hepified. This is a O(Logn) loop 
	while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) 
	{ 
		// Swap this node with its parent 
		minHeap->pos[minHeap->array[i]->v] = (i-1)/2; 
		minHeap->pos[minHeap->array[(i-1)/2]->v] = i; 
		swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]); 
		// move to parent index 
		i = (i - 1) / 2; 
	} 
} 
// A utility function to check if a given vertex 'v' is in min heap or not 
int isInMinHeap(struct MinHeap *minHeap, int v) 
{ 
if (minHeap->pos[v] < minHeap->size) 
	return 1; 
return 0; 
} 
//A function to give directions to the players. At every node the decision to be made is displayed.
void direction(int co[][2],int n)
{
	int j,o=1; // o is the orientation of the player
	for(j=n+1;j>=0;j--)
		{
		co[j+1][0]=co[j][0];	co[j+1][1]=co[j][1];
		}
	//n++;
	co[0][0]=co[0][1]=0;
	for(j=0;j<n;j++)
			{
					//R denotes right, L denotes Left, S denotes straight
					switch(o)
				{
					case 1: if(co[j+1][0]-co[j][0]>0)
								{printf(" R"); o=4;}
							else if(co[j+1][0]-co[j][0]<0)
								{printf(" L"); o=2;}
							else printf(" S");
							break;
					case 2: if(co[j+1][1]-co[j][1]>0)
								{printf(" R"); o=1;}
							else if(co[j+1][1]-co[j][1]<0)
									{printf(" L"); o=3;}
							else printf(" S");
							break;
					case 3: if(co[j+1][0]-co[j][0]<0)
							{printf(" R"); o=2;}
							else if(co[j+1][0]-co[j][0]>0)
									{printf(" L"); o=4;}
								else printf(" S");
							break;
					

                                                             case 4: if(co[j+1][1]-co[j][1]<0)
								{printf(" R"); o=3;}
							else if(co[j+1][1]-co[j][1]>0)
									{printf(" L"); o=1;}
								else printf(" S");
							break;
				}
			 } 
	
}
//To obtain the path from the shortest path tree by the method of backtracking
void getPath(int parent[], int j, int n, int co[][2],int index) 
{ 
		
	// Base Case : If j is source 
	
	if (parent[j] == - 1) 
		return; 
	index--;
	getPath(parent, parent[j],n,co,index); 
	co[index][0]=j/n;	co[index][1]=j%n; 
} 
// A utility function used to print the solution 
void printArr(int dist[], int n, int parent[]) 
{ 
	int i,j;
	printf("Vertex Distance from Source Path\n"); 
	for (i = 1; i < n; i++) 
	{
		
		if(dist[i]==INT_MAX)
			printf("\n(%d,%d) \t\t No Path",i/(int)sqrt(n),i%(int)sqrt(n));
		else
			{
			int co[dist[i]][2];
			for(j=0;j<dist[i];j++)
				co[j][0]=co[j][1]=0;	
			printf("\n(%d,%d) \t\t %d\t\t(0,0)", i/(int)sqrt(n),i%(int)sqrt(n), dist[i]); 
			getPath(parent, i, sqrt(n), co,dist[i]);
			for(j=0;j<dist[i];j++)
				printf("->(%d,%d)",co[j][0],co[j][1]);
			direction(co,dist[i]);
			}
 	}
}
// The main function that calulates distances of shortest paths from src to all vertices. It is a O(ELogV) function 
void dijkstra(Graph graph, int src) 
{ 
	int co[2][2];
	int V = graph->V,v;// Get the number of vertices in graph 
	int dist[V];	 // dist values used to pick minimum weight edge in cut 
	int parent[V];	// Parent array to store shortest path tree. parent[i] stores the parent of node i in the shortest path tree
	parent[0]=-1;
	// minHeap represents set E 
	struct MinHeap* minHeap = createMinHeap(V); 
	// Initialize min heap with all vertices. dist value of all vertices 
	dist[src] = 0;	// Make dist value of src vertex as 0 so that it is extracted first 
	for (v = 0; v < V; ++v) 
	{ 
		if(v!=src)
		{
			dist[v] = INT_MAX; 
		}
			minHeap->array[v] = newMinHeapNode(v, dist[v]); 
			minHeap->pos[v] = v; 
		
	} 
	
	// Initially size of min heap is equal to V 
	minHeap->size = V; 
	// In the following loop, min heap contains all nodes whose shortest distance is not yet finalized. 
	
	while (minHeap->size != 0)//|| u != des) 
	{ 
		// Extract the vertex with minimum distance value 
		MinHeapNode* minHeapNode = extractMin(minHeap); 
		int u = minHeapNode->v; // Store the extracted vertex number 
		 //Traverse through all adjacent vertices of u (the extracted vertex) and update their distance values 
		struct ListNode* pCrawl = graph->array[u].head; 
		while (pCrawl != NULL) 
		{ 
			int v = pCrawl->dest; 
			// If shortest distance to v is not finalized yet, and distance to v through u is less than its previously calculated distance 
			if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX &&(1 + dist[u] < dist[v])) 
			{ 
				dist[v] = dist[u] + 1;	// update distance value in min heap also 
				parent[v]=u;
				decreaseKey(minHeap, v, dist[v]); 
			} 
			
                         pCrawl = pCrawl->next; 
		} 
	} 
	// print the calculated shortest distances 
	//for(v=0;v<V;v++)
	//	printf("%d ",parent[v]);
	
	printf("\n\n");
	//getPath(parent,des,V,co);
	printArr(dist, V,parent); 
} 
int main() 
{ 
	int V;	int i,j; 
	printf("\nEnter size of Grid: ");	scanf("%d",&V);
	V++;
	Graph graph = createGraph(V*V); 
	int a[V][V];
	printf("\nEnter the matrix layout of the grid: (0 if node is open, 1 if node is closed)\n");
	//Input the matrix
	for(i=0;i<V;i++)
		for(j=0;j<V;j++)
			scanf("%d",&a[i][j]);
	//Add edges to the graph based on the grid matrix		
	for(i=0;i<V;i++)
		for(j=0;j<V;j++)
		{
			if(a[i][j])
				{
					//Check whether current node is connected to right node
					if(i+1<V&&a[i+1][j])
						addEdge(graph,V*i+j,V*(i+1)+j);
					//Check whether current node is connected to right node
					if(j+1<V&&a[i][j+1])
						addEdge(graph,V*i+j,V*i+j+1);
				}
		}
	dijkstra(graph, 0); 
	return 0; 
}