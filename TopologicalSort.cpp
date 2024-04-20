#include <omp.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_NODES 1000
#define MAX_EDGES 500*999

// Structure to represent a stack
struct Stack {
    int data[MAX_NODES];
    int top;
};

// Function to initialize the stack
void initialize(struct Stack *s) {
    s->top = -1;
}

// Function to push an element onto the stack
void push(struct Stack *s, int value) {
    if (s->top == MAX_NODES - 1) {
        // printf("Stack overflow\n");
        return;
    }
    s->data[++s->top] = value;
}

// Function to pop an element from the stack
int pop(struct Stack *s) {
    if (s->top == -1) {
        // printf("Stack underflow\n");
        exit(1);
    }
    return s->data[s->top--];
}

// Function to check if the stack is empty
int isEmpty(struct Stack *s) {
    return (s->top == -1);
}

// Function to perform DFS and topological sorting (Sequential Version)
void topologicalSortUtilSequential(int v, int adj[MAX_NODES][MAX_NODES], int visited[MAX_NODES], struct Stack *Stack, int V) {
    visited[v] = 1; // Mark the current node as visited

    // Recur for all adjacent vertices
    for (int i = 0; i < V; i++) {
        if (adj[v][i] && !visited[i])
            topologicalSortUtilSequential(i, adj, visited, Stack, V);
    }

    // Push current vertex to stack which stores the result
    push(Stack, v);
}

// Function to perform Topological Sort (Sequential Version)
void topologicalSortSequential(int adj[MAX_NODES][MAX_NODES], int V) {
    struct Stack Stack;
    initialize(&Stack);
    int visited[MAX_NODES];
    for (int i = 0; i < V; i++)
        visited[i] = 0; // Mark all vertices as not visited

    // Call the recursive helper function to store
    // Topological Sort starting from all vertices one by
    // one
    for (int i = 0; i < V; i++) {
        if (!visited[i])
            topologicalSortUtilSequential(i, adj, visited, &Stack, V);
    }

    // Print contents of stack
    // printf("Topological sorting of the graph (Sequential): ");
    // while (!isEmpty(&Stack)) {
    //     printf("%d ", pop(&Stack));
    // }
}

// Function to perform DFS and topological sorting (Parallel Version)
void topologicalSortUtilParallel(int v, int adj[MAX_NODES][MAX_NODES], int visited[MAX_NODES], struct Stack *Stack, int V) {
    visited[v] = 1; // Mark the current node as visited

    // Recur for all adjacent vertices
    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < V; i++) {
        if (adj[v][i] && !visited[i])
            topologicalSortUtilParallel(i, adj, visited, Stack, V);
    }

    // Push current vertex to stack which stores the result
    #pragma omp critical
    {
        push(Stack, v);
    }
}

// Function to perform Topological Sort (Parallel Version)
void topologicalSortParallel(int adj[MAX_NODES][MAX_NODES], int V, int num_threads) {
    struct Stack Stack;
    initialize(&Stack);
    int visited[MAX_NODES];
    for (int i = 0; i < V; i++)
        visited[i] = 0; // Mark all vertices as not visited

    // Set number of threads
    omp_set_num_threads(num_threads);

    // Call the recursive helper function to store
    // Topological Sort starting from all vertices one by
    // one
    #pragma omp parallel for
    for (int i = 0; i < V; i++) {
        if (!visited[i])
            topologicalSortUtilParallel(i, adj, visited, &Stack, V);
    }

    // // Print contents of stack
    // printf("\nTopological sorting of the graph (Parallel): ");
    // while (!isEmpty(&Stack)) {
    //     printf("%d ", pop(&Stack));
    // }
}

int main() {
    int nodes[] = {100, 200, 300, 400 ,500 ,600, 700}; // Array of node counts to iterate over
    int num_tests = sizeof(nodes) / sizeof(nodes[0]); // Number of tests to run

    for (int test = 0; test < num_tests; test++) {
        int V = nodes[test]; // Current number of nodes
        int E = V * (V - 1) / 2; // Max edges for a complete graph with V nodes

        // Ensure E does not exceed MAX_EDGES
        if (E > MAX_EDGES) {
            printf("Error: Number of edges (%d) for V = %d exceeds MAX_EDGES.\n", E, V);
            continue;
        }

        // Edges
        int edges[MAX_EDGES][2];
        E = 0; // Reset E to use as counter
        for (int i = 0; i < V - 1; i++) {
            for (int j = i + 1; j < V; j++) {
                edges[E][0] = i;
                edges[E][1] = j;
                E++;
            }
        }

        // Graph represented as an adjacency matrix
        int adj[MAX_NODES][MAX_NODES] = {0};
        for (int i = 0; i < E; i++) {
            adj[edges[i][0]][edges[i][1]] = 1;
            // Assuming undirected graph for simplicity, remove if directed
            adj[edges[i][1]][edges[i][0]] = 1; // Add this line if the graph is undirected
        }

        // Number of threads
        int num_threads = 8;

        printf("\nRunning test for V = %d, E = %d\n", V, E);

        // Sequential Version
        double start_seq = omp_get_wtime();
        topologicalSortSequential(adj, V);
        double end_seq = omp_get_wtime();
        printf("\nSequential Execution Time: %lf seconds", end_seq - start_seq);

        // Parallel Version
        double start_parallel = omp_get_wtime();
        topologicalSortParallel(adj, V, num_threads);
        double end_parallel = omp_get_wtime();
        printf("\nParallel Execution Time: %lf seconds\n", end_parallel - start_parallel);

        // Calculate Speedup
        double speedup = (end_seq - start_seq) / (end_parallel - start_parallel);
        printf("Speedup: %lf\n", speedup);

        // Calculate Efficiency
        double efficiency = speedup / num_threads;
        printf("Efficiency: %lf\n", efficiency);
    }

    return 0;
}
