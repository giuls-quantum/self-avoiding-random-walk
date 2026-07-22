#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>

// ---------------------------------------------------------
// Fast, high-quality PRNG (Xorshift32)
// ---------------------------------------------------------
uint32_t prng_state = 123456789;

uint32_t xorshift32() {
    prng_state ^= prng_state << 13;
    prng_state ^= prng_state >> 17;
    prng_state ^= prng_state << 5;
    return prng_state;
}

// Generate a random integer between 0 and max-1
int rand_int(int max) {
    return xorshift32() % max;
}

// ---------------------------------------------------------
// Main SARW Simulation
// ---------------------------------------------------------
int main() {
    // Configuration
    const int MAX_STEPS = 50;
    const int NUM_TRIALS = 100000;
    
    // Seed the custom PRNG
    prng_state = (uint32_t)time(NULL);

    // Grid size for O(1) collision lookups: (2N + 1)
    const int GRID_SIZE = 2 * MAX_STEPS + 1;
    const int ORIGIN = MAX_STEPS;
    
    // Allocate grid on the heap to prevent stack overflow for large N
    bool** visited = (bool**)malloc(GRID_SIZE * sizeof(bool*));
    for (int i = 0; i < GRID_SIZE; i++) {
        visited[i] = (bool*)calloc(GRID_SIZE, sizeof(bool));
    }

    // Directions: Up, Right, Down, Left
    int dx[] = {0, 1, 0, -1};
    int dy[] = {1, 0, -1, 0};

    int successful_walks = 0;
    double total_squared_distance = 0.0;

    printf("Starting %d SARW trials with max steps = %d\n", NUM_TRIALS, MAX_STEPS);

    for (int trial = 0; trial < NUM_TRIALS; trial++) {
        int x = ORIGIN;
        int y = ORIGIN;
        visited[x][y] = true;
        
        bool trapped = false;
        
        for (int step = 0; step < MAX_STEPS; step++) {
            // Check if completely trapped before picking
            if (visited[x][y+1] && visited[x+1][y] && 
                visited[x][y-1] && visited[x-1][y]) {
                trapped = true;
                break;
            }
            
            // Pick a valid random direction
            int dir;
            int nx, ny;
            do {
                dir = rand_int(4);
                nx = x + dx[dir];
                ny = y + dy[dir];
            } while (visited[nx][ny]); // O(1) collision check!
            
            x = nx;
            y = ny;
            visited[x][y] = true;
        }
        
        // If it survived, record metrics
        if (!trapped) {
            successful_walks++;
            int dist_sq = (x - ORIGIN)*(x - ORIGIN) + (y - ORIGIN)*(y - ORIGIN);
            total_squared_distance += dist_sq;
        }
        
        // Reset grid efficiently for the next trial
        for (int i = 0; i < GRID_SIZE; i++) {
            for(int j = 0; j < GRID_SIZE; j++){
                visited[i][j] = false;
            }
        }
    }

    printf("Successful walks: %d / %d\n", successful_walks, NUM_TRIALS);
    if (successful_walks > 0) {
        printf("Mean Squared End-to-End Distance <R^2>: %.3f\n", 
               total_squared_distance / successful_walks);
    }

    // Cleanup
    for (int i = 0; i < GRID_SIZE; i++) free(visited[i]);
    free(visited);

    return 0;
}
