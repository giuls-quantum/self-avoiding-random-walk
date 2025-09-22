// selfa_random_walk.c
//  Created by Giulia Maniccia on 22/09/25

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#define L 100           // grid size
#define x0 50           // starting x
#define y0 50           // starting y
#define Ttot 500        // maximum time per walk
#define TMIN 20         // minimum steps for averages
#define NPROVE 5000     // number of simulations
#define SEED 73756933   // fixed seed. Comment this line and uncomment the next one to switch to a time-based seed
// #define SEED ((unsigned int)time(NULL))

#define RIGHT 0
#define LEFT 1
#define UP 2
#define DOWN 3

typedef unsigned long long int RANDOM;
/************************************** function prototypes **********************************/

int* arraydin(int dim);
double* doublearraydin(int dim);
void intro();
int randomgen(int max);
int verify_neighbor(int *grid, int x, int y, int *xP, int *xM, int *yP, int *yM, bool periodic);
int findspot(int *grid, int x, int y, char **walk, int i, int j, int *xP, int *xM, int *yP, int *yM, bool periodic);
char** walkmatrix();
char* oppositewalk(char *walkexact);
int samewalk(char **walk, char *walkexact, int index);
int countwalk1(char **walk);
int countwalk2(char **walk);
void arrayborderperiodic(int *xP, int *xM, int *yP, int *yM);
void printwalk(const char *datafile, const char *outfile);
void printdistance(const char* datafile, const char* outfile);
void printaverage(const char* datafile, const char* outfile);
void printprob(const char* datafile, const char* outfile);

/************************************* main*************************************/

int main() {
    intro();
    srand48(SEED);
    // setting up arrays
    double *DeltaR = doublearraydin(Ttot+1);
    int *grid = arraydin(L*L); //save grid as a 1D array
    double *lengthN = doublearraydin(Ttot);
    int *xP = arraydin(L);
    int *xM = arraydin(L);
    int *yP = arraydin(L);
    int *yM = arraydin(L);
    arrayborderperiodic(xP,xM,yP,yM);
    char **walk = walkmatrix();
    // setting up variables
    int t=0, i, j, k, spot, freespot, x, y, deltar, counter=0;
    double averageR=0, averageT=0;
    // setting up data files
    FILE* positions = fopen("data/SARW_positions.dat","w");
    fprintf(positions,"#Time\t x\t y\t Distance from origin\n%d\t%d\t%d\t%d\n",0,x0,y0,0);
    FILE* averages = fopen("data/SARW_averages.dat","w");
    fprintf(averages,"#Time\t Distance from origin\n");
    FILE* distance = fopen("data/SARW_distance.dat","w");
    fprintf(distance,"#Time\t DeltaR[t]\n");
    FILE* probability = fopen("data/SARW_probability.dat","w");
    fprintf(probability,"#Walk length\t Probability\n");

    printf("Maximum steps: %d, number of attempts: %d\n", Ttot, NPROVE);
  
    // SIMULATION START
    for(j=0;j<NPROVE;j++) {
        x = x0;
        y = y0;
        grid[y0+L*x0] = 1; //spots crossed are indentified with label 1
        t = 0;
        deltar = 0;

        // WALK START
        for(i=0;i<Ttot;i++) {
            t++;
            // Set true for periodic border conditions, false for standard border conditions
            freespot = verify_neighbor(grid, x, y, xP, xM, yP, yM, true); // verify close cells: free empty spot=0, no available spot=-1
                
            // no free spot available:
            if(freespot < 0) {
                deltar = ((x-x0)*(x-x0) + (y-y0)*(y-y0));
                printf("#*********** Dead end reached **********\nAttempt n. %d:\tTime %d\t Distance from origin %d\n",j,t,deltar);
                if(j==0) fprintf(positions,"%d\t %d\t %d\t %d\n",t,x,y,deltar); //saving only the first walk steps
                if(t<TMIN) counter++;
                fprintf(averages,"%d\t %d\t %d\n",j,t,deltar);
                DeltaR[t] += (double)deltar/(double)NPROVE;
                averageR += deltar;
                averageT += t;
                lengthN[t] += 1.;   //save how many walks of length=t
                break; //ends cycle
                    
            // there is at least one spot availabla:
            } else {
                spot = findspot(grid, x, y, walk, i, j, xP, xM, yP, yM, true); //move randomly in a free spot
                x = spot / L;
                y = spot % L;
                deltar = ((x-x0)*(x-x0) + (y-y0)*(y-y0));
                grid[spot] = 1;
                DeltaR[t] += (double)deltar/(double)NPROVE;
                if(j==0) fprintf(positions,"%d\t %d\t %d\t %d\n",t,x,y,deltar);
            }
        }
        // WALK ENDS

    for(k=0;k<L*L;k++) grid[k] = 0;
    }
    // SIMULATION ENDS

    averageT /= NPROVE;
    averageR /= NPROVE;
    printf("Average time to reach a dead-end: %lf, average distance: %lf\n", averageT, averageR);

    // Estimate probability to reach a certain distance
    for(i=0;i<Ttot;i++) {
        lengthN[i] /= NPROVE;
        fprintf(probability,"%d\t%lf\n",i,lengthN[i]);
        if((i<TMIN)&&(DeltaR[i]>0))fprintf(distance,"%d\t %lf\n",i,DeltaR[i]); // using low t values to avoid discontinuities of dead-end walks
    }

    fclose(positions);
    fclose(averages);
    fclose(distance);
    fclose(probability);
    // Plots
    printwalk("data/SARW_positions.dat", "plots/SelfA_Random_Walk.pdf");
    printaverage("data/SARW_averages.dat", "plots/SARW_averageR.pdf");
    printdistance("data/SARW_distance.dat", "plots/SARW_distance.pdf");
    printprob("data/SARW_probability.dat", "plots/SARW_prob.pdf");
   
    // CLEANUP
    free(grid);
    free(DeltaR);
    free(lengthN);
    for(j=0;j<NPROVE;j++) free(walk[j]);
    free(walk);
    free(xP); free(xM); free(yP); free(yM);
    fclose(positions); fclose(averages); fclose(distance); fclose(probability);
    return 0;
}
 
/***************************************** program end ***********************************/


/***************************************** function definitons ****************************************/

void intro() {
    printf("Self-Avoiding Random Walk on a %dx%d grid\n", L, L);
}
/*******************************************************************************************************/
int* arraydin(int dim) {
    int* arr = (int*)calloc(dim, sizeof(int));
    if (!arr) { perror("calloc"); exit(EXIT_FAILURE); }
    return arr;
}
/*******************************************************************************************************/
double* doublearraydin(int dim) {
    double* arr = (double*)calloc(dim, sizeof(double));
    if (!arr) { perror("calloc"); exit(EXIT_FAILURE); }
    return arr;
}
/*******************************************************************************************************/
int randomgen(int max) {
    return (int)(drand48() * max);
}
/*******************************************************************************************************/
int verify_neighbor(int *grid, int x, int y, int *xP, int *xM, int *yP, int *yM, bool periodic) {
    int right, left, up, down;
    if(periodic) {
        //periodic border conditions = pacman effect
        right = (grid[y + L*xP[x]]==0);
        left  = (grid[y + L*xM[x]]==0);
        up    = (grid[yP[y]+L*x]==0);
        down  = (grid[yM[y]+L*x]==0);
    } else {
        //no periodic border conditions == wall at each border
        right = (x<L-1 && grid[y + L*(x+1)]==0);
        left  = (x>0  && grid[y + L*(x-1)]==0);
        up    = (y<L-1 && grid[y+1+L*x]==0);
        down  = (y>0  && grid[y-1+L*x]==0);
    }
    return (right || left || up || down) ? 0 : -1;
}
/*******************************************************************************************************/
int findspot(int *grid, int x, int y, char **walk, int i, int j, int *xP, int *xM, int *yP, int *yM, bool periodic) {
    int spot, dir;
    do {
        //pick a random direction and check neighboring cells
        dir = randomgen(4);
        if(dir==RIGHT) spot = periodic ? y + L*xP[x] : y + L*(x+1);
        else if(dir==LEFT) spot = periodic ? y + L*xM[x] : y + L*(x-1);
        else if(dir==UP) spot = periodic ? yP[y] + L*x : y+1+L*x;
        else spot = periodic ? yM[y] + L*x : y-1+L*x;
    } while(grid[spot]>0);
    // save in which direction to move
    if(dir==RIGHT)walk[j][i]='R';
    if(dir==LEFT)walk[j][i]='L';
    if(dir==UP)walk[j][i]='U';
    if(dir==DOWN)walk[j][i]='D';
    // return available spot
    return spot;
}
/*******************************************************************************************************/
char** walkmatrix() {
    char **walk = (char**)calloc(NPROVE, sizeof(char*));
    if(!walk) { perror("calloc"); exit(EXIT_FAILURE); }
    for(int p=0;p<NPROVE;p++){
        walk[p] = (char*)calloc(Ttot, sizeof(char));
        if(!walk[p]) { perror("calloc"); exit(EXIT_FAILURE); }
    }
    return walk;
}
/*******************************************************************************************************/
void arrayborderperiodic(int *xP, int *xM, int *yP, int *yM) {
    for(int a=0;a<L;a++){
        xP[a]=(a+1)%L; xM[a]=(a-1+L)%L;
        yP[a]=(a+1)%L; yM[a]=(a-1+L)%L;
    }
}
/*******************************************************************************************************/
void printwalk(const char *datafile, const char *outfile){
    FILE *plot = popen("gnuplot -persist", "w");
    if (!plot) { perror("popen"); return; }
    fprintf(plot, "set terminal pdfcairo\n");
    fprintf(plot, "set output '%s'\n",outfile);
    fprintf(plot, "set xlabel 'x'\n");
    fprintf(plot, "set ylabel 'y'\n");
    fprintf(plot, "set title 'Self-Avoiding Random Walk (Top View)'\n");
    fprintf(plot, "set palette defined (0 '#0000FF', 0.5 '#00FF00', 1 '#FF0000')\n");
    fprintf(plot, "set view 0\n set grid back\n set format z ''\n unset border\n");
    // plots walk: col2 = x, col3 = y, col1 = time
    fprintf(plot, "splot '%s' using 2:3:1 with linespoints pointtype 7 pointsize 0.3 palette title 'Evolution in time'\n", datafile);
    fprintf(plot, "set output\n");
    pclose(plot);
}
/*******************************************************************************************************/
void printaverage(const char *datafile, const char *outfile){
    FILE *plot = popen("gnuplot -persist", "w");
    if (!plot) { perror("popen"); return; }
    fprintf(plot, "set terminal pdfcairo\n");
    fprintf(plot, "set output '%s'\n",outfile);
    fprintf(plot, "set xlabel 'Attempt'\n");
    fprintf(plot, "set ylabel 'Average distance from origin'\n");
    fprintf(plot, "plot '%s' using 1:3 title 'Squared distance'\n", datafile);
    fprintf(plot, "set output\n");
    pclose(plot);
}
/*******************************************************************************************************/
void printdistance(const char *datafile, const char *outfile){
    FILE *plot = popen("gnuplot -persist", "w");
    if (!plot) { perror("popen"); return; }
    fprintf(plot, "set terminal pdfcairo\n");
    fprintf(plot, "set output '%s'\n",outfile);
    fprintf(plot, "set xlabel 'Time'\n");
    fprintf(plot, "set ylabel 'Distance'\n");
    fprintf(plot, "set log\n set grid\n");
    fprintf(plot, "plot x**1.5 title 't^{1.5}', '%s' using 1:2 title 'Squared distance in time'\n", datafile);
    fprintf(plot,"unset log\nset out\n");
    fprintf(plot, "set output\n");
    pclose(plot);
}
/*******************************************************************************************************/
void printprob(const char *datafile, const char *outfile){
    FILE *plot = popen("gnuplot -persist", "w");
    if (!plot) { perror("popen"); return; }
    fprintf(plot, "set terminal pdfcairo\n");
    fprintf(plot, "set output '%s'\n",outfile);
    fprintf(plot, "set xlabel 'Walk length'\n");
    fprintf(plot, "set ylabel 'Probability over simulations'\n");
    fprintf(plot, "plot '%s' using 1:2 title 'Probability of walk of length n'\n", datafile);
    fprintf(plot, "set output\n");
    pclose(plot);
}
