//COMMANDLINE PARAMETERS: int matrixsize; int number of iterations; int heatsource 1; int heatsource 2; int heatsource 3; int heatsource 4; double heat; double radius heatsource



#include <cmath>
#include <chrono>
#include <iostream>
#include <omp.h>


double** init_matrix(int d1, int d2); //reserves heap memory for the matrix and ititializes the matrix to 0.0

void delete_matrix(double** matrix, int d1); //frees up the heap memory for the matrix when it isn't needed any longer to avoid memory leaks

double distance(int x0a, int x1a, int y0a, int y1a); //calculates the distance between point x and y with their respective 4d positions

void hf (double** &G,double** &matrix, int i, int j); //runs the defined heatfuntion on entry  i,j,k,l of the given matrix

void hf_iterator(double** &matrix,double** &matrix2, int size);  //iterates over the whole matrix except the boundaries and calls the heatfunction

void hf_flipflop(double** &matrix,double** &matrix2, int size, int number_of_iterations);

void heater (int i, int j, double dist, double heat, int size, double** matrix); //"heats up" all points in a radius "distance" around point i,j,k,l to "heat", 
                                                                                                //size is the matrixlength in one dimension, matrix specifies the matrix to be altered
int main(int argc, char** argv) {



  if (argc < 7) {                                    //checks whether we have enough commandline paramters given
    printf("Insufficient input parameters!\nCOMMANDLINE PARAMETERS: integer matrixsize, integer number of iterations,");
    printf("integer heatsource 1, integer heatsource 2, integer heatsource 3, integer heatsource 4, doubleprecision floating point heat, doubleprecision floating point radius heatsource \n");
    exit(1);
  }
  
  auto begin = std::chrono::high_resolution_clock::now();

  int number_of_iterations = atoi(argv[2]);
  int heatsourcepos1 = atoi(argv[3]);
  int heatsourcepos2 = atoi(argv[4]);
  double heat = atoi(argv[5]);
  double radius = atoi(argv[6]);
  int d1 = atoi(argv[1]);
  int d2 = atoi(argv[1]);
  double** G = init_matrix(d1, d2);
  double** P = init_matrix(d1, d2);
  //int size = atoi(argv[1]);
  //int threadnumber = atoi(argv[7]);

  //omp_set_num_threads(threadnumber);

  heater(heatsourcepos1, heatsourcepos2, radius, heat, d1, G);
  heater(heatsourcepos1, heatsourcepos2, radius, heat, d1, P);

  hf_flipflop(G, P, d1, number_of_iterations);

  int size = d1;

  if (number_of_iterations%2 == 0){

  for (int index3 = 0; index3<size;index3++){
          printf("\n");
        for (int index4 = 0; index4<size;index4++){
            printf("%f  ", G[index3][index4]);}}}
  else {for (int index3 = 0; index3<size;index3++){
          printf("\n");
        for (int index4 = 0; index4<size;index4++){
            printf("%f  ", P[index3][index4]);}}}

  delete_matrix(G, d1);
  delete_matrix(P, d1);

  auto end = std::chrono::high_resolution_clock::now();

  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << "ns" << std::endl;
}

double** init_matrix(int d1, int d2) // function init_matrix reserves memory for a d1*d2*d3*d4 matrix on the heap and initializes the matrix values with 0.0
{
  double** matrix = (double**) malloc(sizeof(double*) * d1);
  for (int i1 = 0; i1 < d1; i1++) {
    matrix[i1] = (double*) malloc(sizeof(double) * d2);
        for (int i2 = 0; i2 < d2; i2++) {
          matrix[i1][i2] = .0;
        }
      }
  return matrix;
}

void delete_matrix(double** matrix, int d1) // frees up the allocated memory of the matrix when we don't need it anymore
 {
  for (int i1 = 0; i1 < d1; i1++) {
    
    free(matrix[i1]);
  }
  free(matrix);
}

double distance(int x0a, int x1a, int y0a, int y1a) //function distance calculates the distance between two points in 4d, and takes the coordiantes of the first point as integers first, the ones of the other one second
{
    double x0 = (int) x0a;
    double x1 = (int) x1a;
    double y0 = (int) y0a;
    double y1 = (int) y1a;
       
    double sq0 = pow(x0-y0, 2.0);
    double sq1 = pow(x1-y1, 2.0);
    return (sqrt(sq0+sq1));
}

void hf (double** &G,double** &matrix, int i, int j){                         //defining the 4dimensional heatfunction

    matrix[i][j] = G[i][j] + 0.24 * ((-4.0)*G[i][j] 
    + G[i+1][j] + G[i-1][j] + G[i][j+1] + G[i][j-1]
     );

}

void hf_iterator(double** &matrix,double** &matrix2, int size){
//#pragma omp parallel
{
  for (int index1 = 1; index1<size-1; index1++){
    //#pragma omp parallel
    for (int index2 = 1; index2<size-1;index2++){
          hf(matrix,matrix2, index1, index2);
          }}
          }
          }

void hf_flipflop(double** &matrix,double** &matrix2, int size, int number_of_iterations){
    for (int t = 0; t < number_of_iterations; t++){
        if ((t%2) == 0){
            hf_iterator(matrix, matrix2, size);
        }
        else{
            hf_iterator(matrix2, matrix, size);
        }
        
    }}

void heater (int i, int j, double dist, double heat, int size, double** matrix){  //first four parameters are the coordinates of the middle of the heatsphere, dist is the radius of the heatsphere, heat is its heat, size is the length of each matrix dimension and matrix is the Matrix itself
  for (int index1 = 0; index1<size; index1++){
    for (int index2 = 0; index2<size;index2++){
          if (distance(i,j,index1, index2)<dist){
            matrix[index1][index2]=heat;
            }
        }
  }
}