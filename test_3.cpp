#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <iostream>

double**** init_matrix(int d1, int d2, int d3, int d4) // function init_matrix reserves memory for a d1*d2*d3*d4 matrix on the heap and initializes the matrix values with 0.0
{
  double**** matrix = (double****) malloc(sizeof(double***) * d1);
  for (int i1 = 0; i1 < d1; i1++) {
    matrix[i1] = (double***) malloc(sizeof(double**) * d2);
    for (int i2 = 0; i2 < d2; i2++) {
      matrix[i1][i2] = (double**) malloc(sizeof(double*) * d3);
      for (int i3 = 0; i3 < d3; i3++) {
        matrix[i1][i2][i3] = (double*) malloc(sizeof(double) * d4);
        for (int i4 = 0; i4 < d4; i4++) {
          matrix[i1][i2][i3][i4] = .0;
        }
      }
    }
  }
  return matrix;
}


void delete_matrix(double**** matrix, int d1, int d2, int d3) // frees up the allocated memory of the matrix when we don't need it anymore
 {
  for (int i1 = 0; i1 < d1; i1++) {
    for (int i2 = 0; i2 < d2; i2++) {
      for (int i3 = 0; i3 < d3; i3++) {
        free(matrix[i1][i2][i3]);
      }
      free(matrix[i1][i2]);
    }
    free(matrix[i1]);
  }
  free(matrix);
}

double distance(int x0a, int x1a, int x2a, int x3a, int y0a, int y1a, int y2a, int y3a) //function distance calculates the distance between two points in 4d, and takes the coordiantes of the first point as integers first, the ones of the other one second
{
    double x0 = (int) x0a;
    double x1 = (int) x1a;
    double x2 = (int) x2a;
    double x3 = (int) x3a;
    double y0 = (int) y0a;
    double y1 = (int) y1a;
    double y2 = (int) y2a;
    double y3 = (int) y3a;
       
    double sq0 = pow(x0-y0, 2.0);
    double sq1 = pow(x1-y1, 2.0);
    double sq2 = pow(x2-y2, 2.0);
    double sq3 = pow(x3-y3, 2.0);
    return (sqrt(sq0+sq1+sq2+sq3));
}


void hf (double**** G, int i, int j, int k, int l){                         //defining the 4dimensional heatfunction

    G[i][j][k][l] = G[i][j][k][l] + 0.123 * ((-8.0)*G[i][j][k][l] + G[i+1][j][k][l]
     + G[i-1][j][k][l] + G[i][j+1][k][l] + G[i][j-1][k][l] + G[i][j][k+1][l]
      + G[i][j][k-1][l] + G[i][j][k][l+1] + G[i][j][k][l-1]);

}

void hf_iterator(double**** matrix, int size){
  for (int index1 = 1; index1<size-1; index1++){
    for (int index2 = 1; index2<size-1;index2++){
      for (int index3 = 1; index3<size-1;index3++){
        for (int index4 = 1; index4<size-1;index4++){
          hf(matrix, index1, index2, index3, index4);
}}}}}

void heater (int i, int j, int k , int l, double dist, double heat, int size, double**** matrix){  //first four parameters are the coordinates of the middle of the heatsphere, dist is the radius of the heatsphere, heat is its heat, size is the length of each matrix dimension and matrix is the Matrix itself
  for (int index1 = 0; index1<size; index1++){
    for (int index2 = 0; index2<size;index2++){
      for (int index3 = 0; index3<size;index3++){
        for (int index4 = 0; index4<size;index4++){
          if (distance(i,j,k,l,index1, index2, index3, index4)<dist){
            matrix[index1][index2][index3][index4]=heat;
            }
        }
      }
    }
  }
}


int main(int argc, char** argv) {

  auto begin = std::chrono::high_resolution_clock::now();


  /*if (argc < 5) {
    printf("Specify a size for every dimension!\n");
    exit(1);
  }*/

  int number_of_iterations = atoi(argv[2]);
  int heatsourcepos1 = atoi(argv[3]);
  int heatsourcepos2 = atoi(argv[4]);
  int heatsourcepos3 = atoi(argv[5]);
  int heatsourcepos4 = atoi(argv[6]);
  double heat = atoi(argv[7]);
  double radius = atoi(argv[8]);
  int d1 = atoi(argv[1]);
  int d2 = atoi(argv[1]);
  int d3 = atoi(argv[1]);
  int d4 = atoi(argv[1]);
  double**** G = init_matrix(d1, d2, d3, d4);

  //printf("%f \n", G[3][3][3][3]);

  heater(heatsourcepos1, heatsourcepos2, heatsourcepos3, heatsourcepos4, radius, heat, d1, G);

  //printf("%f \n", G[3][3][3][3]);

  for (int iterator = 0; iterator<number_of_iterations; iterator++){
    hf_iterator(G, d1);
  }

  //printf("%f \n", G[3][3][3][3]);

  delete_matrix(G, d1, d2, d3);

  auto end = std::chrono::high_resolution_clock::now();

  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << "ns" << std::endl;


}