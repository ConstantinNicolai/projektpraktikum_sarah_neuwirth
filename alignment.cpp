#include <cmath>
#include <chrono>
#include <iostream>
#include <omp.h>


    struct alignas(64) slot{
          double a;
      };

 

  int main(int argc, char** argv){

      //omp_set_num_threads(atoi(argv[2]));

      int count = atoi(argv[1]);
      //slot* b = (slot*) malloc(sizeof(slot) * count);
      double* b =  (double*) malloc(sizeof(double)*count);

      slot** array = (slot**) malloc(sizeof(slot*)*5);
      for (int j= 0; j<5; j++){
          array[j] = (slot*) malloc(sizeof(slot) * count);
      }
      //slot* c = (slot*) malloc(sizeof(slot) * count);

      auto begin = std::chrono::high_resolution_clock::now();
  
      for (int i =0; i<5; i++){
      for(int i=0; i<count; i++){
          b[i]=pow(i,20);
      }}
      
      auto end = std::chrono::high_resolution_clock::now();

      auto begin1 = std::chrono::high_resolution_clock::now();

      #pragma omp parallel for collapse(2)
      for (int h =0; h<5; h++){ 
      //#pragma omp parallel for
      for(int i=0; i<count; i++){
          array[h][i].a=pow(i,20);
      }}

      
      auto end1 = std::chrono::high_resolution_clock::now();

      std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << "ns" << std::endl;
      std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end1-begin1).count() << "ns" << std::endl;


      

      

      
    return 0;
  }