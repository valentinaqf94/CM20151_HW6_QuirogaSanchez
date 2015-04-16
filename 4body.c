#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//para efectos de este script G = 1


#define A 20.0
#define B 1.0
#define C 30.0
#define D 1.0

#define h 0.001
#define min_t 0
#define max_t 1

float func_x_prime(float t, float x, float y);
float func_y_prime(float t, float x, float y);
void runge_kutta_4_step(float t_old, float x_old, float y_old, float *t_new, float *x_new, float *y_new);
void runge_kutta_4(float t0, float x0, float y0, char *filename);

int main(int argc, char **argv){

  if(argc != 3){
    printf("We need 2 arguments beside the name of the executable: x0 and y0\n");
    exit(1);
  }

  float t0 = 0.0;
  float x0 = atof(argv[1]);
  float y0 = atof(argv[2]);

  char filename[50] = "poblaciones_";
  strcat(filename, argv[1]);
  strcat(filename, "_");
  strcat(filename, argv[2]);
  strcat(filename,".dat");

  runge_kutta_4(t0, x0, y0, filename);

  return 0;
}

float func_x_prime(float t, float x, float y){
  return A*x - B*x*y;
}

float func_y_prime(float t, float x, float y){
  return -C*y + D*x*y;
}

void runge_kutta_4_step(float t_old, float x_old, float y_old, float *t_new, float *x_new, float *y_new){
  
  float kx1 = func_x_prime(t_old, x_old, y_old);
  float ky1 = func_y_prime(t_old, x_old, y_old);


  float t1 = t_old + h/2.0;
  float x1 = x_old +(h/2.0)*kx1;
  float y1 = y_old +(h/2.0)*ky1;

  float kx2 = func_x_prime(t1, x1, y1);
  float ky2 = func_y_prime(t1, x1, y1);
  

  float t2 = t_old + h/2.0;
  float x2 = x_old + (h/2.0)*kx2;
  float y2 = y_old + (h/2.0)*ky2;

  float kx3 = func_x_prime(t2, x2, y2);
  float ky3 = func_y_prime(t2, x2, y2);

  
  float t3 = t_old + h;
  float x3 = x_old + h*kx3;
  float y3 = y_old + h*ky3;

  float kx4 = func_x_prime(t3, x3, y3);
  float ky4 = func_y_prime(t3, x3, y3);


  float kx = (1.0/6.0)*(kx1 + 2.0*kx2 + 2.0*kx3 + kx4);
  float ky = (1.0/6.0)*(ky1 + 2.0*ky2 + 2.0*ky3 + ky4);

  *t_new = t_old + h;
  *x_new = x_old + h*kx;
  *y_new = y_old + h*ky;
}
  
void runge_kutta_4(float t0, float x0, float y0, char *filename){
  
  int n = (int) ((max_t - min_t)/h);
  
  float *t = malloc(sizeof(float)*n);
  float *x = malloc(sizeof(float)*n);
  float *y = malloc(sizeof(float)*n);

  float t_new, x_new, y_new;

  t[0] = t0;
  x[0] = x0;
  y[0] = y0;

  FILE *out = fopen(filename, "w");
  fprintf(out, "%f %f %f\n", t0, x0, y0);

  int i;

  for(i=1;i<n;i++){
    runge_kutta_4_step(t[i-1], x[i-1], y[i-1], &t_new, &x_new, &y_new);
    t[i] = t_new;
    x[i] = x_new;
    y[i] = y_new;
    fprintf(out, "%f %f %f\n", t_new, x_new, y_new);
  }
  fclose(out);
}
