#include <eigen/Eigen>
#include <unsupported/Eigen/FFT>
#include <Optimization/Problems/fftw3.h>


void FFT2DR2R(double *f, double *F, int width, int height){
  fftw_plan p = fftw_plan_r2r_2d(width,height,f,F,FFTW_RODFT01, FFTW_RODFT01,0);
  //  fftw_plan_dft_2d(width,height,f,F,FFTW_FORWARD,FFTW_ESTIMATE);
  //fftw_plan p = fftw_plan_dft_2d(width, height, f, F, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}


void FFT2D(fftw_complex *f, fftw_complex *F, int width, int height){
  fftw_plan p = fftw_plan_dft_2d(width, height, f, F, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void iFFT2D(fftw_complex *f, fftw_complex *F, int width, int height){
  fftw_plan p = fftw_plan_dft_2d(width, height, f, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

#define N 10
int main(){

    //const int N = 10;
    fftw_complex in[N][N], out[N][N];

    for (int i =0; i < N; i++){
        in[i][i][0] = 1.0*i;
        in[i][i][1] = 0;
    }

    //std::cout << "Matrix going in: " << std::endl;
    for (int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][j][0], out[i][j][1]);

    FFT2D(*in,*out,N,N);
    for (int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][j][0], out[i][j][1]);
}
