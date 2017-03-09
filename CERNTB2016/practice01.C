#include <iostream>

double** Func(double** layer, int Nlayer, int Ncell)
{
    for(int i=0;i<Nlayer;++i)
        for(int j=0;j<Ncell;++j)
            layer[i][j] += 1;
    return layer;
}


int main()
{
    //double layer[100][7]={0};
    double** layer;
    layer = new double*[100];
    for(int i=0;i<100; ++i)
        layer[i] = new double[7];
    for(int i=0;i<100;++i)
        for(int j=0;j<7;++j)
            layer[i][j] = 1;
    
    double** LLL;
    LLL = Func(layer, 100, 7);

    for(int i=99;i>-1;--i)
        delete[] layer[i];
    delete[] layer;

    double* a = new double();
    delete a;

    double* a = new double[7];
    delete[] a;
}

