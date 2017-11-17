//////////////////////////////////////////////////////////
// 
// 
// 
// 
//////////////////////////////////////////////////////////

#ifndef makePlots_function_h
#define makePlots_function_h

#include <iostream>
#include <vector>

void Average_Calculation(vector<std::pair<double,double> >** vec_,vector<double> &avg)
{
    for(int i = 0;i < 7;i++)
    {
        vector<std::pair<double,double> >::iterator iter = vec_[i] -> begin();
        vector<std::pair<double,double> >::iterator iend = vec_[i] -> end();
        std::cout<<"vec_ size :"<<vec_[i]->size()<<endl;
        double sum = 0;
        double average = 0;
        for (iter;iter != vec_[i] -> end();iter++)
        {
            sum = (iter -> second) + sum;
            //cout << "the value = "<<(iter->second)<<endl;    
            //cout << "sum = "<<sum<<endl;
        }
        average = sum/(vec_[i] -> size());
        cout<<"average of ring "<<i+1<<" = "<<average<<endl;
        avg.push_back(average);
    }
}

void error_on_mean (vector <std::pair<double,double> >**vec_,vector<double> avg ,vector<double> &err)
{
    for (int i = 0; i < 7; i++)
    {
        vector<std::pair<double,double> >::iterator iter = vec_[i] -> begin();
        vector<std::pair<double,double> >::iterator iend = vec_[i] -> end();
        double average = avg [i];
        cout<<"average = "<<average<<endl;
        double var = 0;
        double error_on_mean = 0;
        for (iter; iter != vec_[i] ->end();iter++)
        {
            var = var + pow(((iter->second)-average),2);
        } 
        cout<<"var: "<<var<<endl;
        std::cout<<"vec_ size :"<<vec_[i]->size()<<endl;
        var = var /(vec_[i]->size()-1);  
        error_on_mean = sqrt(var)/sqrt(vec_[i]->size());
        cout<<"error_on_mean of ring"<<i+1<<" = "<<error_on_mean<<endl;
        err.push_back(error_on_mean);
    }

}
#endif
