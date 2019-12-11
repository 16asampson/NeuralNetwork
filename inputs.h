#ifndef INPUTS_H
#define INPUTS_H

#include<vector>  //process vector
#include<fstream>  // file i/o
#include<iostream> // cerr
#include<list>
#include<cmath>

using namespace std;


struct Example
{
    // x1, x2, x0
    vector<double> x = {0, 0};
    vector<double> y;
};

//purpose: to read in the input cordinates data from the training file
//and input it into a list for processing by the Linear Learner
//return: dataList
vector<Example> readInDataList(const string& fname)
{
    vector<Example> dataList;
    fstream myfile;
    int num;
    //open the training file
    myfile.open(fname);

    if(myfile.fail())
    {
        cerr << "Unable to open file \"" << fname << "\", terminating" << endl;
        exit(-1);
    }

    //get all input from file and put in cordniate list
    while(!myfile.eof())
    {
        Example temp;
        double tempY;
        myfile >> temp.x[0] >> temp.x[1] >> tempY;
        temp.y.push_back(tempY);
        temp.y[0] = temp.y[0];
        dataList.push_back(temp);
    }
    
    myfile.close(); //close the training file
    return dataList;

}

//inplementation of sigmoid function
double sigmoid(double input)
{
    double temp = 1/(1+ exp(-input));
    //cout << "sig input = " << input << endl;
    //cout << "sigmoid = " << temp << endl;
    return temp;
}

//temp functions
//print vector
void PrintVector(vector<double> temp)
{
    for (int x=0; x< temp.size(); x++)
    {
        cout << temp[x] << "    ";
    }
    cout << endl;
}


#endif