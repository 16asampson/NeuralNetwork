#include<iostream>
#include<fstream>
#include<cmath>
#include<list>
#include<stdlib.h> //srand rand
#include<time.h>    //time
#include<vector>

#include "inputs.h"
using namespace std;

class Layer
{   
    //variables
    int nu = 4; //num units
    int ni = 3; //num inputs
            //[ni+1][nu]
    //int W[4][4]; //W[k][j] is the weight applied to inout k for unit j
    vector<vector<double>> W;
    vector<double> input; //input vector
    vector<double> output; //output vector
    
    public:
    //variables
    //double deltaW[4][4]; //unit's weight change
    vector<vector<double>> deltaW;
    double eta;
    vector<double> delta_prev; //error term of inputs

    void init(int units, int inputs, double LRate)
    {
        nu = units;
        ni = inputs;
        eta = LRate;

        //make a grid of size [ni+1][nu]
        vector<vector<double>> temp_weights((ni+1),vector<double>(nu));
        vector<vector<double>> temp_deltaW((ni+1),vector<double>(nu));

        srand(time(NULL));
        //randomize weights
        for (int x=0; x<temp_weights.size();x++)
        {
            for(int y=0; y<temp_weights[x].size();y++)
            {
                temp_weights[x][y] = rand() % 100;
                temp_weights[x][y] = temp_weights[x][y]/100;
            }
        }
        srand(1);
        

        //set size of output vector
        vector<double> temp_output(nu);
        output = temp_output;

        //set size of input vector
        vector<double> temp_input(ni+1);
        input = temp_input;

        //set grids as the global placeholders
        W = temp_weights;
        deltaW = temp_deltaW;

        //set size of delta_previous
        vector<double> temp_delta_prev(ni);
        delta_prev = temp_delta_prev;

        //output random weights
        //cout << "Random weights" << endl;
        //print();

    }

    //print weight table for layer
    void print()
    {
        
        //loop through all units
        for (int u=0; u<nu;u++)
        {
            //loop through all inputs
            for (int x=0; x<(ni+1); x++)
            {
                cout << W[x][u] << "    ";
            }
            cout << endl;
        }
        
        

    }


    //function feed forward
    vector<double> feedForward( vector<double> inp)
    {
        input = inp;
        //this may cause error due to size difference
        //input[ni] = 1;

        input.push_back(1);
        //PrintVector(input);
        for (int u =0; u < nu; u++)
        {
            double temp =0;
            for (int k=0; k<ni; k++)
            {
                temp = temp + W[k][u] * input[k];

            }
            //outputs are 1 and 2 so increment up 1
            output[u] = sigmoid(temp);
        }
        return output;

    }

    vector<double> backProp(vector<double> delta)
    {
        int hold;

        //calculate delta of previous layer
        for (int j=0; j<ni; j++)
        {
            double temp = 0;
            for (int u=0; u<nu; u++)
            {
                temp = temp + delta[u] * W[j][u];
            }
            delta_prev[j] = temp * input[j] * (1- input[j]);

        }
        //change weights for all units in layer
        for (int j=0; j < ni; j++)
        {
            for (int u=0; u<nu; u++)
            {

                deltaW[j][u] = eta*delta[u]*output[u];

                W[j][u] = W[j][u] + deltaW[j][u];

            }
        }
        return delta_prev;
    }
    
    


};

//return a trained Nural network
vector<Layer> Backpropagation(double LearningRate)
{
    //TODO: read input document
    vector<Example> ExList;
    string fExamples = "./hw3data.txt";
    ExList = readInDataList(fExamples);


    //TODO: get LearningRate from user
    double LRate = LearningRate;
    //cin >> LRate;

    //create a nueral network with 3 layers
    vector<Layer> NN(2);

    //initiallize Layers in NN
    NN[0].init(4,2,LRate);
    NN[1].init(1,4,LRate);

    //vector of Ycap
    vector<double> YCap(ExList.size());

    

    //repeate backpropogation for some time
    for (int count=0; count<10000; count++)
    {
        //for all examples 
        for (int e=0;e<ExList.size(); e++)
        {
            vector<double> values = ExList[e].x;
            
            //for each layer in NN feed forward from input to output
            for(int l=0; l<NN.size(); l++)
            {
                values = NN[l].feedForward(values);
            }

            YCap = values;

            vector<double> delta;

            //calculate delta of output layer
            delta.push_back(2*(ExList[e].y[0] - YCap[0])* YCap[0] * (1-YCap[0]));


            //for each layer in NN from output to input layer
            for (int l=NN.size()-1; l>-1;l--)
            {
                //cout << "Back Propogate" << endl;
                delta = NN[l].backProp(delta);
            }
            
        }
    }
    //return the trained NN
    return NN;

}

//returns sum of squares error from validation data
double sumOfSquares(vector<Layer> NN)
{
    //TODO: read input document
    vector<Example> ExList;
    string fExamples = "./hw3valid.txt";
    ExList = readInDataList(fExamples);

    //vector of Ycap
    vector<double> YCap;

    //loop through all examples in list (validation data)
    //get YCap Values
    for (int e=0; e< ExList.size(); e++)
    {
        vector<double> values = ExList[e].x;
            
        //for each layer in NN feed forward from input to output
        for(int l=0; l<NN.size(); l++)
        {
            values = NN[l].feedForward(values);
        }

        //store YCap of every example
        YCap.push_back(values[0]);
    }
    //cout << "Sum of YCap = ";
    //PrintVector(YCap);


    //calculate sum of squares error
    double Error =0;
    for (int e=0; e<ExList.size(); e++)
    {
        double temp = ExList[e].y[0] - YCap[e];
        //cout << "differnece" << temp << endl;
        Error = Error + pow(temp,2);
    }

    return Error;
}

int main()
{
    double LRate = 0;
    //get Learning Rate eta
    cin >> LRate;

    //make and train neural network
    vector<Layer> NN = Backpropagation(LRate);

    cout << "CS-5001: HW3" << endl;
    cout << "Programmer: Austin Sampson" << endl;
    //print inputs weights table
    cout << "TRAINING" << endl;
    cout << "Using learning rate eta = " << LRate << endl;
    cout << "Using 10000 iterations " << endl;

    cout << "OUTPUT" << endl;
    cout << "input Layer" << endl;
    NN[0].print();
    cout << "output Layer" << endl;
    NN[1].print();
    
    cout << "Validation" << endl;
    cout << "Sum-of-Squares Error = " << sumOfSquares(NN) << endl;


    return 0;
}