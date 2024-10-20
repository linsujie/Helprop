#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include<thread>
#include "particle.h"
#include<mutex>
#include"import.h"
#include<string>
#include<sstream>



std::mutex mtx;







void import(int &argc, char* argv[], std::vector<double> &ekin, std::vector<double> &spec, std::vector<double> &param, std::vector<int> &seq){
    std::cout << argc << "   2:  " << *argv[2] << "   yes" << std::endl;

    std::string inputFileName;
    std::string spectrumFileName;
    for(int i=1 ; i<argc ; i++){
        if(argv[i][0] == '-' && argv[i][1] == 'i'){
            inputFileName = argv[i+1];
            // std::cout << inputFileName << std::endl;
        }
        else if(argv[i][0] == '-' && argv[i][1] == 's'){
            spectrumFileName = argv[i+1];
            // std::cout << spectrumFileName << std::endl;
        }
    }

    std::string line;

    std::ifstream inputFile(inputFileName);
    while(std::getline(inputFile, line)){
        std::string x;
        double y;
        std::istringstream iss(line);
        iss >> x >> y;
        if(x == "E0") seq.push_back(1), param.push_back(y);
        else if(x == "B0") seq.push_back(2), param.push_back(y);
        else if(x == "polarity") seq.push_back(3), param.push_back(y);
        else if(x == "angle") seq.push_back(4), param.push_back(y);
        else if(x == "D") seq.push_back(5), param.push_back(y);
        else if(x == "indexA") seq.push_back(6), param.push_back(y);

    }
    
    std::ifstream spectrumFile(spectrumFileName);
    while(std::getline(spectrumFile, line)){
        double x, y;
        std::istringstream iss(line);
        iss >> x >> y;
        ekin.push_back(x);
        spec.push_back(y);
    }

    // std::cout << "ekin:   " << ekin.size() << std::endl;




}


int main(int argc, char* argv[]){


    std::vector<double> ekin;                                               //set spectrum energy bin
    std::vector<double> spec;                                               //boundary differential flux
    std::vector<double> param;                                              //parameters vector
    std::vector<int> seq; 

    import(argc, argv, ekin, spec, param, seq);                //import parameters and spectrum


    std::vector<std::vector<double>> weight;                                //possibility matrix

    int number = 3000;                                                      //test particle number
    int th_num = 1;                                                        //thread count
    for(int i=0;i<ekin.size();i++){
            seq.push_back(7);
            param.push_back(ekin[i]);                                       //get kinetic energy
            std::vector<particle> Particle;
            for(int j=0;j<number/th_num;j++){
                std::vector<std::thread> threads;
                for(int h=0;h<th_num;h++){
                    particle one(param, seq);
                    threads.emplace_back([one, &Particle]() mutable {
                        one.step();
                        std::lock_guard<std::mutex> lock(mtx);
                        Particle.push_back(one);
                    });
                }
                for(int j=0;j<th_num;j++){
                    threads[j].join();
                }
            }

            std::vector<double> bin;
            bin.resize(ekin.size());
            for(int j=0;j<number;j++){
                double eng = Particle[j].Ek;
                for(int k=0;k<ekin.size();k++){
                    if(k==0){
                        double x1 = log(ekin[k+1])/2. - log(ekin[k])/2.;
                        if(log(ekin[k])-x1<=log(eng) && log(eng)<log(ekin[k])+x1) bin[k] += 1./number;
                    }
                    else if(0<k && k<ekin.size()-1){
                        double x0 = log(ekin[k-1])/2. + log(ekin[k])/2.;
                        double x1 = log(ekin[k+1])/2. + log(ekin[k])/2.;
                        if(x0<=log(eng) && log(eng)<x1) bin[k] += 1./number;
                    }
                    else if(k == ekin.size()-1){
                        double x1 = log(ekin[k])/2. - log(ekin[k-1])/2.;
                        if(log(ekin[k])-x1<=log(eng) && log(eng)<log(ekin[k])+x1) bin[k] += 1./number;
                    }
                }
            }
            weight.push_back(bin);
    }


    std::vector<double> Ospec;
    for(int i=0;i<weight.size();i++){
        double value = 0;
        for(int j=0;j<spec.size();j++){
            value += spec[j]*weight[i][j];
        }
        Ospec.push_back(value);
    }


    std::cout << "done" << std::endl;






    

    return 0;
}