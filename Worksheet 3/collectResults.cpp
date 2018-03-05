#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

void analyseResults(int N, const std::vector<double>& compTime, const std::vector<double>& mpiTime) {
    double averageComp(0);
    double averageMPI(0);
    double nStats(compTime.size());
    for (size_t i(0); i<compTime.size(); ++i) {
        averageComp+=compTime[i];
        averageMPI+=mpiTime[i];
    }
    averageComp/=nStats;
    averageMPI/=nStats;
    double varianceComp(0);
    double varianceMPI(0);
    for (size_t i(0); i<compTime.size(); ++i) {
        varianceComp+=(compTime[i]-averageComp)*(compTime[i]-averageComp);
        varianceMPI+=(mpiTime[i]-averageMPI)*(mpiTime[i]-averageMPI);
    }
    varianceComp/=(nStats-1);
    varianceMPI/=(nStats-1);
    std::cout << std::setw(13) << std::left << N << "   ";
    std::cout << std::setprecision(6) << std::setw(13) << std::left << averageComp << "   ";
    std::cout << std::setprecision(6) << std::setw(13) << std::left << varianceComp << "   ";
    std::cout << std::setprecision(6) << std::setw(13) << std::left << averageMPI << "   ";
    std::cout << std::setprecision(6) << std::setw(13) << std::left << varianceMPI << std::endl;
}

int main(int argc,char** argv) {
    if (argc==1) {
        std::cerr << "Missing file input" << std::endl;
        return 1;
    }
    std::ifstream runOutput(argv[1]);
    if (!runOutput.good()) {
        std::cerr << "Open file failed" << std::endl;
        return 1;
    }
    int N(0);
    std::vector<double> compTime;
    std::vector<double> mpiTime;
    std::cout << std::setw(13) << std::left << "N" << "   " << std::setw(13) << std::left << "averageComp" << "   ";
    std::cout << std::setw(13) << std::left << "varianceComp" << "   ";
    std::cout << std::setw(13) << std::left << "averageMPI" << "   " << std::setw(13) << std::left << "varianceMPI" << std::endl;
    
    while (!runOutput.eof()) {
        
        std::string line;
        
        //Date
        std::getline (runOutput,line);
        if (runOutput.eof()) break;
        if (line.size()==0) continue;
        
        // mpiexec -n 64 ./cannon XX YY
        // matrix size
        std::getline (runOutput,line);
        if (line.size()==0) continue;
        int i(line.rfind(','));
        int newN(stoi(line.substr(i+1,line.size()-i-2)));
        if (runOutput.eof()) break;
        
        if (N!=0 && newN!=N) {
            analyseResults(N,compTime,mpiTime);
            compTime.clear();
            mpiTime.clear();
        }
        N=newN;
        
        // computation time
        std::getline (runOutput,line);
        if (line.size()==0) continue;
        i=line.rfind(' ');
        compTime.push_back(stod(line.substr(i+1,line.size()-i-1)));
        if (runOutput.eof()) break;
        
        // MPI time
        std::getline (runOutput,line);
        if (line.size()==0) continue;
        i=line.rfind(' ');
        mpiTime.push_back(stod(line.substr(i+1,line.size()-i-1)));
        if (runOutput.eof()) break;
        
    }
    analyseResults(N,compTime,mpiTime);
    runOutput.close();
    return 0;
}
