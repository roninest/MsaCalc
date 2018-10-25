//#include "csv.h"
#include "../msa.h"
#include "../Gnuplot.h"

#include <iostream>
#include <cctype>

using namespace std;
using namespace msa;


int main() {
    MsaCalc calc;
    calc.dump("Test/example.csv");
    calc.load();
    calc.corr();
    calc.printCorr("correlation.csv");
    calc.printMSA("msa.csv");
//    calc.plot();

    return 0;
}
