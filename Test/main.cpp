//#include "csv.h"
#include "../msa.h"
#include "../Gnuplot.h"

#include <iostream>
#include <cctype>

using namespace std;
using namespace msa;


int main() {
    std::ifstream::sync_with_stdio(false);
    std::ofstream::sync_with_stdio(false);

    MsaCalc calc;
//    calc.dump("Test/example.csv");
//    calc.plot();
    calc.load("Test/example.csv");
//    calc.load();
    calc.print("report.csv");

    return 0;
}
