#define ENABLE_OPENMP
#include "../msa.h"
#include "../Gnuplot.h"

using namespace std;

int main() {
    msa::MsaCalc calc;
    cout << "Dump start" << '\n';
    calc.dump("Test/example.csv");
    cout << "Dump end" << '\n';
//    cout << "Load start" << '\n';
//    calc.load();
//    cout << "Load end" << '\n';
    cout << "Corr start" << '\n';
    calc.correlation();
    cout << "Corr end" << '\n';
    cout << "Print start" << '\n';
    calc.printCorr("correlation.csv");
    calc.printMSA("msa.csv");
    cout << "Print end" << '\n';
//    calc.plot();
}
