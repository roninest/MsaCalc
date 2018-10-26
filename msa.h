//
//
// Created by balrog on 10/10/18.
//
#pragma once

#include "Gnuplot.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <thread>
#include <omp.h>
#include <algorithm>
#ifdef DEBUG
#include <chrono>
#endif

namespace  msa
{

enum class Method {
    CUMULATIVE,
    DIRECT
};


struct Parameter
{
    std::string name = "";
    std::string unit = "";
    double lsl = 0.0f;
    double usl = 0.0f;
    double m = 0.0f;
    double m2 = 0.0f;
    double sm2 = 0.0f;
    double sum = 0.0f;
    double dsm2 = 0.0f;
    double ks = -1.0f;
    size_t n = 0u;

    bool operator<(const Parameter &that)
    {
        return this->name < that.name;
    }
};


struct Cell
{
    double r = 0.0f;
    size_t a = 0u;
    size_t b = 0u;


    bool operator<(const Cell &that)
    {
        return std::abs(this->r) > std::abs(that.r);
    }
};


class MsaCalc
{
public:
    MsaCalc() = default;
    virtual ~MsaCalc() = default;

public:
    size_t check(const std::string &filename)
    {
        std::string line, s;
        size_t columns = 0;
        std::ifstream f;
        f.open(filename.c_str());


        if (!f.is_open()) {
            std::cerr << "Error during file opening: " << filename << std::endl;
        }
        else {
            //if Ok
        }


        if(getline(f, line)) {
            std::stringstream ss(line);
            while(getline(ss, s, ',')) {
                columns++;
            }
        }
        else {
            std::cerr << "File " << filename << " is empty" << std::endl;
        }

        if (f.is_open()) {
            f.close();
        }

        return columns;
    }

    size_t fillSpec(const std::string &filename)
    {
        std::string line, s;
        std::ifstream f;

        size_t columns = check(filename);
        size_t curr = 0;

        f.open(filename);
        if (!f.is_open()) {
            std::cerr << "Error during file opening: " << filename << std::endl;
        }
        else {
            //if Ok
        }

        if (!parameters.empty()) {
            parameters.clear();
        }

        parameters.resize(columns);

        if(getline(f, line)) {
            std::stringstream ss(line);
            curr = 0;
            while(getline(ss, s, ',')) {
                parameters[curr].name = s;
                curr++;
            }

            if (curr != columns) {
                std::cerr << "Row mismatch during 2ns read" << std::endl;
            }
            else {
                //is Ok
            }
        }
        else {
            std::cerr << "File " << filename << " is empty" << std::endl;
        }


        if(getline(f, line)) {
            std::stringstream ss(line);
            curr = 0;
            while(getline(ss, s, ',')) {
                if (!is_float(s)) {
                    curr++;
                    continue;
                }

                parameters[curr].lsl = stod(s);
                curr++;
            }

            if (curr != columns) {
                std::cerr << "Row 2 mismatch with row 1" << std::endl;
            }
            else {
                //is Ok
            }
        }
        else {
            std::cerr << "Row 2 of" << filename << " is empty" << std::endl;
        }


        if(getline(f, line)) {
            std::stringstream ss(line);
            curr = 0;
            while(getline(ss, s, ',')) {
                if (!is_float(s)) {
                    curr++;
                    continue;
                }
                parameters[curr].usl = stod(s);
                curr++;
            }

            if (curr != columns) {
                std::cerr << "Row 3 mismatch with row 1" << std::endl;
            }
            else {
                //is Ok
            }
        }
        else {
            std::cerr << "Row 3 " << filename << " is empty" << std::endl;
        }


        if(getline(f, line)) {
            std::stringstream ss(line);
            curr = 0;
            while(getline(ss, s, ',')) {
                parameters[curr].unit = s;
                curr++;
            }

            if (curr != columns) {
                std::cerr << "Row 4 mismatch with row 1" << std::endl;
            }
            else {
                //is Ok
            }
        }
        else {
            std::cerr << "Row 4 " << filename << " is empty" << std::endl;
        }

        if (f.is_open()) {
            f.close();
        }

        return columns;

    }


    void dump(const std::string &filename)
    {
        std::string line, s;
        std::ifstream f;

        size_t curr = 0;
        size_t row = 0;
        size_t total = 0;
        size_t columns = fillSpec(filename);

        f.open(filename);
        if (!f.is_open()) {
            std::cerr << "Error during file opening: " << filename << std::endl;
        }
        else {
            //if Ok
        }


        while(getline(f,line)) {
            std::stringstream ss(line);
            curr = 0;
            bool ignore = false;

            while(getline(ss, s, ',')) {
                if (ignore) {
                    curr++;
                    continue;
                }
                else if (curr == 25 && s != "PASS_BIN") {
                    curr++;
                    ignore = true;
                    continue;
                }

                if (!is_float(s)) {
                    curr++;
                    continue;
                }

                parameters[curr].n++;
                curr++;
                total++;
            }
        }

        if (f.is_open()) {
            f.close();
        }

        if (!filedump.empty()) {
            filedump.clear();
        }

        filedump.resize(parameters.size());
        for (size_t i = 0; i < parameters.size(); i++) {
            filedump[i].reserve(parameters[i].n);
        }

        f.open(filename);
        if (!f.is_open()) {
            std::cerr << "Error during file opening: " << filename << std::endl;
        }
        else {
            //if Ok
        }


        for (int i = 0; i<4; i++) {
            if (getline(f, line)) {
                row++;
            }
            else {
                std::cerr << "Line " << i+1 << " is empty in file: " << filename << std::endl;
            }
        }


        while(getline(f,line)) {
            std::stringstream ss(line);
            curr = 0;
            bool ignore = false;
            row++;

            while(getline(ss, s, ',')) {
                if (ignore) {
                    curr++;
                    continue;
                }
                else if (curr == 25 && s != "PASS_BIN") {
                    curr++;
                    ignore = true;
                    continue;
                }

                if (!is_float(s)) {
                    curr++;
                    continue;
                }

                filedump[curr].push_back(stof(s));
                curr++;
            }

            if (curr != columns) {
                std::cerr << "Row " << row << " mismatch with row 1" << std::endl;
            }
            else {
                //Ok
            }
        }


        if (f.is_open()) {
            f.close();
        }
    }


    void plot() {
        #pragma omp parallel for
        for (size_t i = 0; i < parameters.size(); i++) {
            if (parameters[i].n > 3) {
                Gnuplot gp;
                gp.send("set terminal png");
                std::string cmd = std::string("set output '") +  parameters[i].name + ".png'";
                gp.send(cmd);
                gp.send("set grid");
                gp.send("plot '-' smooth frequency with boxes"); //with dots
                gp.plot(filedump[i]);
            }
        }
    }

    void corr() {
        if (!correlation.empty()) {
            correlation.clear();
        }

        const size_t n = parameters.size();

        correlation.reserve(( (n/2) * (n-1)/2) );

        #pragma omp parallel for
        for (size_t i = 0; i < n; i++)  {
            #pragma omp parallel for
            for (size_t j = 0; j < i; j++) {
                    if (filedump[i].size() > 10000 && filedump[j].size() > 10000) {
                        double coff = corr(filedump[i], filedump[j], 10000);
                        if (coff) {
                            double m = parameters[i].m;
                            double s = std::sqrt(parameters[i].sm2);
                            if (parameters[i].ks == -1.0f) {
                                parameters[i].ks = ks(filedump[i], m, s);
                            }
                            m = parameters[j].m;
                            s = std::sqrt(parameters[j].sm2);

                            if (parameters[j].ks == -1.0f) {
                                parameters[j].ks = ks(filedump[j], m, s);
                            }

                            correlation.push_back({coff, i, j});
                        }
                }
            }
        }
    }


    double corr(const std::vector<float> &a, const std::vector<float> &b, size_t n) {
        double sumA = 0.0f, sumB = 0.0f, sumAB = 0.0f;
        double sumA2 = 0.0f, sumB2 = 0.0f;

        #pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            // sum of elements of array X.
            sumA = sumA + a[i];

            // sum of elements of array Y.
            sumB = sumB + b[i];

//            sumAmB = sumAmB + (a[i] - b[i]);

            // sum of X[i] * Y[i].
            sumAB = sumAB + a[i] * b[i];

            // sum of square of array elements.
            sumA2 = sumA2 + a[i] * a[i];
            sumB2 = sumB2 + b[i] * b[i];
        }

        const double result = (n*sumAB - sumA*sumB) / std::sqrt( (n*sumA2 - sumA*sumA ) * (n*sumB2 - sumB*sumB) );

        return (abs(result) < 1.0f && abs(result) > 0.8f) ? result : 0.0f;
    }


    void load()
    {
        if (filedump.empty()) {
            std::cerr << "Data is empty" << std::endl;
            return;
        }


        for (size_t i = 0; i < parameters.size(); i++) {
            size_t n = 0u;
            Parameter &param = parameters[i];

            double m = param.m;
            double m2 = param.m2;
            double sm2 = param.sm2;

            if (filedump[i].size() < 10000u){
                for (const auto &it : filedump[i]) {
                    double x = it;
                    m = (x + n*m)/(n + 1);
                    m2 = (x*x + n*m2)/(n + 1);
                    sm2 = m2 - m*m;
                    n++;

                }
            }
            else {
                for (auto it = filedump[i].begin(); it != filedump[i].begin() + 10000u; ++it) {
                    double x = *it;
                    m = (x + n*m)/(n + 1);
                    m2 = (x*x + n*m2)/(n + 1);
                    sm2 = m2 - m*m;
                    n++;
                }

            }


            param.m = m;
            param.m2 = m2;
            param.sm2 = sm2;
        }
    }


    void load(const std::string &filename)
    {
        std::string line, s;
        std::ifstream f;
        size_t curr = 0;
        size_t  row = 0;
        size_t columns = fillSpec(filename);

        f.open(filename);
        if (!f.is_open()) {
            std::cerr << "Error during file opening" << std::endl;
        }
        else {
            //if Ok
        }


        for (int i = 0; i<4; i++) {
            if (getline(f, line)) {
                row++;
            }
            else {
                std::cerr << "Line " << i+1 << " is empty in file: " << filename << std::endl;
            }
        }

        while(getline(f,line)) {
            std::stringstream ss(line);
            curr = 0;
            bool ignore = false;
            row++;

            while(getline(ss, s, ',')) {
                if (ignore) {
                    curr++;
                    continue;
                }
                else if (curr == 25 && s != "PASS_BIN") {
                    curr++;
                    ignore = true;
                    continue;
                }

                if (!is_float(s)) {
                    curr++;
                    continue;
                }

                size_t n = parameters[curr].n;
                double m = parameters[curr].m;
                double m2 = parameters[curr].m2;
                double x = stod(s);
                parameters[curr].m = (x + n*m)/(n + 1);
                parameters[curr].m2 = (x*x + n*m2)/(n + 1);
                parameters[curr].sm2 = m2 - m*m;
                parameters[curr].n++;
                if(method == Method::DIRECT) {
                    parameters[curr].sum += x;
                }
                curr++;
            }

            if (curr != columns) {
                std::cerr << "Row " << row << " mismatch with row 1" << std::endl;
            }
            else {
                //Ok
            }
        }


        if (f.is_open()) {
            f.close();
        }

        //Calculate sigma without error
        if (method == Method::DIRECT) {
            direct(filename);
        }
    }

    void direct(const std::string &filename)
    {
        std::string line, s;
        std::ifstream f;
        size_t curr = 0;
        size_t row = 0;


        f.open(filename.c_str());
        for (int i = 0; i<4; i++) {
            if (getline(f, line)) {
                row++;
            }
            else {
                std::cerr << "Line " << i+1 << " is empty in file: " << filename << std::endl;
            }
        }


        while (getline(f, line)) {
            std::stringstream ss(line);
            curr = 0;
            bool ignore = false;
            row++;

            while (getline(ss, s, ',')) {
                if (ignore) {
                    curr++;
                    continue;
                }
                else if (curr == 25 && s != "PASS_BIN") {
                    curr++;
                    ignore = true;
                    continue;
                }

                if (!is_float(s)) {
                    curr++;
                    continue;
                }

                double m = parameters[curr].m;
                double x = stod(s);
                parameters[curr].dsm2 += (x - m) * (x - m);
                curr++;
            }

        }


        if (f.is_open()) {
            f.close();
        }

    }

    void printCorr(const std::string &filename) {
        std::ofstream f;
        f.open(filename.c_str());

        f << "ParameterA" << ',';
        f << "MA" << ',';
        f << "sA" << ',';
        f << "KS" << ',',
        f << "ParameterB" << ',';
        f << "MB" << ',';
        f << "sB" << ',';
        f << "KS" << ',';
        f << "R" << ',';
        f << std::endl;

#ifdef DEBUG
        std::cout << "ParameterA" << "    ";
        std::cout << "MA" << "    ";
        std::cout << "sA" << "    ";
        std::cout << "ParameterB" << "    ";
        std::cout << "MB" << "    ";
        std::cout << "sB" << "    ";
        std::cout << "R" << "    ";
        std::cout << std::endl;

#endif

        std::sort(correlation.begin(), correlation.end());

        for (const auto &i : correlation) {
            size_t a = i.a;
            size_t b = i.b;
            double r = i.r;
            Parameter &parA = parameters[a];
            Parameter &parB = parameters[b];

            f << parA.name << ',';
            f << parA.m << ',';
            f << std::sqrt(parA.sm2) << ',';
            f << parA.ks << ',';
            f << parB.name << ',';
            f << parB.m << ',';
            f << std::sqrt(parB.sm2) << ',';
            f << parB.ks << ',';
            f << r;
            f << '\n';

#ifdef DEBUG
            std::cout << parameters[a].name << "    ";
            std::cout << parameters[a].m << "    ";
            std::cout << std::sqrt(parameters[a].sm2) << "    ";
            std::cout << parameters[b].name << "    ";
            std::cout << parameters[b].m << "    ";
            std::cout << std::sqrt(parameters[b].sm2) << "    ";
            std::cout << r;
            std::cout << '\n';

#endif
        }

        f << std::flush;

        if (f.is_open()) {
            f.close();
        }

    }

    void printMSA(const std::string &filename) {
        std::ofstream f;
        f.open(filename.c_str());

        f << "Parameter" << ',';
        f << "Unit" << ',';
        f << "LSL" << ',';
        f << "USL" << ',';
        f << "M" << ',';
        f << "s" << ',';
        f << "N" << ',';
        f << "PPKL" << ',';
        f << "PPKH" << ',';
        f << "PPK";
        f << '\n';

#ifdef DEBUG
        std::cout << "Parameter" << "    ";
        std::cout << "Unit" << "    ";
        std::cout << "LSL" << "    ";
        std::cout << "USL" << "    ";
        std::cout << "M" << "    ";
        std::cout << "s" << "    ";
        std::cout << "N" << "    ";
        std::cout << "PPKL" << "    ";
        std::cout << "PPKH" << "    ";
        std::cout << "PPK";
        std::cout << '\n';
#endif

        double sm;
        double m;
        double lsl;
        double usl;
        double ppkl;
        double ppkh;
        double ppk;

        if (method == Method::CUMULATIVE) {
            for (const auto &i : parameters) {
                sm = std::sqrt(i.sm2);
                m = i.m;
                lsl = i.lsl;
                usl = i.usl;
                ppkl = (m - lsl) / 3 / sm;
                ppkh = (usl - m) / 3 / sm;
                ppk = std::min(ppkl, ppkh);
                f << i.name << ',';
                f << i.unit << ',';
                f << i.lsl << ',';
                f << i.usl << ',';
                f << i.m << ',';
                f << sm << ',';
                f << i.n << ',';
                f << ppkl << ',';
                f << ppkh << ',';
                f << ppk;
                f << '\n';
#ifdef DEBUG
                std::cout << i.name << "    ";
                std::cout << i.unit << "    ";
                std::cout << i.lsl << "    ";
                std::cout << i.usl << "    ";
                std::cout << i.m << "    ";
                std::cout << sm << "    ";
                std::cout << i.n << "    ";
                std::cout << ppkl << "    ";
                std::cout << ppkh << "    ";
                std::cout << ppk;
                std::cout << '\n';
#endif
            }

        }
        else if (method == Method::DIRECT) {
            for (const auto &i : parameters) {
                m = i.sum/i.n;
                sm = std::sqrt(i.dsm2/i.n);
                lsl = i.lsl;
                usl = i.usl;
                ppkl = (m - lsl) / 3 / sm;
                ppkh = (usl - m) / 3 / sm;
                ppk = std::min(ppkl, ppkh);
                f << i.name << ',';
                f << i.unit << ',';
                f << i.lsl << ',';
                f << i.usl << ',';
                f << i.m << ',';
                f << sm << ',';
                f << i.n << ',';
                f << ppkl << ',';
                f << ppkh << ',';
                f << ppk;
                f << '\n';

#ifdef DEBUG
                std::cout << i.name << "    ";
                std::cout << i.unit << "    ";
                std::cout << i.lsl << "    ";
                std::cout << i.usl << "    ";
                std::cout << i.m << "    ";
                std::cout << sm << "    ";
                std::cout << i.n << "    ";
                std::cout << ppkl << "    ";
                std::cout << ppkh << "    ";
                std::cout << ppk;
                std::cout << '\n';
#endif
            }
        }

        if (f.is_open()) {
            f.close();
        }
    }

//Setters
public:
    void setMethod(Method method)
    {
        this->method = method;
    }

    Method getMethod()
    {
        return method;
    }

    void clear()
    {
        filedump.clear();
    }


public:
    static bool is_float(const std::string &s)
    {
        constexpr std::array<char, 11> digits = {'.', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

        bool point = false;

        auto it = s.begin();

        if (*it == '-') {
            it++;
        }

        do {
            if (*it == '.') {
                if (point) {
                    return false;
                }
                else {
                    point = true;
                }
            }

            for (auto i : digits) {
                if (i == *it) {
                    break;
                }

                else if (i == '9') {
                    return false;
                }
            }

            it++;
        } while(it != s.end());

        return true;
    }

    static double ks(std::vector<float> &v, double m, double s)
    {
        std::sort(v.begin(), v.end());
        const size_t N = v.size();
        double Dn = 0.0f;

        for (size_t i = 1; i < N; i++) {
            double d = std::abs( phi(v[i], m, s) - static_cast<double>(i)/static_cast<double>(N) );
            if (Dn < d) {
                Dn = d;
            }
        }

        return Dn;
    }

    static double ks(std::vector<double> &v, double m, double s)
    {
        std::sort(v.begin(), v.end());
        const size_t N = v.size();
        double Dn = 0.0f;

        for (size_t i = 1; i < N; i++) {
            double d = std::abs( phi(v[i], m, s) - static_cast<double>(i)/static_cast<double>(N) );
            if (Dn < d) {
                Dn = d;
            }
        }

        return Dn;
    }


//https://www.johndcook.com/blog/cpp_phi/
//The function Φ(x) is the cumulative density function (CDF) of a standard normal (Gaussian) random variable.

    static double phi(double x,  double m, double s)
    {
        // constants
        constexpr double a1 =  0.254829592;
        constexpr double a2 = -0.284496736;
        constexpr double a3 =  1.421413741;
        constexpr double a4 = -1.453152027;
        constexpr double a5 =  1.061405429;
        constexpr double p  =  0.3275911;

        // Save the sign of x
        x = (x - m)/s/M_SQRT2;
        int sign = 1;
        if (x < 0) {
            sign = -1;
        }
        x = abs(x)/M_SQRT2;

        // A&S formula 7.1.26
        double t = 1.0/(1.0 + p*x);
        double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

        return 0.5*(1.0 + sign*y);
    }


protected:
    std::vector<Parameter> parameters;
    std::vector<std::vector<float>> filedump;
    std::vector<Cell> correlation;
    Method method = Method::CUMULATIVE;
};

}
