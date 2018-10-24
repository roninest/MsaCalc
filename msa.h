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
    size_t n = 0u;
    std::string name = "";
    std::string unit = "";
    double lsl = 0.0f;
    double usl = 0.0f;
    double m = 0.0f;
    double m2 = 0.0f;
    double sm2 = 0.0f;
    double sum = 0.0f;
    double dsm2 = 0.0f;
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
//        omp_set_num_threads(4);
        #pragma omp parallel for
        for (size_t i = 0; i < parameters.size(); i++) {
            if (parameters[i].n > 3) {
                Gnuplot gp;
                gp.send("set terminal png");
                std::string cmd = std::string("set output '") +  parameters[i].name + ".png'";
                gp.send(cmd);
//                gp.send("set grid");
                gp.send("plot '-' smooth frequency with dots"); //with dots
                gp.plot(filedump[i]);
            }
            else {
                //No sense to plot empty data
            }
        }
    }

    void corr() {
        long count = 0;

        for (size_t i = 0; i < parameters.size(); i++)  {

            for (size_t j = 0; j < parameters.size(); j++) {
                if (i > j) {
                    if (filedump[i].size() > 10000 && filedump[j].size() > 10000) {
                        double coff = corr(filedump[i], filedump[j], 10000);
                        if (coff < 1.0f && coff > -1.0f) {
                            std::cout << parameters[i].name << ',';
                            std::cout << parameters[j].name << ',';
                            std::cout << coff;
                            std::cout << std::endl;
                        }

                    }

//                    count++;
                }
            }
        }
    }

    double corr(const std::vector<float> &a, const std::vector<float> &b, size_t n) {
        double sumA = 0.0f, sumB = 0.0f, sumAB = 0.0f;
        double sumA2 = 0.0f, sumB2 = 0.0f;
        double sumAmB = 0.0f;

        for (size_t i = 0; i < n; i++) {
            // sum of elements of array X.
            sumA = sumA + a[i];

            // sum of elements of array Y.
            sumB = sumB + b[i];

            sumAmB = sumAmB + (a[i] - b[i]);

            // sum of X[i] * Y[i].
            sumAB = sumAB + a[i] * b[i];

            // sum of square of array elements.
            sumA2 = sumA2 + a[i] * a[i];
            sumB2 = sumB2 + b[i] * b[i];
        }
//        double mA = parameters[a].m;
//        double mB = parameters[b].m;
//        double sA = std::sqrt(parameters[a].sm2);
//        double sB = std::sqrt(parameters[b].sm2);
//
//        double sumA = 0.0f;
//        double sumB = 0.0f;
//
//        const auto vecA = filedump[a];
//        const auto vecB = filedump[b];
//
//        for (size_t i = 0; i < n; i++) {
//            sumA += (vecA[i] - mA);
//            sumB += (vecB[i] - mB);
//        }

//        std::cout << std::endl;

//        double corr = (sumA*sumB) / n/sA/sB;

        return  (n*sumAB - sumA*sumB) / std::sqrt( (n*sumA2 - sumA*sumA ) * (n*sumB2 - sumB*sumB) );
    }

    void load()
    {
        if (filedump.empty()) {
            std::cerr << "Data is empty" << std::endl;
            return;
        }

//        omp_set_num_threads(4);
//        #pragma omp parallel for
        for (size_t i = 0; i < parameters.size(); i++) {
            size_t n = 0u;
            double m = parameters[i].m;
            double m2 = parameters[i].m2;
            double sm2 = parameters[i].sm2;
//            omp_set_num_threads(4);
//            #pragma omp parallel for
            for (const auto &it : filedump[i]) {
                double x = it;
                m = (x + n*m)/(n + 1);
                m2 = (x*x + n*m2)/(n + 1);
                sm2 = m2 - m*m;
                n++;
            }
            parameters[i].m = m;
            parameters[i].m2 = m2;
            parameters[i].sm2 = sm2;
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

    void print(const std::string& filename) {
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

        if (method == Method::CUMULATIVE) {
            for (const auto &i : parameters) {
                sm = std::sqrt(i.sm2);
                m = i.m;
                lsl = i.lsl;
                usl = i.usl;
                double ppkl = (m - lsl) / 3 / sm;
                double ppkh = (usl - m) / 3 / sm;
                double ppk = std::min(ppkl, ppkh);
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
                double ppkl = (m - lsl) / 3 / sm;
                double ppkh = (usl - m) / 3 / sm;
                double ppk = std::min(ppkl, ppkh);
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


protected:
    std::vector<Parameter> parameters;
    std::vector<std::vector<float>>filedump;
    Method method = Method::CUMULATIVE;
};

}
