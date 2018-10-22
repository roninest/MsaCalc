#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>


namespace msa
{


class Gnuplot
{
public:
    Gnuplot()
    {
#ifdef _WIN32
        gnuplotpipe = _popen("gnuplot", "w");
#else
        gnuplotpipe = popen("gnuplot 1>/dev/null 2>/dev/null", "w");
#endif
        if (!gnuplotpipe) {
            std::cerr << ("Gnuplot not found !") << std::endl;
        }
    }

    virtual ~Gnuplot()
    {
        cleanup();
        close();
    }

    void end()
    {
        this->send("e\n");
    }

    void cmd(const std::string &command)
    {
        fprintf(gnuplotpipe, "%s\n", command.c_str());
        render();
    }

    void send(const std::string &command)
    {
        fprintf(gnuplotpipe, "%s\n", command.c_str());
//        render();
    }

    template<typename X, typename Y>
    void sendRAM(const std::vector<std::vector<std::pair<X, Y>>> &arg)
    {
        std::string tmp;

        size_t maxSize = 0;
        for (const auto &it : arg) {
            if (it->size() > maxSize)
                maxSize = it->size();
        }

        for (size_t i = 0; i < maxSize; i++) {
            for (const auto &it : arg) {
                if (it->size() > i) {
                    tmp += std::to_string(it->at(i).first) + " " + std::to_string(it->at(i).second) + " ";
                }
            }

            this->send(tmp);
            tmp.clear();
        }
    }

//    template<typename X, typename Y>
//    void plotRAM(const std::vector<std::vector<std::pair<X, Y>>> &arg)
//    {
//        sendRAM(arg);
//        end();
//    }

    template<typename X, typename Y>
    void sendRAM(const std::vector<std::pair<X, Y>> &arg)
    {
        for (const auto &it : arg) {
            this->send(std::to_string(it->first) + " " + std::to_string(it->second));
        }
    }

    template<class X>
    void sendRAM(const std::vector<X> &arg)
    {
        for (const auto &it : arg) {
//            if (it < 1000 && it > -1000) {
                this->send(std::to_string(it));
//            }
        }
    }

    void sendRAM(const std::vector<std::string> &arg)
    {
        for (const auto &it : arg) {
            this->send(it);
        }
    }



    template<typename X, typename Y>
    void plot(const std::vector<std::pair<X, Y>> &arg)
    {
        sendRAM(arg);
        end();
    }



    template<class X>
    void plot(const std::vector<X> &arg)
    {
        sendRAM(arg);
        end();
    }


    void plot(const std::vector<std::string> &arg)
    {
        sendRAM(arg);
        end();
    }


    template<typename X, typename Y>
    void plot(const std::vector<std::vector<std::pair<X, Y>>> &arg)
    {
        sendRAM(arg);
        end();
    }

//    template<typename X, typename Y>
//    void plotFile(const std::vector<std::vector<std::pair<X, Y>>> &dataList)
//    {
//        std::stringstream tmp;
//        std::string path;
//        std::string fileName = "gnuplot" + std::to_string(window) + "_" + rand(3)
//            + ".dat";
//        std::vector<std::string> plotDatFiles;
//
//        for (auto it = dataList.begin(); it != dataList.end(); ++it) {
//            fileName = "gnuplot" + std::to_string(window) + "_" + rand(3) + ".dat";
//#ifndef _WIN32
//            path = std::string("/tmp/") + fileName;
//#else
//            path = std::string("./tmp/") + fileName;
//#endif
//            plotDatFiles.push_back(path);
//            this->tmpFiles.push_back(path);
//            vec2dat(it, path);
//        }
//
//        cmd("set grid");
//        cmd(std::string("set term ") + "wxt" + std::string(" ") + std::to_string(window));
//
//        //	list2dat(dataList, path);
//        typedef std::vector<std::string> StringVec;
//        typedef StringVec::const_iterator StringVecCIt;
//
//        for (StringVecCIt it = plotDatFiles.begin(); it != plotDatFiles.end();
//             ++it) {
//            if (it == plotDatFiles.begin()) {
//                tmp << "plot " << " '" << it << "' " << "using " << 1 << ":"
//                    << 2 << " with linespoints pt 7 ps 0.5" << std::endl;
//            }
//            else {
//                tmp << "replot " << " '" << it << "' " << "using " << 1 << ":"
//                    << 2 << " with linespoints pt 7 ps 0.5" << std::endl;
//            }
//
//        }
//
//        std::string command = tmp.str();
//        cmd(command);
//
//        tmp.str(std::string());
//        tmp.clear();
//        this->window++;
//    }

//    template<typename X, typename Y>
//    void plotFile(const std::vector<std::pair<X, Y>> &dataVec,
//                  const std::string &param = "")
//    {
//        std::stringstream tmp;
//        std::string path;
//        std::string fileName = "gnuplot" + std::to_string(window) + "_" + rand(3)
//            + ".dat";
//        std::vector<std::string> plotDatFiles;
//
//        fileName = "gnuplot" + std::to_string(window) + "_" + rand(3) + ".dat";
//#ifndef _WIN32
//        path = std::string("/tmp/") + fileName;
//#else
//        path = std::string("./tmp/") + fileName;
//#endif
//        plotDatFiles.push_back(path);
//        this->tmpFiles.push_back(path);
//        vec2dat(dataVec, path);
//
//        cmd("set grid");
//        cmd(std::string("set term ") + "wxt" + std::string(" ") + std::to_string(window));
//
//        //	list2dat(dataList, path);
//        typedef std::vector<std::string> StringVec;
//        typedef StringVec::const_iterator StringVecCIt;
//
//        for (StringVecCIt it = plotDatFiles.begin(); it != plotDatFiles.end();
//             ++it) {
//            if (it == plotDatFiles.begin()) {
//                tmp << "plot " << " '" << it << "' " << "using " << 1 << ":"
//                    << 2 << " with linespoints pt 7 ps 0.5" << std::endl;
//            }
//            else {
//                tmp << "replot " << " '" << it << "' " << "using " << 1 << ":"
//                    << 2 << " with linespoints pt 7 ps 0.5" << std::endl;
//            }
//
//        }
//
//        std::string command = tmp.str();
//        cmd(command);
//
//        tmp.str(std::string());
//        tmp.clear();
//        this->window++;
//    }

//    void plotDat(const std::string &path, size_t col)
//    {
//        std::stringstream tmp;
//
//        cmd("set grid");
//        cmd(std::string("set term ") + "wxt" + std::string(" ") + std::to_string(window));
//
//        for (size_t die = 1; die <= 2 * col; die += 2) {
//            if (die == 1) {
//                tmp << "plot " << " '" << path << "' " << "using " << die << ":"
//                    << (die + 1) << " with linespoints pt 7 ps 0.5"
//                    << std::endl;
//            }
//            else {
//                tmp << "replot " << " '" << path << "' " << "using " << die
//                    << ":" << (die + 1) << " with linespoints pt 7 ps 0.5"
//                    << std::endl;
//            }
//        }
//
//        std::string command = tmp.str();
//        cmd(command);
//
//        tmp.str(std::string());
//        tmp.clear();
//
//        window++;
//    }

    void run(const std::string &script)
    {
        //TODO:
    }

    void render()
    {
        fflush(gnuplotpipe); // flush needed to start render
    }

    void cleanup()
    {
//        if (!tmpFiles.empty()) {
//            for (auto it = tmpFiles.begin(); it != tmpFiles.end(); ++it) {
//                if (fileExist(it))
//                    remove(it->c_str());
//            }
//        }
    }

    void close()
    {
        cmd("exit");
#ifdef _WIN32
        _pclose(gnuplotpipe);
#else
        pclose(gnuplotpipe);
#endif
    }

protected:
    FILE *gnuplotpipe = nullptr;
//    size_t window = 0;
//    std::vector<std::string> tmpFiles;
//    bool ready = false;
};

} //end of namespace Data
