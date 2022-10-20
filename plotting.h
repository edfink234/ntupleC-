#ifndef PLOTTINGH
#define PLOTTINGH

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TError.h"
using std::get;
using std::cout;

using std::string;
using std::vector;
using std::tuple;

//TODO: make this enum class, item 10 Effective Modern C++
enum plotOptions{x_label,y_label,bin_edges};

void _mkdir_recursive(TFile*, const string&);

class Plot
{
private:
    string name;
    string path;
    string write_name;
    string title;
    int __nbins;
    int __x_min;
    int __x_max;
    string __x_label;
    string __y_label;
    vector<double> __bin_edges;
    std::unique_ptr<TH1F> __hist;
    
public:
//    string name;
    Plot();
    ~Plot();
    Plot(string name, string title="", int nbins = 40, int x_min = 0, int x_max = 20, tuple<string, string, vector<double>>&& kwargs = {});
    
    Plot(const Plot&);
    Plot& operator=(const Plot&);
    
    void add(const Plot&);
    void add(const TH1F*);
    void fill(double value, double weight = 1.0);
    void save(TFile* out_file = nullptr);
    void draw();
//    TH1F unfold(...);
    void del();
    TH1F hist();
};

class Plot2D
{
private:
    string name;
    string path;
    string write_name;
    string title;
    int __xbins;
    int __x_min;
    int __x_max;
    int __ybins;
    int __y_min;
    int __y_max;
    string __x_label;
    string __y_label;
    std::unique_ptr<TH2F> __hist;
    
public:
//    string name;
    Plot2D();
    ~Plot2D();
    Plot2D(string name, string title="", int xbins = 40, int x_min = 0, int x_max = 20, int ybins = 40, int y_min = 0, int y_max = 20, tuple<string, string, vector<double>>&& kwargs = {});
    
    Plot2D(const Plot2D&);
    Plot2D& operator=(const Plot2D&);
    
    void add(const Plot2D&);
    void add(const TH2F*);
    void fill(double x, double y, double weight = 1.0);
    void save(TFile* out_file = nullptr);
    void draw();
//    TH1F unfold(...);
    void del();
    TH2F hist();
};

class PlotGroup //1D
{
private:
    vector<Plot> __hists;
public:
    PlotGroup();
    ~PlotGroup();
    PlotGroup(vector<Plot>&&);
    
    PlotGroup(const PlotGroup&);
    PlotGroup& operator=(const PlotGroup&);
    
    void add(const PlotGroup&);
    void fill(double value, double weight=1.0);
    void save(TFile* out_file = nullptr);
    vector<Plot> hists();
};

class PlotGroup2D //2D
{
private:
    vector<Plot2D> __hists;
public:
    PlotGroup2D(vector<Plot2D>&&);
    void add(const PlotGroup2D&);
    void fill(double x, double y, double weight=1.0);
    void save(TFile* out_file = nullptr);
    vector<Plot2D> hists();
};



#endif
