#ifndef PLOTTINGH
#define PLOTTINGH

#include <string>
#include <vector>
#include <tuple>
#include <memory>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TError.h"

//TODO: make these enum classes, item 10 Effective Modern C++
enum plotOptions{x_label,y_label,bin_edges}; //instead of 0, 1, 2 for the kwargs std::tuple

void _mkdir_recursive(TFile*, const std::string&); //helper function for Plot::save and Plot2D::save

/*
 Plot: A utility class that gives a nice interface for a TH1F object
 */

class Plot
{
private:
    std::string name;
    std::string path;
    std::string write_name;
    std::string title;
    int __nbins;
    double __x_min;
    double __x_max;
    std::string __x_label;
    std::string __y_label;
    std::vector<double> __bin_edges;
    std::shared_ptr<TH1F> __hist;
public:
    Plot();
    ~Plot();
    Plot(const std::string name, std::string title="", int nbins = 40, double x_min = 0, double x_max = 20, std::tuple<std::string, std::string, std::vector<double>>&& kwargs = {});
    
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

/*
 Plot2D: A utility class that gives a nice interface for a TH2F object
 */
class Plot2D
{
private:
    std::string name;
    std::string path;
    std::string write_name;
    std::string title;
    int __xbins;
    double __x_min;
    double __x_max;
    int __ybins;
    double __y_min;
    double __y_max;
    std::string __x_label;
    std::string __y_label;
    std::shared_ptr<TH2F> __hist;
    
public:
    Plot2D();
    ~Plot2D();
    Plot2D(std::string name, std::string title="", int xbins = 40, double x_min = 0, double x_max = 20, int ybins = 40, double y_min = 0, double y_max = 20, std::tuple<std::string, std::string>&& kwargs = {});
    
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

/*
 PlotGroup: A utility class that gives a nice interface for a vector of Plot objects
 */
class PlotGroup //1D
{
private:
    std::vector<Plot> __hists;
public:
    PlotGroup();
    ~PlotGroup();
    PlotGroup(std::vector<Plot>&&);
    
    PlotGroup(const PlotGroup&);
    PlotGroup& operator=(const PlotGroup&);
    
    void add(const PlotGroup&);
    void fill(double value, double weight=1.0);
    void save(TFile* out_file = nullptr);
    std::vector<Plot> hists();
};

/*
 PlotGroup2D: A utility class that gives a nice interface for a vector of Plot2D objects
 */
class PlotGroup2D //2D
{
private:
    std::vector<Plot2D> __hists;
public:
    PlotGroup2D();
    ~PlotGroup2D();
    PlotGroup2D(std::vector<Plot2D>&&);
    
    PlotGroup2D(const PlotGroup2D&);
    PlotGroup2D& operator=(const PlotGroup2D&);
    
    void add(const PlotGroup2D&);
    void fill(double x, double y, double weight=1.0);
    void save(TFile* out_file = nullptr);
    std::vector<Plot2D> hists();
};

#endif
