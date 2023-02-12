#include "plotting.h"

#include <iostream>

#include "TH1F.h"

/*
function to recursively create a directory structure in the ROOT file
for the different TH1F's, in practice based on systematics and cuts
*/
void _mkdir_recursive(TFile* out_file, const std::string& full_path)
{
    size_t pos = 0;
    std::string path;
    while ((pos = full_path.find('/',pos)) != std::string::npos)
    {
        path = full_path.substr(0,pos++);
        if (!(out_file->GetDirectory(path.c_str())))
        {
            out_file->mkdir(path.c_str());
        }
    }
    if (!(out_file->GetDirectory(full_path.c_str())))
    {
        out_file->mkdir(full_path.c_str());
    }
}

//                ******************
//                *      Plot      *
//                ******************

/*
Constuctor to create a Plot object that wraps a TH1F, requires the name, title,
number of bins, and the lower and upper bounds x-bounds (x_min and x_max) for
the TH1F, as well three optional arguments for the x-axis title, y-axis title,
and an array of bin edges. The path is everything before the last '/' and the write_name
is everything after the last '/'

e.g.

Plot test("name","title",10,0,10);
Plot test("name","title",10,0,10,std::make_tuple("mass (GeV)","events",std::vector<double>()));
Plot test("name","title",10,0,10,std::make_tuple("mass (GeV)","events",std::vector<double>({1.1,2.2,3.3})));
*/

Plot::Plot(const std::string name, std::string title, int nbins, double x_min, double x_max, std::tuple<std::string, std::string, std::vector<double>>&& kwargs)
:
    name{name},
    path{"/"},
    write_name{name},
    title{title},
    __nbins{nbins},
    __x_min{x_min},
    __x_max{x_max},
    __x_label{std::get<x_label>(kwargs)},
    __y_label{std::get<y_label>(kwargs)},
    __bin_edges{std::get<bin_edges>(kwargs)}
{
    TH1::AddDirectory(kFALSE);
//    __hist->SetDirectory(0);
    size_t temp;
    if ((temp = name.find_last_of('/')) != std::string::npos)
    {
        this->path = name.substr(0,temp);
        this->write_name = name.substr(temp+1);
    }
    if (title.empty())
    {
        title = name;
    }
    if (!(__bin_edges.empty()))
    {
        __hist.reset(new TH1F(name.c_str(), title.c_str(), __bin_edges.size()-1, &__bin_edges[0] ));
    }
    else
    {
        __hist.reset(new TH1F(name.c_str(), title.c_str(),__nbins, __x_min, __x_max));
    }
    __hist->GetXaxis()->SetTitle(__x_label.c_str());
    __hist->GetYaxis()->SetTitle(__y_label.c_str());
}

Plot::~Plot() = default;

Plot::Plot(const Plot& other)
{
    __hist = other.__hist;
    name = other.name;
    path = other.path;
    title = other.title;
    __nbins = other.__nbins;
    __x_min = other.__x_min;
    __x_max = other.__x_max;
    __x_label = other.__x_label;
    __y_label = other.__y_label;
    __bin_edges = other.__bin_edges;
    
}

Plot::Plot() = default;

Plot& Plot::operator=(const Plot& other)
{
    __hist = other.__hist;
    name = other.name;
    path = other.path;
    title = other.title;
    __nbins = other.__nbins;
    __x_min = other.__x_min;
    __x_max = other.__x_max;
    __x_label = other.__x_label;
    __y_label = other.__y_label;
    __bin_edges = other.__bin_edges;
    return *this;
}

void Plot::add(const Plot& plot)
{
    __hist->Add(&(*plot.__hist));
}

void Plot::add(const TH1F* plot)
{
    __hist->Add(plot);
}

void Plot::fill(double value, double weight)
{
    __hist->Fill(value,weight);
}

/*
 Utility function to save Plot object to an output ROOT file at
 the location specified by the path attribute
 */
void Plot::save(TFile* out_file)
{
    if (!out_file)
    {
        std::cerr << "Error in Plot::save, need an output file to save to!\n";
        return;
    }
    if (path != '/')
    {
        _mkdir_recursive(out_file,path);
    }
    out_file->cd(path.c_str());
    __hist->Write(write_name.c_str());
}

void Plot::draw()
{
    __hist->Draw();
}

void Plot::del()
{
    __hist.reset();
}

TH1F Plot::hist()
{
    return *__hist;
}

//                ********************
//                *      Plot2D      *
//                ********************

/*
Constuctor to create a Plot2D object that wraps a TH2F, requires the name, title,
number of bins, and the lower and upper bounds x-bounds and y_bounds
(x_min, y_min and x_max, y_max) for the TH2F, as well
two optional arguments for the x-axis title and the y-axis title,
The path is everything before the last '/' and the write_name
name is everything after the last '/'

e.g.

Plot2D test("name","title",10,0,10,10,0,10);
Plot2D test1("name","title",10,0,10,10,0,10,std::make_tuple("name","title"));
*/

Plot2D::Plot2D(std::string name, std::string title, int xbins, double x_min, double x_max, int ybins, double y_min, double y_max, std::tuple<std::string, std::string>&& kwargs)
:
    name{name},
    path{"/"},
    title{title},
    __xbins{xbins},
    __x_min{x_min},
    __x_max{x_max},
    __ybins{ybins},
    __y_min{y_min},
    __y_max{y_max},
    __x_label{std::get<x_label>(kwargs)},
    __y_label{std::get<y_label>(kwargs)}
{
    TH1::AddDirectory(kFALSE);
    size_t temp;
    if ((temp = name.find_last_of('/')) != std::string::npos)
    {
        this->path = name.substr(0,temp);
        this->write_name = name.substr(temp+1);
    }
    if (title.empty())
    {
        title = name;
    }
    __hist.reset(new TH2F(name.c_str(), title.c_str(), __xbins, __x_min, __x_max, __ybins, __y_min, __y_max));

    __hist->GetXaxis()->SetTitle(__x_label.c_str());
    __hist->GetYaxis()->SetTitle(__y_label.c_str());
}

Plot2D::~Plot2D() = default;

Plot2D::Plot2D(const Plot2D& other)
{
    __hist = other.__hist;
    name = other.name;
    path = other.path;
    title = other.title;
    __xbins = other.__xbins;
    __x_min = other.__x_min;
    __x_max = other.__x_max;
    __ybins = other.__ybins;
    __y_min = other.__y_min;
    __y_max = other.__y_max;
    __x_label = other.__x_label;
    __y_label = other.__y_label;
}

Plot2D::Plot2D() = default;

Plot2D& Plot2D::operator=(const Plot2D& other)
{
    __hist = other.__hist;
    name = other.name;
    path = other.path;
    title = other.title;
    __xbins = other.__xbins;
    __x_min = other.__x_min;
    __x_max = other.__x_max;
    __ybins = other.__ybins;
    __y_min = other.__y_min;
    __y_max = other.__y_max;
    __x_label = other.__x_label;
    __y_label = other.__y_label;
    return *this;
}

void Plot2D::add(const Plot2D& plot)
{
    __hist->Add(&(*plot.__hist));
}

void Plot2D::add(const TH2F* plot)
{
    __hist->Add(plot);
}

void Plot2D::fill(double x, double y, double weight)
{
    __hist->Fill(x,y,weight);
}

void Plot2D::save(TFile* out_file)
{
    if (!out_file)
    {
        std::cerr << "Error in Plot2D::save, need an output file to save to!\n";
        return;
    }
    if (path != '/')
    {
        _mkdir_recursive(out_file,path);
    }
    out_file->cd(path.c_str());
    __hist->Write(write_name.c_str());
}

void Plot2D::draw()
{
    __hist->Draw();
}

void Plot2D::del()
{
    __hist.reset();
}

TH2F Plot2D::hist()
{
    return *__hist;
}

//                ***********************
//                *      PlotGroup      *
//                ***********************

/*
 Constructor to create a group of Plot objects
 */

PlotGroup::PlotGroup(std::vector<Plot>&& hists)
{
    __hists = hists;
}

PlotGroup::PlotGroup() = default;
PlotGroup::~PlotGroup() = default;

PlotGroup::PlotGroup(const PlotGroup& plot_group)
{
    __hists = plot_group.__hists;
}

PlotGroup& PlotGroup::operator=(const PlotGroup& plot_group)
{
    __hists = plot_group.__hists;
    return *this;
}

/*
 Performs TH1::Add on all the Plots in the Plotgroup,
 given that __hists and plots have the same size
 */

void PlotGroup::add(const PlotGroup& plots)
{
    R__ASSERT(__hists.size() == plots.__hists.size());
    
    for (size_t i = 0; i < __hists.size(); ++i)
    {
        __hists[i].add(plots.__hists[i]);
    }
}

void PlotGroup::fill(double value, double weight)
{
    for (auto& hist: __hists)
    {
        hist.fill(value,weight);
    }
}

void PlotGroup::save(TFile* out_file)
{
    for (auto& hist: __hists)
    {
        hist.save(out_file);
    }
}

std::vector<Plot> PlotGroup::hists()
{
    return __hists;
}

//                *************************
//                *      PlotGroup2D      *
//                *************************

PlotGroup2D::PlotGroup2D() = default;
PlotGroup2D::~PlotGroup2D() = default;

PlotGroup2D::PlotGroup2D(const PlotGroup2D& plot_group)
{
    __hists = plot_group.__hists;
}

PlotGroup2D& PlotGroup2D::operator=(const PlotGroup2D& plot_group)
{
    __hists = plot_group.__hists;
    return *this;
}

PlotGroup2D::PlotGroup2D(std::vector<Plot2D>&& hists)
{
    __hists = hists;
}

void PlotGroup2D::add(const PlotGroup2D& plots)
{
    R__ASSERT(__hists.size() == plots.__hists.size());
    
    for (size_t i = 0; i < __hists.size(); ++i)
    {
        __hists[i].add(plots.__hists[i]);
    }
}

void PlotGroup2D::fill(double x, double y, double weight)
{
    for (auto& hist: __hists)
    {
        hist.fill(x,y,weight);
    }
}

void PlotGroup2D::save(TFile* out_file)
{
    for (auto& hist: __hists)
    {
        hist.save(out_file);
    }
}

std::vector<Plot2D> PlotGroup2D::hists()
{
    return __hists;
}
