#include "plotting.h"


void _mkdir_recursive(TFile* out_file, const string& full_path)
{
    size_t pos=0;
    string path;
    while ((pos=full_path.find('/',pos))!=string::npos)
    {
        path = full_path.substr(0,pos++);
//        cout << path << '\n';
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

Plot::Plot(string name, string title, int nbins, int x_min, int x_max, tuple<string, string, vector<double>>&& kwargs)
:
    name{name},
    path{"/"},
    write_name{name},
    title{title},
    __nbins{nbins},
    __x_min{x_min},
    __x_max{x_max},
    __x_label{get<x_label>(kwargs)},
    __y_label{get<y_label>(kwargs)},
    __bin_edges{get<bin_edges>(kwargs)}
{
    size_t temp;
    if ((temp = name.find_last_of('/')) != std::string::npos)
    {
        this->path = name.substr(0,temp);
        this->write_name = name.substr(temp+1);
    }
    if (title.empty())
    {
        title=name;
    }
    if (!(__bin_edges.empty()))
    {
//    https://stackoverflow.com/questions/3093451/is-it-safe-to-pass-a-vector-as-an-array
        __hist.reset(new TH1F(name.c_str(), title.c_str(), __bin_edges.size()-1, &__bin_edges[0] ));
    }
    else
    {
        __hist.reset(new TH1F(name.c_str(), title.c_str(),__nbins, __x_min, __x_max));
    }
    __hist->GetXaxis()->SetTitle(__x_label.c_str());
    __hist->GetYaxis()->SetTitle(__y_label.c_str());
//    cout << name << '\n';
}

Plot::~Plot()
{
    __hist.reset(nullptr);
    
} //?

//https://stackoverflow.com/questions/16030081/copy-constructor-for-a-class-with-unique-ptr
Plot::Plot(const Plot& other) : __hist(new TH1F(*other.__hist))
{
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
    name = other.name;
    path = other.path;
    title = other.title;
    __nbins = other.__nbins;
    __x_min = other.__x_min;
    __x_max = other.__x_max;
    __x_label = other.__x_label;
    __y_label = other.__y_label;
    __bin_edges = other.__bin_edges;
    __hist.reset(new TH1F(*other.__hist));
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
void Plot::save(TFile* out_file)
{
    if (!out_file)
    {
        out_file = __hist->GetDirectory()->GetFile();
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
    __hist.reset(nullptr);
}

TH1F Plot::hist()
{
    return *__hist;
}

//                ********************
//                *      Plot2D      *
//                ********************


Plot2D::Plot2D(string name, string title, int xbins, int x_min, int x_max, int ybins, int y_min, int y_max, tuple<string, string, vector<double>>&& kwargs)
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
    __x_label{get<x_label>(kwargs)},
    __y_label{get<y_label>(kwargs)}
{
    size_t temp;
    if ((temp = name.find_last_of('/')) != std::string::npos)
    {
        this->path = name.substr(0,temp);
        this->write_name = name.substr(temp+1);
    }
    if (title.empty())
    {
        title=name;
    }


    __hist.reset(new TH2F(name.c_str(), title.c_str(),__xbins, __x_min, __x_max, __ybins, __y_min, __y_max));

    __hist->GetXaxis()->SetTitle(__x_label.c_str());
    __hist->GetYaxis()->SetTitle(__y_label.c_str());
//    cout << name << '\n';
}

Plot2D::~Plot2D(){__hist.reset(nullptr);} //?

//https://stackoverflow.com/questions/16030081/copy-constructor-for-a-class-with-unique-ptr
Plot2D::Plot2D(const Plot2D& other) : __hist(new TH2F(*other.__hist))
{
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
    __hist.reset(new TH2F(*other.__hist));
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
        out_file = __hist->GetDirectory()->GetFile();
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
    __hist.reset(nullptr);
}

TH2F Plot2D::hist()
{
    return *__hist;
}

//                ***********************
//                *      PlotGroup      *
//                ***********************

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


PlotGroup::PlotGroup(vector<Plot>&& hists)
{
    __hists = hists;
}

void PlotGroup::add(const PlotGroup& plots)
{
    R__ASSERT(__hists.size() == plots.__hists.size());
    
    for (size_t i=0; i < __hists.size(); ++i)
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

vector<Plot> PlotGroup::hists()
{
    return __hists;
}

//                *************************
//                *      PlotGroup2D      *
//                *************************

PlotGroup2D::PlotGroup2D(vector<Plot2D>&& hists)
{
    __hists = hists;
}

void PlotGroup2D::add(const PlotGroup2D& plots)
{
    R__ASSERT(__hists.size() == plots.__hists.size());
    
    for (size_t i=0; i < __hists.size(); ++i)
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

vector<Plot2D> PlotGroup2D::hists()
{
    return __hists;
}
