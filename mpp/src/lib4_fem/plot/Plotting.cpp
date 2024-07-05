#include "Plotting.hpp"

int Plotting::plotVerbose = 0;
bool Plotting::plotEnabled = false;
PlotStatus Plotting::status = FLUSH;
int Plotting::plotPrecision = 4;

std::unique_ptr<Plotting> Plotting::instance = nullptr;

Plotting &Plotting::Instance() {

  if (!instance) {
    instance = std::make_unique<Plotting>();

    Config::Get("PlotVerbose", plotVerbose);
    Config::Get("PlotEnabled", plotEnabled);
    Config::Get("PlotPrecision", plotPrecision);
  }
  return *instance;
}

void Plotting::SetFilePlotting(bool enabled) { plotEnabled = enabled; }

void Plotting::Precision(int p) { plotPrecision = p; }

void Plotting::RemovePlot(const std::string &name) {
  auto plotEntry = plots.find(name);
  if (plotEntry != plots.end()) { plots.erase(plotEntry); }
}

void Plotting::Clear() { plots.clear(); }

void Plotting::AddData(const string &plotName, const string &dataName, const Vector &data) {
  auto &plot = FindPlot(plotName, data.GetMesh());
  plot.AddData(dataName, data);
}

void Plotting::AddDeformation(const string &plotName, const Vector &deformation) {
  auto &plot = FindPlot(plotName, deformation.GetMesh());
  plot.DeformGrid(deformation);
}

VtuPlot &Plotting::FindPlot(const string &name, const Mesh &mesh) {
  auto plotEntry = plots.find(name);
  if (plotEntry == plots.end()) {
    AddPlot(name, mesh);
    plotEntry = plots.find(name);
  }
  return *(plotEntry->second);
}

std::pair<std::string, PlotStatus> mpp::save_plot(const string &filename) {
  return {filename, Plotting::status};
}