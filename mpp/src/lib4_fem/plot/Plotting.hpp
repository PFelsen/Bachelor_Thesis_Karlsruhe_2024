#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <unordered_map>
#include <utility>

#include "VtuPlot.hpp"

class Plotting {
  static std::unique_ptr<Plotting> instance;
  std::unordered_map<std::string, std::unique_ptr<VtuPlot>> plots{};
protected:
  static string plotPath;
  static bool plotEnabled;
  static int plotVerbose;
  static int plotPrecision;
public:
  static Plotting &Instance();

  static void SetFilePlotting(bool enabled);

  static void Precision(int p);
  static PlotStatus status;

  inline void AddPlot(const string &name, const PlotConfig &plotConfig = PlotConfig()) {
    plots.try_emplace(name, std::make_unique<VtuPlot>(name, plotConfig));
  }

  inline void AddPlot(const std::string &name, const Mesh &mesh,
                      const PlotConfig &plotConfig = PlotConfig()) {
    plots.try_emplace(name, std::make_unique<VtuPlot>(mesh, name, plotConfig));
  }

  inline VtuPlot &GetPlot(const string &name, const PlotConfig &plotConfig = PlotConfig()) {
    auto plotEntry = plots.find(name);
    if (plotEntry == plots.end()) {
      AddPlot(name, plotConfig);
      return *(plots.find(name))->second << name;
    }
    return *(plotEntry->second);
  }

  void RemovePlot(const std::string &plotName);
  void Clear();

  inline void AddPlot(const std::string &plotName, const Vector &v,
                      const PlotConfig &plotConfig = PlotConfig()) {
    AddPlot(plotName, v.GetMesh(), plotConfig);
  };

  inline void AddFullVertData(const std::string &plotName, const std::string &dataName,
                              const Vector &data) {
    AddData(plotName, dataName, data);
  }

  void AddData(const std::string &plotName, const std::string &dataName, const Vector &data);

  void AddDeformation(const std::string &plotName, const Vector &deformation);

  inline void SavePlot(const std::string &name, const std::string &filename) {
    const auto plotEntry = plots.find(name);
    if (plotEntry != plots.end()) { plotEntry->second->PlotFile(filename); }
  };

  inline void SavePlot(const std::string &name) { SavePlot(name, name); };

  VtuPlot &FindPlot(const string &name, const Mesh &mesh);
};

inline VtuPlot &GetPlot(const string &name, const PlotConfig &plotConfig = PlotConfig()) {
  return Plotting::Instance().GetPlot(name, plotConfig);
}

namespace mpp {
inline void PlotData(const string &name, const string &fileName, const Vector &data,
                     PlotConfig plotConfig = PlotConfig()) {
  VtuPlot vtuPlot(fileName, plotConfig);
  vtuPlot.AddData(name, data);
  vtuPlot.PlotFile();
}

inline VtuPlot &plot(const string &name, const PlotConfig &plotConfig = PlotConfig()) {
  return GetPlot(name, plotConfig);
}

inline VtuPlot &deformed_plot(const string &name, const Vector &deformation,
                              const PlotConfig &plotConfig = {
                                  .plotBackendCreator = std::make_unique<DisplacementBackend>,
                              }) {
  auto &p = GetPlot(name, plotConfig);
  p.DeformGrid(deformation);
  return p;
}

std::pair<std::string, PlotStatus> save_plot(const string &filename);
constexpr PlotStatus endp = PlotStatus::CLEAR;

inline void plot_mesh(const Mesh &mesh, const PlotConfig &plotConfig = PlotConfig()) {
  auto name = mesh.Name();
  if (name.empty()) name = "Mesh";
  Plotting::Instance().AddPlot(name, mesh, plotConfig);
  Plotting::Instance().SavePlot(name);
  Plotting::Instance().RemovePlot(name);
}

static string intAsString(int i) {
  char buffer[256];
  sprintf(buffer, "%04d", i);
  string str(buffer);
  return string("-") + str;
}
} // namespace mpp

#endif // PLOTTING_HPP
