#ifndef SYNCHRONIZATION_HPP
#define SYNCHRONIZATION_HPP

#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Point.hpp"

class IPlotBackend;
class Mesh;
class ProcSet;
class VtuPlot;

using BackendPlotData = std::tuple<std::vector<Point>, std::vector<double>, uint32_t>;

class DataSynchronization {
protected:
  std::string domainName;
  IPlotBackend &backend;

  explicit DataSynchronization(IPlotBackend &backend) : backend(backend) {}
public:
  virtual ~DataSynchronization() = default;

  const IPlotBackend &GetBackend() const { return backend; }

  /**
   * @brief Adds commonized cell data to the plotting system.
   *
   * @param data the cell data to add. A Point represents the center of a cell.
   * The data of that Point starts at index of Point *  uint32_t (data
   * size)
   */
  virtual void AddCellData(BackendPlotData &data) = 0;

  /**
   * @brief Adds commonized point data to the plotting system.
   *
   * @param data the point data to add. The data of a Point starts at index of
   * Point * uint32_t (data size)
   */
  virtual void AddPointData(BackendPlotData &data) = 0;

  /**
   * @brief Get the current domain name
   *
   * @return std::string& the domain name
   */
  std::string &GetDomainName() { return domainName; }

  /**
   * @brief Set the domain name which will be used for the following
   * transactions
   *
   * @param domainName the domain name
   */
  void SetDomainName(const std::string &domainName) { this->domainName = domainName; }
};

class MeshSynchronization {
private:
  friend VtuPlot;
protected:
  IPlotBackend &backend;

  explicit MeshSynchronization(IPlotBackend &backend) : backend(backend) {}

  /**
   * @brief Adds cell related information to the @link IPlotBackend @endlink.
   *
   * @param mesh the mesh to take the cell infos from
   */
  virtual void addCells(const Mesh &mesh) = 0;
public:
  virtual ~MeshSynchronization() = default;

  /**
   * @brief Adds commonized deformations of the mesh to the plotting system. If
   * there may are multiple displacements for the same point on different
   * processes, one must add it to 'possibleDuplicates'
   *
   * @param deformation the deformation to add. The first Point
   * represents the Point to displace, the second Point describes the
   * displacement
   * @param possibleDuplicates contains points which may also exist on other
   * processes with different displacement. Those processes must be contained in
   * the @link ProcSet @endlink
   */
  virtual void AddDeformation(const std::vector<std::pair<Point, Point>> &deformation,
                              const std::vector<std::pair<Point, ProcSet>> &possibleDuplicates) = 0;
};

#endif
