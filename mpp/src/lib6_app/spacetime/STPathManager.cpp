#include "CyclicPreconditioner.hpp"
#include "STPathManager.hpp"
#include "SpaceTimePreconditioner.hpp"

std::unique_ptr<Preconditioner> createPreconditioner(const std::string &pathChoice,
                                                     STAssemble &assemble) {
  std::unique_ptr<Preconditioner> PC;
  if (pathChoice == "none") {
    std::string pc_name = "PointBlockGaussSeidel";
    Config::Get("Preconditioner", pc_name);
    PC = std::unique_ptr<Preconditioner>(GetPC(pc_name));
  } else if (pathChoice == "adaptive_direct"){
    PC = std::make_unique<STMultiGridPC>(std::make_unique<AdaptiveToZeroDirectPathStrategy>(), assemble);
  } else if (pathChoice == "direct"){
      PC = std::make_unique<STMultiGridPC>(std::make_unique<DirectPathStrategy>(), assemble);
  } else if (pathChoice == "single" || pathChoice == "space_time") {
    PC = std::make_unique<STMultiGridPC>(std::make_unique<SpaceThenTimePathStrategy>(), assemble);
  } else if(pathChoice == "time_space"){
    PC = std::make_unique<STMultiGridPC>(std::make_unique<TimeThenSpacePathStrategy>(), assemble);
  } else if(pathChoice == "time"){
    PC = std::make_unique<STMultiGridPC>(std::make_unique<OnlyTimePathStrategy>(), assemble);
  } else if(pathChoice == "space"){
    PC = std::make_unique<STMultiGridPC>(std::make_unique<OnlySpacePathStrategy>(), assemble);
  }
  string preconditioner;
  Config::Get("Preconditioner", preconditioner);
  if (preconditioner == "CyclicPreconditioner") {
    vector<string> pc_names;
    vector<int> cycle_indices;
    Config::Get("CyclicPCNames", pc_names);
    Config::Get("CyclicPCIndices", cycle_indices);
    std::vector<std::shared_ptr<Preconditioner>> pc;
    pc.reserve(pc_names.size());
    for (const string &name: pc_names) {
      if (name == "MG") {
        pc.push_back(std::move(PC));
      } else {
        pc.push_back(std::shared_ptr<Preconditioner>(GetPC(name)));
      }
    }
    return std::make_unique<CyclicPreconditioner>(pc, cycle_indices);
  } else {
    return std::move(PC);
  }
}







