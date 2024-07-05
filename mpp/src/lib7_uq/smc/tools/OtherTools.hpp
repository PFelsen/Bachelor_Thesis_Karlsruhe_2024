#ifndef MLUQ_OTHER_TOOLS_H
#define MLUQ_OTHER_TOOLS_H

#include <vector>

void normalizeWeights(std::vector<double> &weights);

std::vector<double> mean_in_column(std::vector<std::vector<double>> samples);

std::vector<double> variance_in_column(std::vector<std::vector<double>> samples,
                                       std::vector<double> sample_mean);

#endif // MLUQ_OTHER_TOOLS_H