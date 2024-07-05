import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--out', '-o', required=True)
parser.add_argument('opt_list', nargs=argparse.REMAINDER)

args = parser.parse_args()

generated = """
#include "Logging.hpp"

void MasterLogging::PrintBuildInfo() {
    if (!onMaster)
        return;

    PrintInfo("Build", 1,
"""
num_opts = len(args.opt_list)
for i  in range(0, num_opts - 2, 2):
    generated += f"""        PrintInfoEntry(\"{args.opt_list[i]}\", std::string(\"{args.opt_list[i + 1]}\")),
"""

generated += f"""        PrintInfoEntry(\"{args.opt_list[num_opts - 2]}\", std::string(\"{args.opt_list[num_opts - 1]}\")));
}}
"""

os.makedirs(os.path.dirname(args.out), exist_ok=True)
with open(args.out, 'w') as f:
    f.write(generated)
