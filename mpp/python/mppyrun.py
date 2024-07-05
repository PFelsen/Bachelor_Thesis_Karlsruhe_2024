"""
This file is copied into
the build directory by cmake
in order to be able to execute
test collections
"""
import os
import sys

sys.path.append(os.path.join('..', 'python'))
sys.path.append(os.path.join('..', 'mpp', 'python'))
sys.path.append(os.path.join('..', '..', 'mpp', 'python'))

import argparse

from mppy import Mpp

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Python interface to M++')

    # Run options
    parser.add_argument('-m', '--mute', type=int, default=1,
                        help='mutes output')
    parser.add_argument('-exe', '--executable', type=str, default='M++',
                        help='name of executable')
    parser.add_argument('-n', '--kernels', type=int, default=4,
                        help='number of kernels used by M++')

    # Test options
    parser.add_argument('--mpi_tests', type=int, default=0,
                        help='runs all registered tests in mpi_tests')
    parser.add_argument('--mpp_tests', type=int, default=0,
                        help='runs all registered tests in ctest')
    parser.add_argument('--bench_tests', type=int, default=0,
                        help='runs all registered tests in mpp_bench')
    parser.add_argument('-f', "--fail_state", type=int, default=1,
                        help='forwards the test fail to system')

    args = parser.parse_args()

    project_name = 'M++'
    with open('CMakeCache.txt', 'r') as file:
        for index, line in enumerate(file.readlines()):
            if line.find('CMAKE_PROJECT_NAME') != -1:
                project_name = line[line.rfind('=') + 1:line.rfind('\n')]


    if not bool(args.mute):
        print(f"{os.getcwd()=}")
        print(f"{os.path.isfile(os.getcwd()+'/M++')=}")
        print(f"{os.getcwd()+'/M++'=}")

    if args.mpi_tests and os.path.isfile(os.getcwd()+"/mpi_tests.txt"):
        mpp = Mpp(project_name, executable=args.executable,
                  build_dir=os.getcwd().split("/")[-1],
                  kernels=args.kernels, mute=bool(args.mute), fail=args.fail_state)
    else:
        mpp = Mpp(project_name, executable=args.executable, kernels=args.kernels, mute=bool(args.mute), fail=args.fail_state)
    rc = 0
    if args.mpi_tests:
        rc = mpp.run_mpi_tests(args.kernels)
    elif args.mpp_tests:
        rc = mpp.run_mpp_tests()
    elif args.bench_tests:
        rc = mpp.run_bench_tests(args.kernels)
    sys.exit(0 if rc == 0 else 1)
