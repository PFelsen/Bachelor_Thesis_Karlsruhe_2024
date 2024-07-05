import os
import time
import matplotlib.pyplot as plt
import json

from .utilities.directory_manager import DirectoryManager, remove_file, clean_directory
from .utilities.log_utilities import LogParser
from .utilities.subprocess_utilities import SubprocessManager
from .utilities.test_utilities import print_test_results
from .utilities.vtu_utilities import VtuPlot


class Mpp:
    def __init__(self, project_name='mpp', executable='M++', kernels=4, mute=False,
                 fail=False, build_dir='build', cmake_args=None):

        self.project_name = project_name
        self.executable = executable
        self.kernels = kernels
        self.mute = mute

        self.cmake_args = ['cmake', '..']
        if type(cmake_args) is list:
            self.cmake_args = self.cmake_args + cmake_args
        else:
            if cmake_args is not None:
                self.cmake_args.append(cmake_args)

        self.dm = DirectoryManager(self.project_name, build_dir=build_dir, mute=self.mute)
        self.sm = SubprocessManager(self.mute)
        self.lp = LogParser(self.dm.PROJECT_LOG_DIR)

        self.cwd = self.dm.PROJECT_BUILD_DIR
        self.make_args = ['make', '-j', str(self.kernels)]

        self.data = []

    def create_working_dir(self):
        if not os.path.isdir(self.dm.PROJECT_BUILD_DIR):
            os.mkdir(self.dm.PROJECT_BUILD_DIR)

    def vtu_plot(self, **kwargs):
        p = VtuPlot(**kwargs)
        p.set_wd(self.dm.PROJECT_VTU_DATA_DIR)
        plt.axis('off')
        return p

    def run_cmake(self):
        if not self.mute:
            print('\n================ running cmake ================\n')
        return self.run_subprocess(self.cmake_args, self.dm.PROJECT_BUILD_DIR, self.mute)

    def run_make(self):
        if not self.mute:
            print('\n================ running  make ================\n')
        return self.run_subprocess(self.make_args, self.dm.PROJECT_BUILD_DIR, self.mute)

    def kill(self):
        if not self.mute:
            print('\n================ kill every mpp ================\n')
        self.run_subprocess(['killall', self.executable], self.cwd, self.mute)
        return 0

    def build(self):
        self.create_working_dir()
        rc_cmake = self.run_cmake()
        rc_make = self.run_make()
        if rc_cmake == 0 and rc_make == 0:
            return 0
        else:
            return 1

    def run(self, kernels=None, config=None, kwargs=None):
        if not self.mute:
            print('\n================ running  mpp ================\n')
        run_parameters = self.chain_run_parameters(kernels, config, kwargs)
        return self.run_subprocess(run_parameters, self.cwd, self.mute)

    def chain_run_parameters(self, kernels, config, kwargs, executable=None):
        executable = self.executable if executable is None else executable
        kernels = self.kernels if kernels is None else kernels
        kwargs = {} if kwargs is None else kwargs
        if type(kernels) is dict:
            kernels = {str(k).zfill(2): v for k, v in kernels.items()}
            with open(os.path.join(self.cwd, "hosts"), "w") as f:
                for cid in reversed(sorted(kernels.keys())):
                    f.write("ma-pde{} slots={}\n".format(cid, kernels[cid]))
            kernel_count = str(sum(kernels.values()))
            run_parameters = ['mpirun', '--use-hwthread-cpus', '-x', 'LD_LIBRARY_PATH',
                              '-hostfile', 'hosts', '-np', kernel_count,
                              executable]
        else:
            run_parameters = ['mpirun', '--use-hwthread-cpus',
                              '-np', str(kernels), executable]
        if config:
            #if config.endswith(".json"):
            #    config = os.path.join(self.dm.MPP_JSON_DIR, config)
            run_parameters.append(config)
        for key, arg in kwargs.items():
            kwarg = key + '=' + str(arg)
            run_parameters.append(kwarg)
        return run_parameters

    def run_subprocess(self, args, cwd, mute=False):
        return self.sm.run_subprocess(args, cwd, mute)

    def reset_data(self):
        self.lp.reset_data()
        self.data = []

    def parse_log(self, log_file='logfile.log', additional_entries=None):
        self.data = self.lp.parse_log(log_file, additional_entries)
        return self.data
    
    def parse_json(self, json_file='logfile.json', directory=None):
        directory = self.dm.MPP_JSON_DIR if directory is None else directory
        if json_file == "all":
            for json_file in sorted(os.listdir(directory)):
                with open(os.path.join(directory, json_file)) as file:
                    self.data.append(json.load(file))
        else:
            with open(os.path.join(directory, json_file)) as file:
                self.data.append(json.load(file))
        return self.data

    def clean_vtk(self):
        clean_directory(self.dm.PROJECT_VTK_DATA_DIR, True)

    def clean_vtu(self):
        clean_directory(self.dm.PROJECT_VTU_DATA_DIR, True)

    def clean_log(self):
        clean_directory(self.dm.PROJECT_LOG_DIR, False)

    def clean_json(self):
        clean_directory(self.dm.PROJECT_JSON_DIR, False)

    def clean_python_plots(self):
        clean_directory(self.dm.PROJECT_PY_DATA_DIR, False)

    def clean_data(self):
        self.clean_vtk()
        self.clean_vtu()
        self.clean_log()
        self.clean_json()
        self.clean_python_plots()

    def clean_build(self):
        clean_directory(self.dm.PROJECT_BUILD_DIR, True)

    def clean_cmake_cache(self):
        remove_file(self.dm.cmake_cache)

    def run_mpi_tests(self, kernels):
        print('MPI Tests in project {}'.format(self.dm.PROJECT_BUILD_DIR))
        rc_total = 0

        with open(self.dm.mpi_tests_file_path, 'r') as file:
            lines = file.readlines()
            total_num_of_tests = len(lines)
            for index, line in enumerate(lines):
                # remove /, \n
                absolute_path = line[1:-1]
                self.cwd = "/" + absolute_path[:absolute_path.rfind('/') + 1]
                self.executable = absolute_path[absolute_path.rfind('/') + 1:]

                print(f'        Start {index + 1:>3}: {self.executable}')
                start = time.time()
                rc = self.run(kernels)
                log = self.sm.latest_log
                duration = time.time() - start
                print_test_results(duration, index, self.executable,
                                   rc, log, total_num_of_tests)
                rc_total += rc
        return rc_total

    def run_bench_tests(self, kernels):
        print('MPP Benchmarks in project {}'.format(self.dm.PROJECT_BUILD_DIR))
        rc_total = 0

        with open(self.dm.mpp_bench_file_path, 'r') as file:
            lines = file.readlines()
            total_num_of_tests = len(lines)
            for index, line in enumerate(lines):
                # remove /, \n and cwd
                self.executable = (line[:-1].replace(self.cwd, ''))[1:]
                name = self.executable[self.executable.rfind('/') + 1:]

                bench_kwargs = {
                    "--benchmark_out": name,
                    "--benchmark_out_format": "json"
                }

                print(f'        Start {index + 1:>3}: {name}')
                start = time.time()
                rc = self.run(1, None, bench_kwargs)
                log = self.sm.latest_log
                duration = time.time() - start
                print_test_results(duration, index, name,
                                   rc, log, total_num_of_tests)
                rc_total += rc
        return rc_total

    def run_mpp_tests(self):
        pass


mpp = Mpp(mute=True)


def run_mpi_tests(kernels):
    mpp.run_mpi_tests(kernels)


def run_bench_tests(kernels):
    mpp.run_bench_tests(kernels)


def run_mpp_tests():
    mpp.run_mpp_tests()

