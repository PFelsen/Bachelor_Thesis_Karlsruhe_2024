from os import unlink, listdir

from os.path import abspath, join
from os.path import dirname, realpath
from os.path import isfile, isdir

from shutil import rmtree


class DirectoryManager:
    def __init__(self, project_name='mpp', build_dir='build', mute=True):
        self.PY_FILE_DIR = dirname(realpath(__file__))
        self.MPP_ROOT_DIR = abspath(join(self.PY_FILE_DIR, '../../../'))

        if project_name in ('mpp', 'spacetime', 'tutorial', 'M++'):
            self.PROJECT_ROOT_DIR = self.MPP_ROOT_DIR
        else:
            self.PROJECT_ROOT_DIR = abspath(join(self.MPP_ROOT_DIR, '..'))

        self.MPP_BUILD_DIR = abspath(join(self.MPP_ROOT_DIR, build_dir))
        self.PROJECT_BUILD_DIR = abspath(join(self.PROJECT_ROOT_DIR, build_dir))

        self.cmake_cache = join(self.PROJECT_BUILD_DIR, 'CMakeCache.txt')

        self.mpp_bench_file_path = join(self.PROJECT_BUILD_DIR, 'mpp_bench.txt')
        self.mpi_tests_file_path = join(self.PROJECT_BUILD_DIR, 'mpi_tests.txt')

    @property
    def MPP_LOG_DIR(self):
        return abspath(join(self.MPP_BUILD_DIR, 'log'))

    @property
    def MPP_JSON_DIR(self):
        return abspath(join(self.MPP_BUILD_DIR, 'json'))

    @property
    def MPP_DATA_DIR(self):
        return abspath(join(self.MPP_BUILD_DIR, 'data'))

    @property
    def MPP_VTK_DATA_DIR(self):
        return abspath(join(self.MPP_DATA_DIR, 'vtk'))

    @property
    def MPP_VTU_DATA_DIR(self):
        return abspath(join(self.MPP_DATA_DIR, 'vtu'))

    @property
    def MPP_PY_DATA_DIR(self):
        return abspath(join(self.MPP_DATA_DIR, 'py'))

    @property
    def PROJECT_LOG_DIR(self):
        return abspath(join(self.PROJECT_BUILD_DIR, 'log'))

    @property
    def PROJECT_JSON_DIR(self):
        return abspath(join(self.PROJECT_BUILD_DIR, 'json'))

    @property
    def PROJECT_DATA_DIR(self):
        return abspath(join(self.PROJECT_BUILD_DIR, 'data'))

    @property
    def PROJECT_VTK_DATA_DIR(self):
        return abspath(join(self.PROJECT_DATA_DIR, 'vtk'))

    @property
    def PROJECT_VTU_DATA_DIR(self):
        return abspath(join(self.PROJECT_DATA_DIR, 'vtu'))

    @property
    def PROJECT_PY_DATA_DIR(self):
        return abspath(join(self.PROJECT_DATA_DIR, 'py'))


def remove_file(file):
    if isfile(file):
        try:
            unlink(file)
        except Exception as e:
            print(e)


def clean_directory(directory, recursive):
    if isdir(directory):
        for file in listdir(directory):
            file_path = join(directory, file)
            try:
                if isfile(file_path):
                    unlink(file_path)
                elif isdir(file_path) and recursive:
                    rmtree(file_path)  # subdirs
            except Exception as e:
                print(e)


if __name__ == "__main__":
    dm = DirectoryManager(project_name='mpp', mute=False)
