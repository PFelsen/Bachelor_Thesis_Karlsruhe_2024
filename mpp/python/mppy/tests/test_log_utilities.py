from unittest import TestCase, main
import os
import sys
import pandas as pd
sys.path.append("..")
from utilities import log_utilities



class TestLogParser(TestCase):
    def setUp(self):
        py_file_dir = os.path.dirname(os.path.realpath(__file__))

        self.key_value_log = 'key_value_log'
        self.time_log = 'time_log'

        self.log_parser = log_utilities.LogParser()
        self.log_parser.log_dir = os.path.dirname(os.path.realpath(__file__))

    def tearDown(self) -> None:
        self.log_parser.reset_data()

    def test_read_log(self):
        self.log_parser.read_log(self.key_value_log)

    def test_parse_log(self):
        df = pd.DataFrame(self.log_parser.parse_log(self.key_value_log))
        # Todo This is a test for Bug reproduction: This should be Laplace2D
        # Todo -> conversion of strings with digits is not working !
        # Todo -> This has to be fixed with a new; cleaner parser
        # Todo -> wasn't done as it as a minor bug
        self.assertEqual(df.loc[0, 'Problem'], 2.0)  # Laplace2D != 2.0

        self.assertEqual(df.loc[0, 'Mesh Name'], 'UnitSquare')
        self.assertEqual(df.loc[0, 'Level'], 6)
        self.assertEqual(df.loc[0, 'Mesh width'], [0.015625, 0.015625])
        self.assertEqual(df.loc[0, 'initLevels'], [2, 3, 4])
        self.assertEqual(df.loc[0, 'initSamples'], [5])
        self.assertEqual(df.loc[0, "Time Series Name"], "UniformTimeSeries")
        self.assertEqual(df.loc[0, "t0"], 0)
        self.assertEqual(df.loc[0, "T"], 2)
        self.assertEqual(df.loc[0, "dt"], 0.03125)
        self.assertEqual(df.loc[0, "int and float"], [4.0, 2.0, 1.3, -4.0])
        self.assertEqual(df.loc[0, "list with spaces"], [4.0, 2.0, 1.3, -4.0])
        self.assertEqual(df.loc[0, "list with other spaces"], [4.0, 2.0, 1.3, -4.0])
        self.assertEqual(df.loc[0, "Table with spaces"], [4.0, 2.0, 1.3, -4.0])
        self.assertEqual(df.loc[0, "space-time Cells"], [4194304, 4194304, 4194304, 4194304, 4194304, 0])
        self.assertEqual(df.loc[0, "Computation Time"], "19.28 seconds")

    def test_time_log(self):
        df = pd.DataFrame(self.log_parser.parse_log(self.time_log))
        self.assertEqual(df.loc[0, 'n'], [0, 1, 2, 3, 4])
        self.assertEqual(df.loc[0, 't'], [0, 0.03125, 0.0625, 0.09375, 0.125])
        self.assertEqual(df.loc[0, 'Energy'], [4, 2.8140099, 2.2524386, 1.9247726, 1.7068081])
        self.assertEqual(df.loc[0, 'Mass'], [1, 1, 1, 1, 0.99999999])
        self.assertEqual(df.loc[0, 'Error'], [0, 0.95824483, 1.250116, 1.4509821, 1.6126332])

    def test_parse_log_second_time(self):
        self.log_parser.parse_log(self.key_value_log)
        self.log_parser.parse_log(self.key_value_log)
        df = pd.DataFrame(self.log_parser.data)
        self.assertEqual(list(df.loc[0:1, 'Mesh Name']), ['UnitSquare', 'UnitSquare'])
        self.assertEqual(list(df.loc[0:1, 'Level']), [6, 6])
        self.assertEqual(list(df.loc[0:1, 'Mesh width']), [[0.015625, 0.015625], [0.015625, 0.015625]])
        self.assertEqual(list(df.loc[0:1, 'initLevels']), [[2, 3, 4], [2, 3, 4]])
        self.assertEqual(list(df.loc[0:1, 'initSamples']), [[5], [5]])
        self.assertEqual(list(df.loc[0:1, "Time Series Name"]), ["UniformTimeSeries", "UniformTimeSeries"])
        self.assertEqual(list(df.loc[0:1, "t0"]), [0, 0])
        self.assertEqual(list(df.loc[0:1, "T"]), [2, 2])
        self.assertEqual(list(df.loc[0:1, "dt"]), [0.03125, 0.03125])
        self.assertEqual(list(df.loc[0:1, "Computation Time"]), ["19.28 seconds", "19.28 seconds"])


if __name__ == '__main__':
    main()
