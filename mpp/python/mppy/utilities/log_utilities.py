import re
import os


class LogParser:
    """
    Parser for mpp log files
    """

    float_re = r'([+-]?[0-9]+[.]?[0-9]*[eE]?[+-]?[0-9]*)'
    string_re = r'([a-zA-Z]+[a-zA-Z 0-9(_)+-]*)'
    dots_re = r': [\.]+ '
    seconds_re = r'([0-9]+[.]?[0-9]*[ ]*seconds)'
    minutes_re = r'([0-9]+[:]?[0-9]+[.]?[0-9]*[ ]*minutes)'
    hours_re = r'([0-9]+[:]?[0-9]+[:]?[0-9]+(?:\.[0-9]+)?[ ]*hours)'

    def __init__(self, log_dir=None):
        self.log_dir = log_dir
        self.data = {}

        # Finds floats and integers
        self.float_pattern = re.compile(self.float_re)

        # Example: [4.0, 2, 1.3]
        #self.vec_pattern = re.compile(r'\[[ ]*?(' + self.float_re + r'),?[ ]*?\]')
        self.vec_find_pattern = re.compile('[ \[]*' + self.float_re + '[, \]*]')

        # Finds strings with
        self.string_pattern = re.compile(self.string_re)

        # example: Error(u-u_h)=0.1234
        self.equation_pattern = re.compile(self.string_re + r'\s?=\s?' + self.float_re)

        # example: 01:23:02 hours
        self.time_pattern = re.compile(self.seconds_re + r'|' + self.minutes_re + r'|' + self.hours_re)

    def reset_data(self):
        self.data = {}

    def read_log(self, log_file='log'):
        if os.path.isfile(os.path.join(self.log_dir, log_file)):
            with open(os.path.join(self.log_dir, log_file)) as file:
                content = file.readlines()
            return content
        else:
            raise FileNotFoundError(log_file)

    def parse_log(self, log_file=None, additional_entries=None):
        if log_file == "all":
            for file in sorted(os.listdir(self.log_dir)):
                if file.endswith('.log'):
                    self.parse_single_log(file, additional_entries)
        else:
            self.parse_single_log(os.path.join(self.log_dir, log_file), additional_entries)
        return self.data

    def parse_single_log(self, log_file=None, additional_entries=None):
        log_content = [line.strip() for line in self.read_log(log_file)]

        if log_content[-1].find("end program") == -1:
            print("{} is of an unsuccessful computation.".format(log_file))
            return self.data

        new_data = {"tmp_lin_steps": []}

        for idx, line in enumerate(log_content):
            self._parse_line(line, new_data)

        if new_data["tmp_lin_steps"]:
            value = max(new_data["tmp_lin_steps"])  # always choose max amount of linear steps
            key = 'Linear Steps'
            new_data[key] = value
            new_data.pop("tmp_lin_steps")

        if additional_entries is not None:
            if type(additional_entries) is not dict:
                raise ValueError('additional_entries should be dictionary')
            else:
                for key in additional_entries:
                    new_data[key] = additional_entries[key]

        self.append_data(new_data)
        return self.data

    def _parse_line(self, line: str, data: dict):
        if line.find('end program after') != -1:
            value = "Time not found."
            time_matches = self.time_pattern.findall(line)[0]
            for time_match in time_matches:
                if time_match:
                    value = time_match
            data['Computation Time'] = value
            return
        elif "=" in line:
            for equation in self.equation_pattern.findall(line):
                key = str(self.string_pattern.findall(equation[0])[0])
                value = float(self.float_pattern.findall(equation[1])[0])
                self.add_to_dic(key, value, data)
            return
        elif line.find('GMRES: d(') != -1 or line.find('PCG: d(') != -1:
            value = int(self.float_pattern.findall(line)[0])
            data["tmp_lin_steps"].append(value)
            return

        splitted = re.split(self.dots_re, line, maxsplit=1)
        if len(splitted) != 2:
            return None
        key, value = [k.strip() for k in splitted]

        vec = self.vec_find_pattern.findall(value)
        if len(vec) > 0:
            value = [float(v.strip()) for v in vec if v.strip()]
            self.add_to_dic(key, value, data)
        else:
            floats = self.float_pattern.findall(value)
            if len(floats) > 0:
                try:
                    self.add_to_dic(key=key, value=float(floats[0]), data_dict=data)
                except ValueError:
                    self.add_to_dic(key=key, value=str(floats[0]), data_dict=data)
            else:
                string = self.string_pattern.findall(value)
                if len(string) > 0:
                    self.add_to_dic(key=key, value=string[0], data_dict=data)

    @staticmethod
    def add_to_dic(key, value, data_dict):
        if key is None:
            return data_dict
        key = key.strip()
        if key not in data_dict:
            data_dict[key] = value
        elif type(data_dict[key]) == list:
            data_dict[key].append(value)
        else:
            data_dict[key] = [data_dict[key], value]
        return data_dict

    def append_data(self, new_data):
        previous_runs = 0
        for key in self.data:
            previous_runs = len(self.data[key])  # all values should have same length
            if key not in new_data:
                new_data[key] = None
        for key in new_data:
            if key not in self.data:
                self.data[key] = [None for i in range(previous_runs)]
            self.data[key].append(new_data[key])

