class test_colors:
    OK = '\033[92m' #GREEN
    WARNING = '\033[93m' #YELLOW
    FAIL = '\033[91m' #RED
    RESET = '\033[0m' #RESET COLOR


def prRed(text):
    print(test_colors.FAIL + text + test_colors.RESET)


def prGreen(text):
    print(test_colors.OK + text + test_colors.RESET)


def print_test_results(duration, index, name, rc, log, total_num_of_tests):
    if rc == 0:
        message = f'{index + 1:>3}/{total_num_of_tests:<3}' \
                  f' Test  {index + 1:>3}: {name} '
        message = message.ljust(64, '.')
        message = message + '   Passed   {:3.2f} sec'.format(duration)
        prGreen(message)
    else:
        message = f'{index + 1:>3}/{total_num_of_tests:<3}' \
                  f' Test {index + 1:>3}: {name} '
        message = message.ljust(64, '.')
        message = message + '   Failed   {:3.2f} sec'.format(duration)
        prRed(message)
        prRed(log)