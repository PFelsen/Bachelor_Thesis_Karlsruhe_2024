from subprocess import Popen, PIPE, STDOUT
import os


class SubprocessManager:
    def __init__(self, mute=False):
        self.mute = mute
        self.latest_log = ""
    def run_subprocess(self, args, cwd, mute=None):
        mute = mute if mute is not None else self.mute

        self.latest_log = ""

        if not os.path.exists(cwd):
            raise OSError('working directory ' + cwd + ' does not exist')

        if not isinstance(args, list):
            raise TypeError('expected args as list')

        process = Popen(args, cwd=cwd, stdout=PIPE, stderr=STDOUT)
        while True:
            stdout = process.stdout.readline()
            self.latest_log += stdout.decode('utf-8')
            if not mute:
                if stdout:
                    print(stdout.decode('utf-8').strip('\n'))
            if stdout == b'' and process.poll() is not None:
                break
        process.wait()
        process.stdout.close()

        rc = process.poll()
        return rc
