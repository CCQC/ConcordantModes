import re
import subprocess
import time


class Submit(object):
    def __init__(self, disp_list):
        self.disp_list = disp_list

    def run(self):
        pipe = subprocess.PIPE

        process = subprocess.run(
            "qsub displacements.sh", stdout=pipe, stderr=pipe, shell=True
        )
        self.out_regex = re.compile(r"Your\s*job\-array\s*(\d*)")
        self.job_id = int(re.search(self.out_regex, str(process.stdout)).group(1))

        self.job_fin_regex = re.compile(r"taskid")
        while True:
            qacct_proc = subprocess.run(
                ["qacct", "-j", str(self.job_id)], stdout=pipe, stderr=pipe
            )
            qacct_string = str(qacct_proc.stdout)
            job_match = re.findall(self.job_fin_regex, qacct_string)
            if len(job_match) == len(self.disp_list):
                break
            time.sleep(10)

        output = str(process.stdout)
        error = str(process.stderr)
        pass
