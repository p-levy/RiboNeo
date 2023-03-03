import subprocess
import logging
import sys

def exec_command(cmd):
    logger = logging.getLogger()
    logger.info(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, executable='/bin/bash')
    output, error = p.communicate()
    if p.returncode != 0:
        for line in output.decode("utf-8").split("\n") if output else "":
            logger.error(line.rstrip())
        for line in error.decode("utf-8").split("\n") if error else "":
            logger.error(line.rstrip())
        sys.exit()
