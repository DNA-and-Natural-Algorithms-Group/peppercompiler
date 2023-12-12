import sys
import subprocess
import pkg_resources

def main():
    binary_path = pkg_resources.resource_filename(__name__,'_spuriousSSM')

    subprocess.run([binary_path]+sys.argv[1:])
