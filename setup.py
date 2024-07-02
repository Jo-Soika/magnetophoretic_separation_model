import subprocess
import os

def install_packages():
    if os.path.exists('requirements.txt'):
        subprocess.check_call(['pip', 'install', '-r', 'requirements.txt'])

install_packages()