from setuptools import setup
import os, sys

platform = os.getenv('PLATFORM', sys.platform)

if platform == "win32":
    package_dir={"bin": "bin_windows/bin"}
    scripts = [
        "bin_windows/libfftw3f-3.dll",
        "bin_windows/focus.exe",
        "bin_windows/sginfo.exe",
    ]
elif platform == "darwin":
    package_dir={"bin": "bin_osx/bin"}
    scripts = [
        "bin_osx/focus", 
        "bin_osx/sginfo",
    ]
elif platform == "linux":
    package_dir={"bin": "bin_linux/bin"}
    scripts = [
        "bin_linux/focus", 
        "bin_linux/sginfo",
    ]
else:
    raise RuntimeError


setup(
    package_dir=package_dir,
    data_files = [
        ('bin', scripts)
    ],
)
