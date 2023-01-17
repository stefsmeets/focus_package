from setuptools import setup
import os, sys

platform = os.getenv('PLATFORM', sys.platform)

if platform == "win32":
    data_files = ('Scripts', [
        "bin_windows/libfftw3f-3.dll",
        "bin_windows/focus.exe",
        "bin_windows/sginfo.exe",
        "bin_windows/_kriber.x",
        "bin_windows/_dls76.x",
    ])
elif platform == "darwin":
    data_files = ('bin', [
        "bin_osx/focus", 
        "bin_osx/sginfo",
        "bin_osx/_kriber.x",
        "bin_osx/_dls76.x",
    ])
elif platform == "linux":
    data_files = ('bin', [
        "bin_linux/focus", 
        "bin_linux/sginfo",
        "bin_linux/_kriber.x",
        "bin_linux/_dls76.x",
    ])
else:
    raise RuntimeError


setup(
    data_files = [
        data_files
    ],
)
