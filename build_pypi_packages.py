import subprocess as sp
import shutil
import os

sp.run('python -m build --sdist --no-isolation'.split())

with open('focus_tools/__init__.py') as f:
    exec(f.read())

src = f'dist/focus_package-{__version__}-py3-none-any.whl'

platforms = (
    ('win32', 'win_amd64'),
    ('linux', 'manylinux1_x86_64'),
    ('darwin','macosx_10_9_x86_64'),
    ('darwin','macosx_11_0_arm64'),
)

for plat, tag in platforms:
    env = {
        **os.environ,
        "PLATFORM": plat,
    }

    sp.run(
        'python -m build --wheel --no-isolation'.split(), 
        env=env
    )

    dst = src.replace('any', tag)

    shutil.move(src, dst)

