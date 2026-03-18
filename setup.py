import setuptools
from setuptools.command.build_py import build_py as _build_py
import os
import subprocess
import shutil
import platform
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()


class build_py(_build_py):
    """Custom build_py that runs CMake to build bundled shared libs

    This runs before the Python build step so wheels/sdists include the
    compiled libraries (if CMake is available). Failures do not abort the
    overall build but are printed so packagers can see the issue.
    """

    def run(self):
        here = os.path.abspath(os.path.dirname(__file__))
        data_dir = os.path.join(here, 'wavespec', '__data')

        if os.path.isdir(data_dir):
            build_dir = os.path.join(data_dir, 'build')
            try:
                # Configure and build with CMake if available
                subprocess.check_call(['cmake', '-S', data_dir, '-B', build_dir])
                subprocess.check_call(['cmake', '--build', build_dir, '--config', 'Release'])
            except Exception as e:
                raise RuntimeError(f'CMake build for bundled libraries failed: {e}')

            systype = platform.system()
            if systype == 'Linux':
                ext = 'so'
            elif systype == 'Darwin':
                ext = 'dylib'
            elif systype == 'Windows':
                ext = 'dll'
            else:
                ext = 'so'

            # copy resulting libs into runtime locations expected by Globals.py
            cfg_dirs = ['Release', 'RelWithDebInfo', 'Debug', 'MinSizeRel']
            targets = [
                (
                    [
                        os.path.join(build_dir, 'liblombscargle', f'liblombscargle.{ext}'),
                        os.path.join(build_dir, 'liblombscargle', f'lombscargle.{ext}'),
                        *[
                            os.path.join(build_dir, 'liblombscargle', cfg, f'liblombscargle.{ext}')
                            for cfg in cfg_dirs
                        ],
                        *[
                            os.path.join(build_dir, 'liblombscargle', cfg, f'lombscargle.{ext}')
                            for cfg in cfg_dirs
                        ],
                        os.path.join(build_dir, f'liblombscargle.{ext}'),
                        os.path.join(build_dir, f'lombscargle.{ext}'),
                        os.path.join(data_dir, 'liblombscargle', 'build', f'liblombscargle.{ext}'),
                        os.path.join(data_dir, 'liblombscargle', 'build', f'lombscargle.{ext}'),
                    ],
                    os.path.join(data_dir, 'liblombscargle', f'liblombscargle.{ext}'),
                ),
                (
                    [
                        os.path.join(build_dir, 'libfilter', f'libfilter.{ext}'),
                        os.path.join(build_dir, 'libfilter', f'filter.{ext}'),
                        *[
                            os.path.join(build_dir, 'libfilter', cfg, f'libfilter.{ext}')
                            for cfg in cfg_dirs
                        ],
                        *[
                            os.path.join(build_dir, 'libfilter', cfg, f'filter.{ext}')
                            for cfg in cfg_dirs
                        ],
                        os.path.join(build_dir, f'libfilter.{ext}'),
                        os.path.join(build_dir, f'filter.{ext}'),
                        os.path.join(data_dir, 'libfilter', 'build', f'libfilter.{ext}'),
                        os.path.join(data_dir, 'libfilter', 'build', f'filter.{ext}'),
                    ],
                    os.path.join(data_dir, 'libfilter', f'libfilter.{ext}'),
                ),
            ]

            missing = []
            for src_candidates, dest in targets:
                src = next((p for p in src_candidates if os.path.isfile(p)), None)
                if src is None:
                    missing.append(dest)
                    continue
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copy2(src, dest)
                os.chmod(dest, 0o755)
                print(f'Copied bundled lib {src} -> {dest}')

            if missing:
                raise RuntimeError(
                    'Bundled native libraries were not produced at packaging time: '\
                    + ', '.join(missing)
                )

            # prevent build/venv artifacts under package tree from being
            # discovered as package data during wheel assembly
            shutil.rmtree(build_dir, ignore_errors=True)
            shutil.rmtree(os.path.join(data_dir, 'libfilter', 'build'), ignore_errors=True)
            shutil.rmtree(os.path.join(data_dir, 'liblombscargle', 'build'), ignore_errors=True)
            for p in glob.glob(os.path.join(data_dir, '.venv*')):
                shutil.rmtree(p, ignore_errors=True)

            # remove stale top-level copies - runtime expects libs in subdirs
            for stale in [
                os.path.join(data_dir, 'libfilter.so'),
                os.path.join(data_dir, 'liblombscargle.so'),
                os.path.join(data_dir, 'libfilter.dylib'),
                os.path.join(data_dir, 'liblombscargle.dylib'),
                os.path.join(data_dir, 'libfilter.dll'),
                os.path.join(data_dir, 'liblombscargle.dll'),
            ]:
                if os.path.isfile(stale):
                    os.remove(stale)

        # continue with normal build
        super().run()


setuptools.setup(
    name="wavespec",
    version="0.0.6",
    author="Matthew Knight James",
    author_email="mattkjames7@gmail.com",
    description="Some spectral analysis tools for analyzing waves in data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mattkjames7/wavespec",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX",
    ],
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'DateTimeTools',
    ],
    include_package_data=False,
    package_data={
        'wavespec': [
            '__data/liblombscargle/*.so',
            '__data/libfilter/*.so',
            '__data/liblombscargle/*.dylib',
            '__data/libfilter/*.dylib',
            '__data/liblombscargle/*.dll',
            '__data/libfilter/*.dll',
        ],
    },
    cmdclass={
        'build_py': build_py,
    },
)



