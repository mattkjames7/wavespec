import setuptools
from setuptools.command.build_py import build_py as _build_py
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
from setuptools.dist import Distribution
import os
import subprocess
import shutil
import platform
import glob
import re

with open("README.md", "r") as fh:
    long_description = fh.read()


def getversion():
    """Read the package version from wavespec/__init__.py."""
    thispath = os.path.abspath(os.path.dirname(__file__))
    initfile = os.path.join(thispath, 'wavespec', '__init__.py')

    with open(initfile, 'r') as f:
        lines = f.readlines()

    version = 'unknown'
    for line in lines:
        if '__version__' in line:
            parts = line.split('=')
            version = parts[-1].strip().strip('"').strip("'")
            break
    return version


version = getversion()


def _get_macos_arch():
    """Return a single target macOS arch when one is clearly specified."""
    # Explicit project override takes precedence.
    arch = os.environ.get('WAVESPEC_MACOS_ARCH', '').strip()
    if arch:
        return arch

    # Respect explicit CMake architecture if provided by caller/CI.
    cmake_arch = os.environ.get('CMAKE_OSX_ARCHITECTURES', '').strip()
    if cmake_arch:
        # CMake accepts semicolon-separated arch values.
        parts = [p.strip() for p in cmake_arch.split(';') if p.strip()]
        if len(parts) == 1:
            return parts[0]

    # Parse ARCHFLAGS like: "-arch arm64".
    archflags = os.environ.get('ARCHFLAGS', '')
    matches = re.findall(r'-arch\s+([A-Za-z0-9_]+)', archflags)
    unique = sorted(set(matches))
    if len(unique) == 1:
        return unique[0]

    return None


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
                cmake_cmd = ['cmake', '-S', data_dir, '-B', build_dir]
                if platform.system() == 'Darwin':
                    arch = _get_macos_arch()
                    if arch:
                        cmake_cmd.append(f'-DCMAKE_OSX_ARCHITECTURES={arch}')
                subprocess.check_call(cmake_cmd)
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


class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        # This package includes platform-specific native libraries.
        self.root_is_pure = False

        if platform.system() == 'Darwin' and not self.plat_name:
            arch = _get_macos_arch()
            if arch in ('arm64', 'x86_64'):
                major = platform.mac_ver()[0].split('.')[0] or '14'
                self.plat_name = f'macosx_{major}_0_{arch}'


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        # Force platlib wheel layout for bundled native shared libraries.
        return True


setuptools.setup(
    name="wavespec",
    version=version,
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
        'bdist_wheel': bdist_wheel,
    },
    distclass=BinaryDistribution,
)



