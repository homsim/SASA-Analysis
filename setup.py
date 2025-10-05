"""SASA-Analysis

Repeatedly executes LAMMPS instances to probe the solvent-accessible surface-area
of a given macromolecule in order to create an energy landscape of that molecule.
"""


DOCSTRING = __doc__.split("\n")

from setuptools import setup, find_packages, Extension

def get_extensions():
    """Build C extensions, importing numpy only when needed."""
    try:
        import numpy
        sasa_ext = Extension(
            'sasa_ext.sasa_ext',
            sources=[
                'sasa_ext/sasa_module.c',
                'sasa_ext/sasa_core.c'
            ],
            include_dirs=[
                numpy.get_include(),
                'sasa_ext'
            ],
            extra_compile_args=['-O3', '-std=c99'],
            language='c'
        )
        return [sasa_ext]
    except ImportError:
        # If numpy is not available during setup, return empty list
        # Extension will be built later when numpy is available
        print("Warning: NumPy not found during setup. Extension will be built later.")
        return []

requirements = [
    'numpy>=1.15.0',
    'ovito>=3.0.0',
    'tqdm',
]

setup(
    name="sasa_lammps",
    version="0.3.0",
    author="SASA-Analysis Contributors",
    description="LAMMPS-based solvent accessible surface area analysis",
    long_description=__doc__,
    install_requires=requirements,
    include_package_data=True,
    packages=find_packages(),
    ext_modules=get_extensions(),
    setup_requires=['numpy>=1.15.0'],
    zip_safe=False,
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
