# Minimal setup.py for C extension building only
# Main configuration is in pyproject.toml

from setuptools import setup, Extension

def get_extensions():
    """Build C extensions with NumPy 2.0 compatibility."""
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
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
            language='c'
        )
        return [sasa_ext]
    except ImportError:
        # If numpy is not available during setup, return empty list
        # Extension will be built later when numpy is available
        print("Warning: NumPy not found during setup. Extension will be built later.")
        return []

# Setup only handles C extensions
# All other configuration is in pyproject.toml
setup(
    ext_modules=get_extensions(),
)