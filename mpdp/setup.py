from pathlib import Path
import re
from setuptools import setup, find_packages

#!/usr/bin/env python3
"""
Minimal, robust setup.py for packaging the project in this folder.

Usage:
    python setup.py sdist bdist_wheel
    pip install .
"""


ROOT = Path(__file__).parent

def read_text(path: Path, default: str = "") -> str:
    try:
        return path.read_text(encoding="utf8")
    except Exception:
        return default

# Attempt to get long description from README if present
long_description = read_text(ROOT / "README.md", default=read_text(ROOT / "README.rst", default=""))

# Try to extract a __version__ from the package's __init__.py without importing
def find_version(package: str) -> str:
    init_path = ROOT / package / "__init__.py"
    content = read_text(init_path)
    m = re.search(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", content, re.M)
    return m.group(1) if m else "0.0.0"

# Replace 'mpdp' below with the actual package directory name if different
PACKAGE_NAME = "mpdp"
VERSION = find_version(PACKAGE_NAME)

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description="A short description of the mpdp package",
    long_description=long_description,
    long_description_content_type="text/markdown" if (ROOT / "README.md").exists() else "text/plain",
    author="",
    author_email="",
    url="",
    packages=find_packages(exclude=("tests", "docs")),
    include_package_data=True,
    install_requires=[
        # "requests>=2.0",
    ],
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL3 License",
        "Operating System :: OS Independent",
    ],
    license="GPL3",
    zip_safe=False,
)