#from distutils.core import setup
from setuptools import setup
import os

import hywaves

def _strip_comments(l):
    return l.split('#', 1)[0].strip()

def _pip_requirement(req):
    if req.startswith('-r '):
        _, path = req.split()
        return reqs(*path.split('/'))
    return [req]

def _reqs(*f):
    return [
        _pip_requirement(r) for r in (
            _strip_comments(l) for l in open(
                #os.path.join(os.getcwd(), 'requirements', *f)).readlines()
                os.path.join(os.getcwd(), *f)).readlines()
        ) if r]

def reqs(*f):
    """Parse requirement file.
    Returns:
        List[str]: list of requirements specified in the file.
    Example:
        reqs('default.txt')          # requirements/default.txt
        reqs('extras', 'redis.txt')  # requirements/extras/redis.txt
    """
    return [req for subreq in _reqs(*f) for req in subreq]

def install_requires():
    """Get list of requirements required for installation."""
    return reqs('requirements.txt')

setup(
    name             = 'hywaves',
    version          = hywaves.__version__,
    description      = hywaves.__description__,
    long_description = open('README.md').read(),
    keywords         = hywaves.__keywords__,
    author           = hywaves.__author__,
    author_email     = hywaves.__contact__,
    url              = hywaves.__url__,
    license          = 'LICENSE.txt',
    python_requires  = ">=3.6",
    install_requires = install_requires(),
    packages         = ['hywaves', 'hywaves.statistical', 'hywaves.wswan',
                        'hywaves.wswan.plots', 'hywaves.statistical.plots',
                       ],
    include_package_data = True,
    package_data     = {
        'hywaves.wswan' : ['hywaves/wswan/resources/swan_bin/*'],
    },
    scripts          = [],
)

