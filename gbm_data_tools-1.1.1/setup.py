import os
import shutil
import ssl
import sys
from urllib.request import urlretrieve

from setuptools import setup

# Import the version information and other package variables from the gbm package init file
with open('gbm/__init__.py') as f:
    exec(f.read())

# attempt to install GEOS
def geos_install():
    print("Checking if GEOS library is installed...")
    # custom install path
    home_path = os.path.expanduser('~')
    geos_path = os.environ.get('GEOS_DIR', 'PATH_NOT_SET')

    # check if already installed
    paths = ['/usr', '/usr/local', '/sw', '/opt', '/opt/local',
             home_path, geos_path]
    installed = any([os.path.exists(os.path.join(path, 'include', 'geos_c.h')) \
                     for path in paths])
    if installed:
        print("yes\n")
        return
    print("no\n")

    # create our custom path for GEOS install
    if os.path.basename(geos_path) is 'PATH_NOT_SET':
        return
    if not os.path.exists(geos_path):
        os.makedirs(geos_path)
    print("Custom install of GEOS library in {}".format(geos_path))

    cwd = os.getcwd()
    os.chdir(geos_path)
    
    # download basemap
    print("Downloading basemap...", file=sys.stderr)
    file = os.path.join(geos_path, 'v1.2.1rel.tar.gz')
    ssl._create_default_https_context = ssl._create_unverified_context
    urlretrieve('https://github.com/matplotlib/basemap/archive/v1.2.1rel.tar.gz',
                file)

    # unzip, navigate to GEOS
    print("Decompressing basemap...")
    os.system('tar xvf v1.2.1rel.tar.gz')
    os.chdir(os.path.join(geos_path, 'basemap-1.2.1rel', 'geos-3.3.3'))

    # install GEOS
    print("Installing GEOS library...")
    os.system('./configure --prefix=$GEOS_DIR')
    os.system('make; make install')

    # cleanup
    print("Cleaning up...")
    os.chdir('../..')
    os.system('rm -r basemap-1.2.1rel/')
    os.system('rm v1.2.1rel.tar.gz')
    os.chdir(cwd)
    print("Done.")

# attempt to install 'extra' requirement basemap 
def basemap_install():
    geos_install()
    return 'basemap @ git+https://github.com/matplotlib/basemap@v1.2.1rel'

# the test data files
def test_files():
    test_data_dir = 'gbm/test/data'
    test_data = os.listdir(test_data_dir)
    test_data = [os.path.join('data', file) for file in test_data if not \
        file.endswith('160509374')]
    test_burst = os.listdir(os.path.join(test_data_dir, '160509374'))
    test_burst = [os.path.join('data', '160509374', file) for file in test_burst]
    test_data.extend(test_burst)
    return test_data


# the documentation files
def doc_files():
    doc_dir = 'doc/_build/html'
    docs = []
    for root, dirs, files in os.walk(doc_dir):
        for name in files:
            docs.append(os.path.join(root, name))
    return docs


# the notebooks
def tutorial_files():
    tutorial_dir = 'tutorials/'
    tutorials = os.listdir(tutorial_dir)
    tutorials = [os.path.join(tutorial_dir, tutorial) for tutorial in tutorials]
    tutorials = [tutorial for tutorial in tutorials if tutorial.endswith('ipynb')]
    return tutorials


# make the home directory
def mkdir_home():
    # if it exists, do nothing
    if os.path.exists(home_path):
        return

    os.makedirs(home_path)
    os.mkdir(help_path)
    os.mkdir(tutorial_path)
    os.mkdir(data_path)

# copy API docs
def copy_docs(docs):
    path_prefix = os.path.commonprefix(docs)
    rel_paths = [os.path.relpath(doc, path_prefix) for doc in docs]
    for doc, rel_path in zip(*(docs, rel_paths)):
        try:
            shutil.copyfile(doc, os.path.join(help_path, rel_path))
        except:
            os.makedirs(os.path.join(help_path, os.path.dirname(rel_path)))
            shutil.copy(doc, os.path.join(help_path, rel_path))


def requirement_control():
    if sys.version_info.minor == 5:
        requirements = ['pyproj==1.9.6', 'numpy>=1.17.3, <= 1.18.4', 
                        'scipy>=1.1.0', 'matplotlib>=2.2.3, <=3.2.1',
                        'astropy>=3.1', 'healpy>=1.12.4, <=1.13.0']
    else:
        requirements = ['pyproj>=1.9.6', 'numpy>=1.17.3', 
                        'scipy>=1.1.0', 'matplotlib>=2.2.3, <=3.2.1',
                        'astropy>=3.1', 'healpy>=1.12.4']
    return requirements
        

# Import the version information from the gbm package init file
with open('gbm/__init__.py') as f:
    exec(f.read())

# make all the directories and copy the files over
mkdir_home()
copy_docs(doc_files())
[shutil.copy(tutorial_file, tutorial_path) for tutorial_file in tutorial_files()]

test_data = test_files()

setup(
    name='gbm-data-tools',
    version=__version__,
    scripts=['scripts/gbm-demos', 'scripts/gbm-docs'],
    packages=['gbm', 'gbm.background', 'gbm.binning', 'gbm.data',
              'gbm.lookup', 'gbm.plot', 'gbm.simulate', 'gbm.spectra', 
              'gbm.test'],
    url='',
    license='license.txt',
    author='Fermi GBM',
    author_email='',
    description='The Fermi GBM Data Tools',
    package_data={'gbm': ['McIlwainL_Coeffs.npy'], 'gbm.test': test_data},
    include_package_data=True,
    python_requires='>=3.5',
    install_requires=requirement_control(),
    extras_require={'basemap': [basemap_install()]},
    dependency_links=['https://github.com/matplotlib/basemap/archive/v1.2.1rel.tar.gz']
)
