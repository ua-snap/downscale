# INSTALLATION STEPS
--
## Pre-install Notes:
This assumes that you have a functioning python3 installation, if not see: `https://github.com/EarthScientist/etc/blob/master/python_362_without_root.sh`

## clone the repo
`git clone git@github.com:ua-snap/downscale.git`

## enter the directory
`cd ./downscale`

## make a new virtual environment named venv
`~/.localpython/bin/python3.5 -m venv venv`

## source in (activate) the new virtualenv `venv` so we use that python installation
`source venv/bin/activate`

## install some stuff
```sh
pip install --upgrade pip
pip install ipython
pip install -r requirements.txt

## install the package
python setup.py install
```
## You now have downscale installed.

