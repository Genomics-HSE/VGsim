cd src

rm -r build
rm -r dist
rm -r VGsim_test.egg*

sudo apt-get install python3-setuptools
python3 -m pip install --upgrade setuptools

python3 setup.py clean
python3 setup.py build
python3 setup.py install

cd ../