cd src

rm -r build
rm -r dist
rm -r VGsim_test.egg*

python3 setup.py clean
python3 setup.py build
python3 setup.py install

cd ../