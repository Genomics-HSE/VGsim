# ViralSim

File with rates is required.

Example:
```
python setup.py build_ext --inplace
./ViralSim.py example/example.rt -it 100000000 -pm example/example.pp example/example.mg -su example/example.su

for n in {1..5}; do ./ViralSim.py example/example.rt -it 1000000 -pm example/example.pp example/example.mg -su example/example.su; done >> errors.txt
```
