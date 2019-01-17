# SOPHiA 

```python
cd Sophia
python setup.py install
```
if you receive an error due to lzma.h (OSX):
```bash
export HTSLIB_CONFIGURE_OPTIONS=--disable-lzma
pip install pysam
```

To run and generate data files and visualizations:
```python
python generate_report.py
```
