cd Sophia
python setup.py install

if you receive an error due to lzma.h:
    export HTSLIB_CONFIGURE_OPTIONS=--disable-lzma
    pip install pysam