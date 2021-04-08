# rr_cache

RetroRules Cache

## Memory management

### File mode
This is the default mode. All cache data are stored into files on disk and loaded in memory each time the tool is used. In this mode, fingerprint in memory is equal to the size of cache files loaded in memory multiplied by the number of processes which are running at the same time. Option can be specified by `--store-mode file`.

### DB mode
In order to save memory space, cache data can be loaded once in a database (redis) so that the memory space taken is equal to one instance of the cache, whatever the number of processes whic are running. Option can be specified by `--store-mode <db_host>`, where `db_host` is the hostname on which redis server is running.


## Install
### From Conda
```sh
[sudo] conda install -c brsynth -c conda-forge rr_cache
```

## Use

### Load rrCache in memory
**Full cache into files**
```python
from rr_cache import rrCache

cache = rrCache(db='file')
print(cache.cid_src)
```

**Full cache into Redis DB**
For multiple instances of rrCache simultaneously, rrCache can be loaded into one single Redis database:
```python
from rr_cache import rrCache

cache = rrCache(db='localhost')
print(cache.cid_src)
```
`localhost` means that rrCache will look for a redis database locally. If there is not, it will start a brand new redis server. `localhost` could be replaced by any hostname that hosts the Redis database.

**A part of cache**
For less loading time and memory footprint, a part of the cache can be loaded:
```python
from rr_cache import rrCache

cache = rrCache(attrs='cid_strc')
print(cache.cid_src)
```

### (Re-)generate the cache
**From Python code**
```python
from rr_cache import rrCache

rrCache.generate_cache(outdir)
```

**From CLI**
After having installed rptools Python module:
```sh
python -m rr_cache --gen_cache <folder>
```


## Test
Please follow instructions below ti run tests:
```
cd tests
pytest -v
```
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)).


## Authors

* **Joan Hérisson**
* **Melchior du Lac**

## Acknowledgments

* Thomas Duigou


## Licence
rr_cache is released under the MIT licence. See the LICENCE file for details.
