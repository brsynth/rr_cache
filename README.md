# rr_cache - RetroRules Cache

## Requirements
rr_cache has the following dependencies:
- `brs_utils`
- `requests`
- `rdkit`
- `colored`

These dependencies can be installed through conda package manager with:
```sh
conda install -c conda-forge brs_utils requests rdkit colored
```
or by creating a dedicated environment:
```sh
conda env create -f environment.yaml
```

## Install
### From Conda
```sh
[sudo] conda install -c brsynth -c conda-forge rr_cache
```

## Use

All cache data are stored into files on disk and loaded in memory each time the tool is used. Memory fingerprint is equal to the size of cache files loaded in memory multiplied by the number of processes which are running at the same time.

### Load rrCache in memory
```python
from rr_cache import rrCache

cache = rrCache()
print(cache.cid_strc)
```

**A part of cache**
For less loading time and memory footprint, a part of the cache can be loaded:
```python
from rr_cache import rrCache

cache = rrCache(['cid_strc'])
cache.get_compound('MNXM2')
```
***From CLI***
```sh
python -m rr_cache --reaction-rules <RuleID_1> <RuleID_2>...
```
If `RuleID_i` is(are) set, prints out the reaction rule(s) or all reaction rules if nothing is specified.

```sh
python -m rr_cache --reactions <RxnID_1> <RxnID_2>...
```
If `RxnID_i` is(are) set, prints out the reaction(s) or all reactions if nothing is specified.

```sh
python -m rr_cache --compounds <CmpdID_1> <CmpdID_2>...
```
If `CmpdID_i` is(are) set, prints out the compound(s) or all compounds if nothing is specified.

### (Re-)generate the cache
**From Python code**
```python
from rr_cache import rrCache

rrCache.generate_cache(outdir)
```

**From CLI**
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
