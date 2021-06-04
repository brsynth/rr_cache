# Release history

## 0.8.0
#### 4 Jun 2021
- feat: add --compounds CLI arguments
- feat: add get_compound(), get_list_of_compound()
- feat: get_reaction(), get_list_of_reactions()
- feat: get_reaction_rules(), get_list_of_reaction_rules()

## 0.7.1
#### 10 May 2021
- put back checksum controls

## 0.7.0
#### 10 May 2021
- remove deprecated infos
- BREAK: rename 'rr_full_reactions' into 'template_reactions'
- add cache-dir argument for generating or reading the cache

## 0.6.1
#### 4 May 2021
- fix: change sha for 'cid_strc' and 'rr_full_reactions' files according to previous changes

## 0.6.0
#### 3 May 2021
- add both deprecated reactions IDs and new ones in 'rr_full_reactions'

## 0.5.0
#### 3 May 2021
- add manual substitutions of IDs with no structure (from 'MNXM_replacement_20190524.csv')

## 0.4.0
#### 16 April 2021
- add '--reactions' option to print out rrCache content

## 0.3.0
#### 16 April 2021
- add '--reaction-rules' option to print out rrCache content

## 0.2.0
- make rrCache attributes private
- add get() and set() methods to read/write rrCache attributes

## 0.1.0
- split the loading of cache from the builder

## 0.0.3
- change printout informations

## 0.0.2
- change input command-line arguments
- change printout informations

## 0.0.1
- first commit
