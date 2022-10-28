# Changelog

## 0.1.0 (2022-10-28)


### Features

* add release please workflo ([6735e75](https://www.github.com/hydra-genetics/filtering/commit/6735e75406362e5bc62f5db5eede67556f8c1da7))
* add test and functions used to filter on both VEP and FORMAT. ([6eb13ec](https://www.github.com/hydra-genetics/filtering/commit/6eb13ec1e2a48f007156e786d982df649990ab52))
* Added rules: background_filter, artifact_filter ([a2ece03](https://www.github.com/hydra-genetics/filtering/commit/a2ece033b9b1c805962e4018f41aa0135da9716d))
* **ci:** add workflow for semantic commits ([b00de0d](https://www.github.com/hydra-genetics/filtering/commit/b00de0def422f62c84bfcc1f316cc4a5685d30c4))
* **ci:** pull-request template ([986b46c](https://www.github.com/hydra-genetics/filtering/commit/986b46cbd8760b7599ae2f3482941d3e575d05ff))
* filtering config is now configurable ([11a6343](https://www.github.com/hydra-genetics/filtering/commit/11a6343489c8d04e4a366963bc4bed87f6bdbb56))
* filtering rule now working ([c43ca2a](https://www.github.com/hydra-genetics/filtering/commit/c43ca2a62aa404912c7621f3d4cf44b5a1d475e8))
* functions used to convert a string to a nested list. ([11ffd7d](https://www.github.com/hydra-genetics/filtering/commit/11ffd7d0ca6c754404fc0d2b15231aed33a8d978))
* hard filter vcf file based on bed file using bcftools ([17b66fb](https://www.github.com/hydra-genetics/filtering/commit/17b66fb09cc327b4ca8b33df49607c645cb0e5a3))
* make config.yaml location more flexible ([3f79a47](https://www.github.com/hydra-genetics/filtering/commit/3f79a475d70809d7147015967fdf67d1485b479a))
* make configfile/confgilefiles argument mandatory ([73e6f30](https://www.github.com/hydra-genetics/filtering/commit/73e6f3022d71fa9056d6420cb0e6a19ee46a2452))
* make it possible to configure how NA data should be handled. ([e2c9cdb](https://www.github.com/hydra-genetics/filtering/commit/e2c9cdb8c47063ce32d9345c9a368c957cde4891))
* make tag more specific ([f0e6717](https://www.github.com/hydra-genetics/filtering/commit/f0e6717d1a472d937ddcb2ddeb26628535c61e9a))
* Moved rule: filter_vcf_on_format ([03851f0](https://www.github.com/hydra-genetics/filtering/commit/03851f0c0f1280e48a946c1153a15bbcf6db8b64))
* new rule: add_multi_snv_in_codon ([afe81cc](https://www.github.com/hydra-genetics/filtering/commit/afe81cc6d316df3eb2b086b6fe1ef2cae762f314))
* New rules: sort_vcf, bgzip_vcf, tabix_vcf ([c89886d](https://www.github.com/hydra-genetics/filtering/commit/c89886dd2e453110da0724a21821cbf7378c97ce))
* possible to extract INFO and exist exist and !exist as filtering. ([387f672](https://www.github.com/hydra-genetics/filtering/commit/387f672ce378c4789b7a27ea6d3b52729095bacc))
* print variant and filter that generate a TypeError ([3aec16e](https://www.github.com/hydra-genetics/filtering/commit/3aec16ee5277ff3a42a46f27264155de00beee3e))
* rule name split into soft and hard. Also some bugfixes ([736c616](https://www.github.com/hydra-genetics/filtering/commit/736c616909bcb5c1eefb2e2d28f80471efdd63e4))
* update functions and add tests. ([020429f](https://www.github.com/hydra-genetics/filtering/commit/020429fbd9527343216bb29bb11768ff1bbc5aa6))
* update README with correct dag graph ([aacf425](https://www.github.com/hydra-genetics/filtering/commit/aacf4257534e807105f9822849a1e4c7d073b6ec))
* update snakemake-version ([207e652](https://www.github.com/hydra-genetics/filtering/commit/207e65214b6845f5a58346489bd941ead0281cf9))
* update to latest hydra-genetics tools ([0dffe51](https://www.github.com/hydra-genetics/filtering/commit/0dffe51ace6afae857941d2674c3c4024680c58c))
* update to tools version 0.6.0 ([92b6022](https://www.github.com/hydra-genetics/filtering/commit/92b6022ad1ece221be0280bd2cb2f5196ffe71aa))


### Bug Fixes

* add_multi_snv_in_codon now takes a boolean flag instead of string ([5979fd4](https://www.github.com/hydra-genetics/filtering/commit/5979fd4031d8201e330a5e568ec394315478d3b7))
* added to schema ([b00fb8f](https://www.github.com/hydra-genetics/filtering/commit/b00fb8fd78dc1885e56bb9f6e0872ba228216100))
* change from match to search, will make it possible to locate substring in the middle of a string ([9750f26](https://www.github.com/hydra-genetics/filtering/commit/9750f267679d944b0fdae6759d4ee1108ed3649d))
* change run column name to flowcell in units. ([8b1f770](https://www.github.com/hydra-genetics/filtering/commit/8b1f7704502946267272e902c58e93401203ea36))
* changed formatting in vcf files ([d6c5688](https://www.github.com/hydra-genetics/filtering/commit/d6c56883bc155d096231690c065342c12fd969af))
* changed formatting in vcf files ([6f476ad](https://www.github.com/hydra-genetics/filtering/commit/6f476adbcc9d921cc92572a126785b8092cc2fde))
* comment spelling error ([bd5fd7e](https://www.github.com/hydra-genetics/filtering/commit/bd5fd7e81bff65748dcfb9d436d7d2631b47318f))
* common include firts ([dab16c6](https://www.github.com/hydra-genetics/filtering/commit/dab16c662af844b38482cf8c8449e6c6f91ef6be))
* config bugfix ([67ad5a9](https://www.github.com/hydra-genetics/filtering/commit/67ad5a9229f0dbb053de7681f3ffd68e303c69ad))
* handle case where tuple is return and no index is specified. ([355b81a](https://www.github.com/hydra-genetics/filtering/commit/355b81a302bfafb1c238e6e0820e6c1159ebc389))
* harmonized config file ([efe746c](https://www.github.com/hydra-genetics/filtering/commit/efe746c155c08413a1a7b80ca3a0c5c93a4cf694))
* missing conda env dependency ([61fde94](https://www.github.com/hydra-genetics/filtering/commit/61fde941e5b64dfb5a2e02f79fc4e699e957087f))
* missing samtools in conda. Also added missing containers in config ([d2495d3](https://www.github.com/hydra-genetics/filtering/commit/d2495d32cf30e15196d0efe839682fd298786032))
* output file names ([b94ee6a](https://www.github.com/hydra-genetics/filtering/commit/b94ee6a9be0c7ac3dd9e87f9eb4b429698797a59))
* removed broken bcftools wrapper. Fixed test output ([6465fd8](https://www.github.com/hydra-genetics/filtering/commit/6465fd8da3d6ea0f1c6e383f72726c2446462764))
* removed unneeded output file ([2340dde](https://www.github.com/hydra-genetics/filtering/commit/2340dde55991334ee313ec8ceaf300e471c75997))
* rename workflow file ([c0ce2a3](https://www.github.com/hydra-genetics/filtering/commit/c0ce2a3372b8a72e53eb1c487d681cbfb635c3bd))
* schemas ([ece5e2d](https://www.github.com/hydra-genetics/filtering/commit/ece5e2daeb6c5ba2d5a82ba1561da48f49f12ad6))
* tests and field handling ([afd4bac](https://www.github.com/hydra-genetics/filtering/commit/afd4bacabe5c2ba06edb8a79e39c1ba78b46971c))
* typo in config. Removed test files ([fc64015](https://www.github.com/hydra-genetics/filtering/commit/fc64015f089ccf292ca3e160d415872cf36f83cb))
* updated conda in rule ([b9cc5d8](https://www.github.com/hydra-genetics/filtering/commit/b9cc5d880eba2b3296c4c9cd6ca86f8740f6530a))
* wrong module name ([1ef8064](https://www.github.com/hydra-genetics/filtering/commit/1ef80643a445da20f50a0b2f7d0a782eef9796a3))


### Documentation

* update config info ([9c20b2d](https://www.github.com/hydra-genetics/filtering/commit/9c20b2d142a7fb277e5c46976e2c0433876a2599))
