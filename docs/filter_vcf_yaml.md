# Filter_vcf yaml format

This hydra-genetics script uses yaml to define how to filter a vcf file.

## Filter format

The format is straightforward and is as follows:

```yaml
filters:
  FILTER_NAME1:                            # Name section
    description: "description of filter 1" # Description section
    expression: "expression1"              # Expression section
    soft_filter_flag: "flag1"              # Flag name section
    soft_filter: "True"                    # Soft/Hard filter setting 

  FILTER_NAME2:                            # Name section
    description: "description of filter 2" # Description section
    expression: "expression2"              # Flag name section
    soft_filter: "False"                   # Soft/Hard section
```

<br />

- **Name section**: a simple string used solely to distinguish different filters easily.
- **Description section**: a description used to describe the filtering that will be applied.
- **Expression section**: the actual filtering that will be performed. How to write this expression will be describe further down this page. **REQUIRED**
- **Flag name section**: flag added to the filter section. **REQUIRED** if soft_filter=True
- **Soft/Hard section**: this variable can be set to "True" for soft filtering, i.e the variant will have it's filter column updated. With false all variants matching the filter will be removed from the output file.  **REQUIRED**

## Expression format

Consist of a data extraction definition and a comparison operator/function, multiple expression can be combined using logical operators. Parentheses can be used to group expressions. 

Expression format:


- DATA_EXTRACTION [<|>|=|!=] VALUE
- VALUE [<|>|=|!=] DATA_EXTRACTION
- exist[regex, DATA_EXTRACTION]
- !exist[regex, DATA_EXTRACTION]

### Data extraction

The data extraction expect the following format: DATA_SOURCE:NA_HANDLING:FIELD:COLUMN

<br />

**DATA_SOURCE**: specify from which column or part of column that the data will be extracted:

| Value | Description |
| --- | --- |
| VEP | vep field in the info column used as data source |
| FORMAT | all fields in the format column |
| INFO | all fields in the info column |
| QUAL | all data in QUAL column |

<br />

**NA_HANDLING**: **OPTIONAL** : defines how an expression will behave if no data is found for the data_source or data_source, field combination. The following options are available: 

| Value | Description |
| --- | --- |
| NA_TRUE | when None is found the expression will return True and filter variant |
| NA_FALSE | when None is found the expression will return False and not filter variant |
| NA_ERROR | when None is found the an error will be raised if value not found, terminating the script |

 **Default** will be NA_FALSE, i.e a filter will not remove variant

<br />

 **FIELD**: any field in info, format or vep string. Not required for QUAL since the column doesn't contain key value content.

<br />

**COLUMN**: **OPTIONAL** : some fields in info and format will return a list of values and COLUMN makes it possible to specify which of these values that will be extracted, if COLUMN isn't set the entire list will be returned.

Example:


 - VEP:SYMBOL
 - INFO:Artifact:0 > 3
 - FORMAT:NA_TRUE:SB_mutect2:1
 

### Logical and comparison operators
Common logical and comparison operators are supported, can be combined with parentheses.

Logical operators:


 - and 
 - or

 **NOTE**: don't use capital letters


### Comparison operators

| Operator | Description |
| --- | --- |
| > | greater then |
| < | less then |
| = | equal to |
| != |  not equal to |

<br />

Example:


- FORMAT:SB_mutect2:1 > 400
- FORMAT:AD:1 > 10 and PLCH2 = VEP:SYMBOL
- QUAL < 75 and INFO:Artifact:0 > 3

### Comparison functions
In addition to the comparison operators a function named "exist" is available, can be used with "!" (i.e !exist). This function takes two arguments, a regex defining what to look for and a data source definition. The function is written like this `exist[regex, DATA_EXTRACTION]`

<br />

Example:


- exist[XM_[0-9]+, VEP:Feature]
- !exist[XM_[0-9]+, VEP:Feature]
- exist[XM_[0-9]+, VEP:Feature] and exist[mutect2, INFO:CALLERS] 

## Examples

**Soft filter**
```yaml
filters:
  vaf:
    description: "Soft filter variants with low vaf (AF lower than 0.01)"
    expression: "(FORMAT:AF:0 < 0.01)"
    soft_filter_flag: "AF_lt_0.01"
    soft_filter: "True"
  artifacts:
    description: "Soft filter variants found in more than 3 normal samples"
    expression: "(INFO:Artifact:0 > 3 or INFO:Artifact:1 > 3)"
    soft_filter_flag: "Artifact_gt_3"
    soft_filter: "True"
  background:
    description: "Soft filter position with where backgound distribution overlaps variant (lower than 4 SD from median)"
    expression: "INFO:PositionNrSD < 4 and !exist[1-hotspot, INFO:Hotspot]"
    soft_filter_flag: "Background_lt_4"
    soft_filter: "True"
  germline:
    description: "Soft filter germline SNVs based on GnomAD"
    expression: "(VEP:gnomAD_AF > 0.005)"
    soft_filter_flag: "Germline"
    soft_filter: "True"
  qual_soft:
    expression: "QUAL < 75"
    soft_filter: "True"
    soft_filter_flag: "qual"
```

<br />

**Hard filter**
``` yaml
filters:
  intron:
    description: "Hard filter intronic variants"
    expression: "(exist[intron_variant, VEP:Consequence] and !exist[splice, VEP:Consequence] and VEP:SYMBOL != MET and VEP:SYMBOL != TERT and !exist[COSV[0-9]+, VEP:Existing_variation])"
    soft_filter: "False"
  vaf:
    description: "Hard filter variants with low vaf (AF lower than 0.01)"
    expression: "(FORMAT:AF:0 < 0.01)"
    soft_filter: "False"
  artifacts:
    description: "Hard filter variants found in more than 3 normal samples"
    expression: "((INFO:Artifact:0 > 3 and INFO:ArtifactNrSD:0 < 5) or (INFO:Artifact:1 > 3 and INFO:ArtifactNrSD:1 < 5))"
    soft_filter: "False"
  background:
    description: "Hard filter position with where backgound distribution overlaps variant (lower than 4 SD from median)"
    expression: "(INFO:PositionNrSD < 4 and !exist[1-hotspot, INFO:Hotspot])"
    soft_filter: "False"
  germline:
    description: "Hard filter germline"
    expression: "(VEP:gnomAD_AF > 0.005)"
    soft_filter: "False"
```
