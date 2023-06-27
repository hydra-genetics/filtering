# Hydra-genetics filter module

![dag plot](images/all.png){: style="height:40%;width:40%"}

<br />

The filter module works a litle differently compared to other hydra-genetics modules. It's not really a
pipeline, It's more a collection of rules, then a pipeline, that can be triggered by inserting ´{TYPE_OF_FILTER}.{TAGS}´ 
between the filename and vcf.gz, ex sample1.filter.MY_TAG.vcf.gz. 

Possible ´{TYPE_OF_FILTER}´ß are:
 - exclude: utilizes bcftools (bcftools_filter_exclude_region)
 - include: utilizes bcftools (bcftools_filter_include_region)
 - filter: custom script (filter_vcf)

 ´{MY_TAG}´ must be configured in the config file under the following intries:

 - bcftools_filter_include_region
 - bcftools_filter_exclude_region
 - filter_vcf

 which one depends on the rule you want to use

 ex:
 '''yaml
   

  bcftools_filter_include_region:
    EXONS: design/exons.bed

  filter_vcf:
    VEP_FILTERING: filter/vep_filtering.yaml

 '''
