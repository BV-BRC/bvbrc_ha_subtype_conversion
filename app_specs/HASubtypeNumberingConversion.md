
# Application specification: HASubtypeNumberingConversion

This is the application specification for service with identifier HASubtypeNumberingConversion.

The backend script implementing the application is [App-HASubtypeNumberingConversion.pl](../service-scripts/App-HASubtypeNumberingConversion.pl).

The raw JSON file for this specification is [HASubtypeNumberingConversion.json](HASubtypeNumberingConversion.json).

This service performs the following task:   HA Subtype Numbering Conversion

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| input_source | Source of input (feature_list, fasta_data, fasta_file, genome_group) | enum  | :heavy_check_mark: |  |
| input_fasta_data | Input sequence in fasta formats | string  |  |  |
| input_fasta_file | Input sequence as a workspace file of fasta data | wsid  |  |  |
| input_feature_group | Input sequence as a workspace feature group | wsid  |  |  |
| input_feature_list | Input sequence as a list of feature ids | string  |  |  |
| types | Selected types in the submission | string  | :heavy_check_mark: |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |

