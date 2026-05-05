# GPhase detailed parameter reference

This page lists the full parameters for the pipeline scripts exposed by the main `gphase` program. The main README keeps only the most important parameters.

## `gphase pipeline`

Script: `pipeline/gphase_pipeline.sh`

Usage:
```
gphase pipeline -f <fa_file> -g <gfa> -c <collapse_num_file> -m <map_file> --n_chr <n_chr> --n_hap <n_hap> -p <output_prefix>
```

### Required parameters

| Parameter | Description |
| --- | --- |
| `-f <fa_file>` | FASTA file containing genome/unitig sequences. |
| `-g <gfa>` | GFA file representing the assembly graph. |
| `-c <collapse_num_file>` | Collapse number file for unitigs, usually from `popcnv/06.genes.round.cn`. |
| `-m <map_file>` | Hi-C/Pore-C/Omni-C/CiFi mapping file. Supported suffixes: `.bam`, `.pairs`. |
| `-p <output_prefix>` | Prefix for output files. Only `[a-zA-Z0-9.]` is allowed. |
| `--n_chr <n_chr>` | Number of chromosomes. Must be a positive integer. |
| `--n_hap <n_hap>` | Number of haplotypes. Must be a positive integer. |

### Optional parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `-e <enzyme_site>` | `GATC` | Restriction enzyme cutting site. |
| `--nor_hic <mode>` | `ratio` | Normalization mode for 3C link connections. Choices: `no`, `ratio`, `length`. |

### Preprocessing parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--cluster_q <cluster_q>` | `1` | Filtered MAPQ value for clustering. Enabled only when the input mapping file is BAM. |
| `--scaffold_q <scaffold_q>` | `0` | Filtered MAPQ value for scaffolding. Enabled only when the input mapping file is BAM. |

`--cluster_q` must be greater than or equal to `--scaffold_q`.

### Chromosome clustering parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--split_gfa_n <split_gfa_n>` | `5` | Number of common neighbors when splitting GFA. Valid range: `2-5`. |
| `--chr_pm <partig_chr_pm>` | `0.95` | Similarity of partig when clustering chromosomes. Valid range: `0.8 <= x < 1`. |
| `--r_max <r_max>` | `3` | Maximum value of parameter `R` during Louvain clustering. |
| `--t_len_T <t_len_T>` | `7` | Threshold for filtering the total length of the cluster. Set to `0` to disable filtering. |
| `--a_len_T <a_len_T>` | `7` | Threshold for filtering the average unitig length within a cluster. Set to `0` to disable filtering. |

### Haplotype clustering parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--hap_pm <partig_hap_pm>` | `0.7` | Similarity of partig when clustering haplotypes. Valid range: `0.6 <= x < 1`. |
| `--rescue` | off | Rescue the subgraph. |
| `--expand` | off | Expand the haplotype clustering process. |
| `--reassign_number <reassign_number>` | `1` | Number of reassignment steps. Valid values: `1`, `2`, `3`. |

### Haplotype scaffolding parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--thread <thread>` | `12` | Number of parallel processes. |
| `--no_contig_ec` | off | Disable contig error correction in YaHS. |
| `--no_scaffold_ec` | off | Disable scaffold error correction in YaHS. |
| `--min_len <min_len>` | `50` | Minimum scaffold length in kb in HapHiC sort. Valid range: `0-1000`. |
| `--mutprob <mutprob>` | `0.6` | Mutation probability. Valid range: `0.1-0.9`. |
| `--ngen <ngen>` | `20000` | Number of generations for GA. |
| `--npop <npop>` | `200` | Population size for GA. |
| `--processes <processes>` | `32` | Processes for fast sorting. |

## `gphase contact-pair`

Script: `pipeline/contact_pair_pipeline.sh`

Usage:
```
gphase contact-pair <genome.fa> <reads1.fq.gz> [reads2.fq.gz ...] [options]
```

### Required parameters

| Parameter | Description |
| --- | --- |
| `<genome.fa>` | Genome FASTA file. |
| `<reads*.fq.gz>` | One or more raw contact-pair reads in FASTQ format. |

### Optional parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `-o <prefix>` | `contact_pair` | Output directory/prefix. Final BAM path: `<prefix>/map.concatemer2pe.bam`. |
| `-t <threads>` | `16` | Number of threads. |
| `-x <preset>` | `map-ont` | minimap2 preset. Choices: `map-ont`, `map-hifi`. |
| `-q <mapq>` | `0` | MAPQ cutoff for `concatemer2pe.py`. |
| `-i <percent_identity>` | `0` | Percent identity cutoff. Valid range: `0-100`. |
| `-l <alignment_length>` | `0` | Alignment length cutoff. |
| `-h` | - | Show help and exit. |

Output:
```
<prefix>/map.concatemer2pe.bam
```

## `gphase popcnv`

Script: `pipeline/popCNV_pipeline.sh`

Usage:
```
gphase popcnv -f <fa_file> -p <output_prefix> -t <threads> -r <reads1.fastq.gz [reads2.fastq.gz ...]>
```

### Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `-f <fa_file>` | required | FASTA file containing unitig sequences. |
| `-p <output_prefix>` | required | Prefix for output files. |
| `-t <threads>` | `32` | Number of threads. |
| `-r <reads_files>` | required | One or more FASTQ/FASTQ.GZ read files. |
| `-h`, `--help` | - | Show help and exit. |

The file used as `gphase pipeline -c` is:
```
popcnv/06.genes.round.cn
```

## `gphase ppl`

Script: `pipeline/PPL_pipeline.sh`

Usage:
```
gphase ppl -g <genome.fa> -f <reads.fq> [options]
```

### Required parameters

| Parameter | Description |
| --- | --- |
| `-g <genome.fa>` | Genome FASTA file. |
| `-f <reads.fq>` | FASTQ file. |

### Optional parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `-j <jar_path>` | `bin/PPL-0.1.1.jar` | PPL jar file. |
| `-s <site>` | `^GATC` | Restriction enzyme site. |
| `-o <prefix>` | `PPL` | Output prefix. |
| `-t <threads>` | `12` | Number of threads. |
| `-q <mapq>` | `1` | MAPQ cutoff. |
| `-h` | - | Show help and exit. |

Output pairs file:
```
PPL/map.PPL.pairs
```
