#kmer_freq

Program that produce the most frequent DNA-k-mers of arbitrary length sorted by frequency in a given FASTQ file.

###usage
    kmer_freq --filename <path> --kmersize <size> --topcount <count>
    
    test run with `make run`
    `PAGING` flag can be switched for page-mapped file reads using POSIX mmap().


###roadmap
- switch to batch'd realloc() line buffer for page-mapped FASTQ file read.
- determistic offset skip for page-mapped FASTQ files.

