# red

This is a github fork of RED (REpeat Detector) by Hani Z. Girgis for investigation of possible inclusion in Cactus's preprocessor (ideally as a replacement for lastz).

## Homepage

[http://toolsmith.ens.utulsa.edu/](http://toolsmith.ens.utulsa.edu/)

## Manuscript

[https://doi.org/10.1186/s12859-015-0654-5](https://doi.org/10.1186/s12859-015-0654-5)

## Version used to start this fork

[http://toolsmith.ens.utulsa.edu/red/data/DataSet1Src.tar.gz](http://toolsmith.ens.utulsa.edu/red/data/DataSet1Src.tar.gz) with md5sum 37574e07eb8d012cb6975e7163e9b62b

## Help Screen

```
This is Red (REpeat Detector) designed and developed by Hani Zakaria Girgis, PhD.

Version: 05/22/2015

Argument pairs of the form: -flag value are required.
Valid argument pairs:
        -gnm input genome directory, required.
                Files with ".fa" extension in this directory are used for completing the table of the adjusted counts.
                These Files are scanned for repeats.
        -dir directory including additional input sequences, optional.
                Files with ".fa" extension in this directory are NOT used for completing the table.
                These Files MUST have different names from those in the genome directory.
                These Files are scanned for repeats.
        -len word length equals k defining the k-mer. The default is floor(log_4(genome size)).
        -ord order of the background Markov chain. The default is floor(k/2)-1.
        -gau half width of the mask. The default is based on the GC content.
                20 if the GC content > 33% and < 67%, 40 otherwise.
        -thr the threshold score of the low adjusted scores of non-repeats. The default is 2.
        -min the minimum number of the observed k-mers. The default is 3.
        -tbl file where the table of the adjusted counts is written, optional.
        -sco directory where scores are saved, optional.
                Score files have the ".scr" extension.
        -cnd directory where candidate regions are saved, optional.
                Candidates files have the ".cnd" extension.
        -rpt directory where repeats locations are saved, optional.
                Repeats files have the ".rpt" extension.
        -msk directory where masked sequences are saved, optional.
                Masked sequences files have the ".msk" extension.
        -frm the format of the output: 1 (chrName:start-end) or 2 (chrName      start   end).
                The output format are zero based and the end is exclusive.
        -hmo file where the HMM is saved, optional.

Examples:
        The following command runs Red with the defaults and generates the masked sequences.
        Red -gnm genome_directory -msk output_directory

        The following command runs Red with the defaults and generates the masked sequences and the locations of repeats.
        Red -gnm genome_directory -msk output_directory -rpt output_directory
```

