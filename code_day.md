code_day
========================================================
author: Ramón Gallego
date: 20180827
autosize: true

What does this pipeline do?
========================================================

Once you have retrieved all your sequences from the Illumina machine, you will need to asign them to the sample they came from,  and minimize the noise and cross-contamination inherent to this platform.

This pipeline will need the follwing **input** files:

- One or more pairs of fastq files
- A metadata file, with one row per sample, containing the information needed to demultiplex the sequences

```
# A tibble: 93 x 6
   sample_id Tag    pri_index_name sec_index_seq file1        file2       
   <chr>     <chr>  <chr>          <chr>         <chr>        <chr>       
 1 SA07A.1   Tag_12 Lib_B          GATGAC        Lib-B_S2_L0… Lib-B_S2_L0…
 2 SA07B.1   Tag_12 Lib_D          GATGAC        Lib-D_S4_L0… Lib-D_S4_L0…
 3 SA07C.1   Tag_4  Lib_B          GCGCTC        Lib-B_S2_L0… Lib-B_S2_L0…
 4 TR07A.1   Tag_8  Lib_G          CTCGCA        Lib-G_S7_L0… Lib-G_S7_L0…
 5 TR07B.1   Tag_2  Lib_A          ACAGCA        Lib-A_S1_L0… Lib-A_S1_L0…
 6 TR07C.1   Tag_10 Lib_B          TGTATG        Lib-B_S2_L0… Lib-B_S2_L0…
 7 LL07A.1   Tag_7  Lib_G          ATATCG        Lib-G_S7_L0… Lib-G_S7_L0…
 8 LL07B.1   Tag_1  Lib_A          ACACAC        Lib-A_S1_L0… Lib-A_S1_L0…
 9 LL07C.1   Tag_5  Lib_G          ACATGT        Lib-G_S7_L0… Lib-G_S7_L0…
10 PO07A.1   Tag_9  Lib_A          TCGCAT        Lib-A_S1_L0… Lib-A_S1_L0…
# ... with 83 more rows
```
- A parameters file, which will be the bridge between the pipeline and the metadata file 

What does this pipeline give you?
========================================================

In return, this pipeline will give you the composition of your samples. The **output** files are:

- ASV_Table: It is in a long format, with three columns: Hash, sample, nReads


```
# A tibble: 5 x 3
  sample  Hash                                     nReads
  <chr>   <chr>                                     <int>
1 LL08B.3 a66d2cafe17c8203af1b9a46ae27ba84707bcd20      8
2 LL08C.2 9c6622bc6a7f13f9592c93c84d5576b438579098      5
3 SA08A.2 ad2a209d5762efa8de63ad8ba6b53b4a337341de      2
4 LL08C.1 ad2a209d5762efa8de63ad8ba6b53b4a337341de      3
5 LL08B.3 c84c078afb75141b53650e01b2f1cf1e1cb8c3bc      5
```

- Hash_key: In case you used hashing, the key to convert sequences to hashes


```
# A tibble: 4 x 2
  Hash                                     Sequence                       
  <chr>                                    <chr>                          
1 a66d2cafe17c8203af1b9a46ae27ba84707bcd20 AGAACTTGGTAATCCTTGATATTGTGGTAG…
2 9c6622bc6a7f13f9592c93c84d5576b438579098 AAAAGGTAATAGTATGGTGTTAGTGGCTAT…
3 ad2a209d5762efa8de63ad8ba6b53b4a337341de AAAAAGGTAATAGTATGGTGTTAGTGGCTA…
4 c84c078afb75141b53650e01b2f1cf1e1cb8c3bc ACTGTTTATCCTCCAGAACTTGGTAATCCT…
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-4](code_day-figure/unnamed-chunk-4-1.png)