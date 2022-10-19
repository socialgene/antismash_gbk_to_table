# antismash_gbk_to_table

## Installation

```
pip install https://github.com/socialgene/antismash_gbk_to_table/archive/main.tar.gz
```

## Usage

```
wget https://antismash-db.secondarymetabolites.org/output/GCF_004339725.1/GCF_004339725.1.zip
mkdir GCF_004339725
unzip -n -d GCF_004339725 GCF_004339725.1.zip

antismash_gbk_to_table --input ./GCF_004339725/*.region*.gbk --output GCF_004339725_table.tsv --header true
```

## Results

```
cat GCF_004339725_table.tsv
```

| locus           | tool      | type         | start   | end     | strand | category | contig_edge | neighbouring | product                | protocluster_number | candidate_cluster_number | candidate_cluster_numbers | region_number | kind            |
| --------------- | --------- | ------------ | ------- | ------- | ------ | -------- | ----------- | ------------ | ---------------------- | ------------------- | ------------------------ | ------------------------- | ------------- | --------------- |
| NZ_SMGK01000001 | antismash | protocluster | 22733   | 45843   | 1      |          | False       |              | LAP                    | 1                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | proto_core   | 32733   | 35843   | 1      |          |             |              | LAP                    | 1                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | cand_cluster | 22733   | 45843   | 1      |          | False       |              | LAP                    |                     | 1                        |                           |               | single          |
| NZ_SMGK01000001 | antismash | region       | 22733   | 45843   | 1      |          | False       |              | LAP                    |                     |                          | 1                         | 1             |                 |
| NZ_SMGK01000001 | antismash | cand_cluster | 558000  | 645023  | 1      |          | False       |              | NRPS; T1PKS; NRPS-like |                     | 1                        |                           |               | neighbouring    |
| NZ_SMGK01000001 | antismash | region       | 558000  | 645023  | 1      |          | False       |              | NRPS; T1PKS; NRPS-like |                     |                          | 1; 2; 3                   | 2             |                 |
| NZ_SMGK01000001 | antismash | protocluster | 558000  | 610087  | 1      |          | False       |              | NRPS                   | 1                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | proto_core   | 578000  | 590087  | 1      |          |             |              | NRPS                   | 1                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | cand_cluster | 558000  | 610087  | 1      |          | False       |              | NRPS; T1PKS            |                     | 2                        |                           |               | chemical_hybrid |
| NZ_SMGK01000001 | antismash | protocluster | 559358  | 607310  | 1      |          | False       |              | T1PKS                  | 2                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | proto_core   | 579358  | 587310  | -1     |          |             |              | T1PKS                  | 2                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | protocluster | 603188  | 645023  | 1      |          | False       |              | NRPS-like              | 3                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | proto_core   | 623188  | 625023  | 1      |          |             |              | NRPS-like              | 3                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | cand_cluster | 603188  | 645023  | 1      |          | False       |              | NRPS-like              |                     | 3                        |                           |               | single          |
| NZ_SMGK01000001 | antismash | cand_cluster | 1172805 | 1235383 | 1      |          | False       |              | T3PKS; terpene         |                     | 1                        |                           |               | neighbouring    |
| NZ_SMGK01000001 | antismash | region       | 1172805 | 1235383 | 1      |          | False       |              | T3PKS; terpene         |                     |                          | 1; 2; 3                   | 3             |                 |
| NZ_SMGK01000001 | antismash | protocluster | 1172805 | 1213854 | 1      |          | False       |              | T3PKS                  | 1                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | proto_core   | 1192805 | 1193854 | 1      |          |             |              | T3PKS                  | 1                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | cand_cluster | 1172805 | 1213854 | 1      |          | False       |              | T3PKS                  |                     | 2                        |                           |               | single          |
| NZ_SMGK01000001 | antismash | protocluster | 1213380 | 1235383 | 1      |          | False       |              | terpene                | 2                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | proto_core   | 1223380 | 1225383 | 1      |          |             |              | terpene                | 2                   |                          |                           |               |                 |
| NZ_SMGK01000001 | antismash | cand_cluster | 1213380 | 1235383 | 1      |          | False       |              | terpene                |                     | 3                        |                           |               | single          |
| NZ_SMGK01000002 | antismash | protocluster | 34229   | 76787   | 1      |          | False       |              | NRPS-like              | 1                   |                          |                           |               |                 |
| NZ_SMGK01000002 | antismash | proto_core   | 54229   | 56787   | -1     |          |             |              | NRPS-like              | 1                   |                          |                           |               |                 |
| NZ_SMGK01000002 | antismash | cand_cluster | 34229   | 76787   | 1      |          | False       |              | NRPS-like              |                     | 1                        |                           |               | single          |
| NZ_SMGK01000002 | antismash | region       | 34229   | 76787   | 1      |          | False       |              | NRPS-like              |                     |                          | 1                         | 1             |                 |
| NZ_SMGK01000002 | antismash | protocluster | 640246  | 662158  | 1      |          | False       |              | terpene                | 1                   |                          |                           |               |                 |
| NZ_SMGK01000002 | antismash | proto_core   | 650246  | 652158  | 1      |          |             |              | terpene                | 1                   |                          |                           |               |                 |
| NZ_SMGK01000002 | antismash | cand_cluster | 640246  | 662158  | 1      |          | False       |              | terpene                |                     | 1                        |                           |               | single          |
| NZ_SMGK01000002 | antismash | region       | 640246  | 662158  | 1      |          | False       |              | terpene                |                     |                          | 1                         | 2             |                 |
| NZ_SMGK01000008 | antismash | protocluster | 82608   | 93489   | 1      |          | False       |              | bacteriocin            | 1                   |                          |                           |               |                 |
| NZ_SMGK01000008 | antismash | proto_core   | 87608   | 88489   | 1      |          |             |              | bacteriocin            | 1                   |                          |                           |               |                 |
| NZ_SMGK01000008 | antismash | cand_cluster | 82608   | 93489   | 1      |          | False       |              | bacteriocin            |                     | 1                        |                           |               | single          |
| NZ_SMGK01000008 | antismash | region       | 82608   | 93489   | 1      |          | False       |              | bacteriocin            |                     |                          | 1                         | 1             |                 |
