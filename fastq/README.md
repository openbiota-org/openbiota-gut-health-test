# Your sequencing files go here

Make one folder per sample and put that sample's two FASTQ files in it:

```
fastq/
└── my-sample/
    ├── my-sample_1.fastq.gz
    └── my-sample_2.fastq.gz
```

The files from a sequencing provider are usually named `SOMETHING_1.fastq.gz`
and `SOMETHING_2.fastq.gz`, or `SOMETHING_R1.fastq.gz` and
`SOMETHING_R2.fastq.gz`. Keep those names: the sample is named by the file,
so the two above are the sample `my-sample`. The folder's name is up to you.

Then, from the repository folder:

```
bin/openbiota run
```

The report is written to `results/my-sample/my-sample_report.pdf`.

Nothing you put in this folder is ever added to git: the whole folder is
ignored apart from this note.
