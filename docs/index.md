<p class="eyebrow">Documentation</p>

# OpenBiota Gut Health Test

**Raw shotgun sequencing reads in. A gut-health report you can actually read out.**

The OpenBiota software (`openbiota` on the command line) takes the paired-end
FASTQ files behind a shotgun stool metagenome and turns them — every read, no
subsampling — into the OpenBiota Gut Health Test: a gut-health index,
species-level community composition, copies-per-100-genomes capacity for 25
metabolite pathways, resemblance scores against 33 published disease-associated
microbiome patterns, a 485-target pathogen screen and an estimated microbiome
age. Page one is the whole answer in graphics; the pages behind it are the
evidence, and every number traces back to the reads that produced it.

<figure markdown>
  ![Page 1 of the report](images/summary-page.png){ width="720" }
  <figcaption>Page 1: the gut-health dial, headline findings, metabolite pathways on one shared scale, the community donut and disease-pattern resemblance.</figcaption>
</figure>

!!! note "Research use only"
    The report measures *genetic capacity* and *community resemblance*. It does
    not measure metabolites and it does not diagnose disease. A match to a
    disease signature is a clue for further investigation, not a diagnosis.

## Where to start

<div class="grid cards" markdown>

-   **[Quickstart](QUICKSTART.md)**

    From a fresh clone to your first PDF: tools, reference data, one command.

-   **[Installation](INSTALL.md)**

    macOS and Linux, the two MetaPhlAn environments, verifying the install.

-   **[Usage](USAGE.md)**

    Every `openbiota` command and flag, caching, performance tuning.

-   **[Reading the report](OUTPUT.md)**

    What each section means, the seven-band scale, `results.json`.

-   **[Method](METHOD.md)**

    How the measurements are made and why they are made that way.

-   **[Validation](VALIDATION.md)**

    Sensitivity, specificity, the reference cohort, and a like-for-like
    comparison against a commercial report.

-   **[Organism detection](EXPANDED_DETECTION.md)**

    Nine detection methods over the largest current catalogues, merged into
    one organism list with one share of the whole per organism, competitive
    confirmation and strain placement.

-   **[Measured capability](MEASURED_CAPABILITY.md)**

    What the detection found on simulated communities with known composition,
    against the targets it was built to; stated as measured, not as claimed.

</div>

## The project

OpenBiota is free for any noncommercial purpose under the
[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0);
commercial use requires a licence from OpenBiota. The website is
<a href="../index.html">openbiota.com</a>; the code lives on
[GitHub](https://github.com/openbiota-org/openbiota-gut-health-test). Add a
published study, improve an analysis, or make a finding easier to understand —
every contribution brings more of the world's microbiome research into a report
anyone can generate.
