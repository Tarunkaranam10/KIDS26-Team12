# Deterministic 450K-to-EPIC feature bridge

## Frozen candidate mask

`scripts/build_probe_bridge.py` constructs the candidate feature space from Zhou Lab InfiniumAnnotationV1 archived 2022-09 hg19 manifests. The current run produced:

| Stage | CpG probes |
|---|---:|
| HM450 manifest | 482,421 |
| EPIC-v1 manifest | 862,927 |
| Probe IDs shared by both | 449,521 |
| Shared autosomal probes with identical hg19 CpG coordinate and `MASK_general=FALSE` on both arrays | 384,640 |

Output: `config/shared_autosomal_probes.txt`

SHA-256: `f359e43bc171c544e55e60000905bb0bba1c44ee0b37f3144c8ec733e052e5fc`

The metadata sidecar records URLs, manifest hashes, counts, creation time, and rules. Downloaded annotation inputs remain under ignored `data/raw/annotation/`.

## Rules and rationale

1. Intersect exact `cg` probe identifiers in HM450 and the intended EPIC manifest.
2. Retain chromosomes 1-22. Sex-chromosome removal limits sex/lineage shortcuts and follows major cross-platform classifier precedents; sex-specific analysis would require a separate declared model.
3. Require the same hg19 CpG coordinate in both manifests. A mismatched mapping is excluded rather than guessed.
4. Require `MASK_general=FALSE` in both manifests. This versioned composite mask captures manifest-annotated mapping, cross-reactivity, SNP/confounding, and reliability exclusions. The exact underlying mask columns remain in the source manifests and sidecar provenance.
5. Do not perform outcome-based selection here. Training-fold variance filtering occurs inside each model-training fold.
6. For a sample, preserve missingness; training-fold imputation parameters are applied at inference. Record observed-probe coverage and abstain below a frozen coverage/OOD gate.

This mask avoids multiple uncoordinated blacklist versions and excessive cumulative filtering. It is a technical candidate mask, not a guarantee that every probe is biologically transportable.

## Required freeze checks

- Confirm whether PBTP used EPIC v1 or EPIC v2. The generated list is for **EPIC v1** and must not be reused for EPIC v2.
- Confirm PBTP probe identifiers and genome-build metadata against the exact array manifest.
- Compare beta distributions and missingness on public samples processed through the intended frozen pipeline.
- Lock manifest URLs/hashes, output hash, and software versions before viewing PBTP outcomes.
- Keep raw genotype-sensitive/SNP probes outside ordinary methylation predictors. Any experimental BAF project requires raw channels, allele orientation, heterozygosity evidence, and separate genomic validation.

## Reproduction

```sh
python scripts/build_probe_bridge.py
```

The script is deterministic for the pinned remote bytes and refuses silent coordinate/mask ambiguity. If the annotation source changes bytes, the new hashes make the change visible and require a new feature-space version.
