# Input formats

## TXT counts

Each block starts with `Time Stamp:` and contains `Elapsed: <hour> hours` followed by eight tab-separated rows A–H, each with 12 numeric values. See the deterministic writer in `demo.synthetic_records()` for a complete working fixture. Exactly 0,2,...,24 h and all 96 wells are required. The original text is stored with the experiment.

## PlateMap worksheet

One row per well, 96 unique wells, A1–H12 (zero-padded well names also accepted).

| Column | Values |
|---|---|
| plate_id | A single physical-plate identifier; can be blank for isolated file examples |
| well | A1–H12 |
| cell_line | EV or R273H |
| treatment | Treatment identifier/name; DMSO and Cells only have dedicated roles |
| replicate | Technical replicate label |
| treatment_class | Screen drug, Named drug, or Control/reference |
| dmso_reference | Excel boolean true for DMSO reference wells |
| include_in_plots | Optional source metadata; false does not exclude counts from analysis |
| notes | Optional source notes |

RSL3 retains its positive-control identity and may also be included in drug rankings. Screen drug wells default to the established 5 μM design, RSL3 to 1 μM; other named-treatment doses are not invented. Each cell type must have its own DMSO controls. Only PlateMap is read; any Data worksheet is ignored.

## Chinese library layout

The worksheet 实验板逐孔信息 requires 样本板号, 细胞, 样本孔位, 孔类型, MCE编号, 药物名称. Optional fields include 原库板号, 原库孔位, 复孔编号. E-suffixed plates contain EV and H-suffixed plates contain H/R273H. Matching prefixes and drug identifiers form a pair.

Local data, user layouts, reports, settings and databases are excluded by .gitignore. Check staged files before sharing any customized checkout.
