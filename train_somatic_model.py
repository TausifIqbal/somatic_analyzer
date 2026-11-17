import cyvcf2
import pandas as pd
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    roc_auc_score,
        )

from xgboost import XGBClassifier

truth_vcf_path = "data/truth_somatic_snv.vcf.gz"
truth_variants = set()

for v in cyvcf2.VCF(truth_vcf_path):
    truth_variants.add((v.CHROM, v.POS, v.REF, v.ALT[0]))

print(f"Loaded {len(truth_variants):,} ground-truth somatic SNVs.")

cand_vcf_path = "data/filtered_candidates.vcf.gz"
cand_vcf = cyvcf2.VCF(cand_vcf_path)

print("VCF Sample Columns:", cand_vcf.samples)

tumor_idx = cand_vcf.samples.index("TUMOR")
norm_idx = cand_vcf.samples.index("NORMAL")

print(f"TUMOR index: {tumor_idx} | NORMAL index: {norm_idx}")

records =[]

for v in cand_vcf:
    # keep biallelic SNVs only
    if len(v.REF) != 1 or len(v.ALT) != 1 or len(v.ALT[0]) != 1: 
        continue

    # Extract AF (Allele Frequency) and AD (Allelic Depth: [ref_count, alt_count])
    af = v.format("AF")
    ad = v.format("AD")
    if af is None or ad is None:
        continue

    tumor_af = float(af[tumor_idx][0])
    normal_af = float(af[norm_idx][0])

    tumor_ref_dp, tumor_alt_dp =  int(ad[tumor_idx][0]), int(ad[tumor_idx][1])
    normal_ref_dp, normal_alt_dp = int(ad[norm_idx][0]), int(ad[norm_idx][1])

    tumor_dp = tumor_ref_dp + tumor_alt_dp
    normal_dp = normal_ref_dp + normal_alt_dp

    def get_info_float(key):
        val = v.INFO.get(key)
        if val is None:
            return np.nan
        return float(val[0]) if isinstance(val, tuple) else float(val)

    tlod = get_info_float("TLOD") 
    nlod = get_info_float("NLOD")
    ecnt = get_info_float("ECNT")
    hcnt = get_info_float("HCNT")

    is_somatic = 1 if (v.CHROM, v.POS, v.REF, v.ALT[0]) in truth_variants else 0

    records.append({
        "CHROM": v.CHROM,
        "POS": v.POS,
        "tumor_af": tumor_af,
        "normal_af": normal_af,
        "af_diff": tumor_af - normal_af,
        "af_ratio": (tumor_af + 1e-4) / (normal_af + 1e-4),
        "tumor_dp": tumor_dp,
        "normal_dp": normal_dp,
        "tumor_alt_dp": tumor_alt_dp,
        "normal_alt_dp": normal_alt_dp,
        "dp_ratio": (tumor_dp+1)/(normal_dp +1),
        "tlod": tlod,
        "nlod": nlod,
        "ecnt": ecnt,
        "hcnt": hcnt,
        "label": is_somatic,
    })
df = pd.DataFrame(records)

print("\nDataset shape:", df.shape)
print("Class distribution (1 = Somatic, 0 = Non-Somatic):\n", df["label"].value_counts())

feature_cols = [ 
                "tumor_af", "normal_af", "af_diff", "af_ratio",
                "tumor_dp", "normal_dp", "dp_ratio",
                "tlod", "nlod", "ecnt", "hcnt"
                ]

test_chroms = {"chr1", "chr2"}
train_df = df[~df["CHROM"].isin(test_chroms)]
test_df = df[df["CHROM"].isin(test_chroms)]

X_train, y_train = train_df[feature_cols], train_df["label"]
X_test, y_test = test_df[feature_cols], test_df["label"]

print("f\nTraining variants: {len(X_train):,} | Testing Varints (chr1-2): {len(X_test):,}")


model = XGBClassifier(
        n_estimators=300,
        max_depth=6,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        random_state=42,
        n_jobs=-1,
        eval_metric="logloss"
        )
model.fit(X_train, y_train)

y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:,1]

tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)
precision = tp / (tp + fp)
f1 = 2 * (precision * sensitivity) / (precision + sensitivity)
auroc = roc_auc_score(y_test, y_prob)

print("\n=== Test Set Performance (Held-out chr1 & chr2) ===")
print(f"True Positives (TP):  {tp:,}")
print(f"False Positives (FP): {fp:,}")
print(f"True Negatives (TN):  {tn:,}")
print(f"False Negatives (FN): {fn:,}")
print("-" * 45)
print(f"Sensitivity (Recall): {sensitivity:.4f}")
print(f"Specificity:          {specificity:.4f}")
print(f"Precision:            {precision:.4f}")
print(f"F1 Score:             {f1:.4f}")
print(f"ROC-AUC:              {auroc:.4f}")
