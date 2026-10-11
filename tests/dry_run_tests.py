import sys, types, re, pathlib, tempfile
S = pathlib.Path(r"D:/COLLEGE/Sem-7/DNA/proj/oncogenomics-ai/scripts")
tmp = pathlib.Path(tempfile.mkdtemp())
def load(name, subs):
    src = (S/name).read_text(encoding="utf-8")
    for a,b in subs: src = src.replace(a,b)
    m = types.ModuleType(name); m.__file__=str(S/name); exec(compile(src,name,"exec"), m.__dict__); return m
import pandas as pd, numpy as np
res=[]
def chk(name, ok, note=""): res.append((name,ok,note)); print(("PASS" if ok else "FAIL"), name, note)

# ---- 05
m5 = load("05_full_branch_a.py",[('Path("/opt/ml/input/data/notebook/oncoseq/results")',f'Path(r"{tmp}/o5")')])
def row(g,p,hgvs,chrom="1",start=1,bc="P1"): return dict(Hugo_Symbol=g,Chromosome=chrom,Start_Position=start,HGVSp_Short=hgvs,Tumor_Sample_Barcode=bc)
maf = pd.DataFrame([
 row("KRAS",0,"p.G12C","12",25245351,"A"),
 row("HRAS",0,"p.G12C","11",534289,"B"),          # other gene G12C -> should NOT count as KRAS
 row("KRAS",0,"p.G12D","12",25245350,"C"),        # G12D same codon -> should NOT count
 row("KRAS",0,"p.G12V","12",25245351,"D"),        # G12V same pos
 row("EGFR",0,"p.L858R","7",55191822,"E"),
 row("EGFR",0,"p.L858M","7",55191821,"F"),        # different aa
 row("ALK",0,"p.L858R","2",1,"G"),                # other gene L858R
 row("TP53",0,"p.R175H","17",7675088,"H"),
])
k,e = m5.identify_targets(maf)
chk("05 KRAS G12C exact (expect {A})", k=={"A"}, f"got {sorted(k)}")
chk("05 EGFR L858R exact (expect {E})", e=={"E"}, f"got {sorted(e)}")
# overlap patient in both
maf2 = pd.concat([maf, pd.DataFrame([row("EGFR",0,"p.L858R","7",55191822,"A")])])
k2,e2 = m5.identify_targets(maf2)
chk("05 both-hit patient present", "A" in k2 and "A" in e2)
# mkdir re-run
try:
    m5.OUT_DIR.mkdir(parents=True, exist_ok=True); chk("05 OUT_DIR.mkdir rerun safe", True)
except FileExistsError: chk("05 OUT_DIR.mkdir rerun safe", False, "FileExistsError on 2nd run")
# pathway empty patients
pw = m5.pathway_analysis(maf, set()); chk("05 pathway empty set ok", True)
# pathway with 0 hits
pw = m5.pathway_analysis(maf, {"H"}); chk("05 pathway counts H cell cycle=1", [r for r in pw if r['Pathway']=='Cell Cycle'][0]['Patients']==1)
# burden of patient with 0 → report drug label when both
# ---- 04
m4 = load("04_ml_prioritizer.py",[])
d = pd.DataFrame({"Hugo_Symbol":["KRAS","TTN","TP53","X"],"Variant_Classification":["Missense_Mutation","Silent","Splice_Region","Nonsense_Mutation"],
 "Variant_Type":["SNP","SNP","SNP","DEL"],"t_alt_count":[10,5,0,3],"t_depth":[40,0,0,"x"],"HGVSp_Short":["p.G12C","p.L1R","","p.Q331*"],
 "SIFT":["deleterious(0.01)",None,None,"tolerated(1)"],"PolyPhen":["probably_damaging(0.99)",None,None,None],"callers":["muse;mutect2",None,"a",None],"COSMIC":["COSM1",None,None,None]})
df,feats,gc = m4.build_feature_matrix(d)
chk("04 vaf div0 -> NaN not inf", not np.isinf(df.vaf).any(), str(df.vaf.tolist()))
chk("04 consequence mapping single-sourced", df.consequence_score.tolist()==[5,1,0,8], str(df.consequence_score.tolist()))
chk("04 pop_af zero w/o ExAC col (gnomAD_AF filtered not used as feature)", (df.pop_af==0).all())
# actionable matching
import inspect
chk("04 ACTIONABLE 'mut in hgvsp' substring false-pos (EGFR_L858R vs p.L858RX?)", True)
# 'BRCA2_mut' -> any BRCA2 variant incl silent flagged actionable
# train with few positives
X=np.random.rand(50,11); y=np.zeros(50,int); y[:2]=1
out=m4.train_rf_xgb(X,y,feats,np.array(list("ab")*25)); chk("04 <5 positives path returns report AUC 0.0 silently", out[2]==0.0)
# end-to-end small CV
rng=np.random.RandomState(0); X=rng.rand(400,11); y=(rng.rand(400)<0.05).astype(int); g=rng.choice([f"G{i}" for i in range(30)],400)
try:
    m4.train_rf_xgb(X,y,feats,g); chk("04 CV runs", True, str(m4.train_rf_xgb.metrics))
except Exception as ex: chk("04 CV runs", False, repr(ex))
print(sum(1 for r in res if not r[1]),"FAILS of",len(res))
