import json, pandas as pd, matplotlib.pyplot as plt

# load UniProt JSON
with open("data/psbA_TM_uniprot.json") as f:
    data = json.load(f)
tm_uniprot = [(f["location"]["start"]["value"], f["location"]["end"]["value"]) 
              for f in data["features"] if f["type"]=="Transmembrane"]

# load DeepTMHMM GFF
tmhmm = []
with open("biolib_results/TMRs.gff3") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.split("\t")
        if parts[1]=="TMhelix":
            tmhmm.append((int(parts[2]), int(parts[3])))

# make dataframe
df = pd.DataFrame([
    *[(x[0],x[1],"UniProt") for x in tm_uniprot],
    *[(x[0],x[1],"DeepTMHMM") for x in tmhmm]
], columns=["start","end","source"])

# plot
fig, ax = plt.subplots(figsize=(8,2))
for i,(row) in enumerate(df.iterrows()):
    y = 0 if row[1]["source"]=="UniProt" else 1
    ax.plot([row[1]["start"], row[1]["end"]],[y,y],lw=8,
            label=row[1]["source"] if row[1]["source"] not in [l.get_label() for l in ax.lines] else None)
ax.set_yticks([0,1])
ax.set_yticklabels(["UniProt","DeepTMHMM"])
ax.set_xlabel("Residue position")
ax.legend()
plt.tight_layout()
plt.savefig("figures/compare_deepTMHMM_to_uniprot.png", dpi=300)
plt.show()
