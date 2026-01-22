import time, wandb, os, shutil

wandb.init(project="raxml", name="search1tree-log")

sources = {
    "parsimony_w_order": "raxml_input/aa_treeConstrained_fiveParsimony.raxml.log",
}

while True:
    for tag, src in sources.items():
        if os.path.exists(src):
            dst = f"{tag}.raxml.log"
            shutil.copyfile(src, dst)
            wandb.save(dst, policy="now")
    time.sleep(30)
