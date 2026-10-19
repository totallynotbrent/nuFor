# post-build gate: catch broken math rendering in the built pages.
# two failure signatures, both seen in real breakage:
# 1. literal latex still visible: raw $..$ spans or \commands in text
# 2. katex fallback: commands flattened letter-by-letter into plain
#    words ("fracpartial rhoYs"), no backslashes - pattern 1 misses it
import os
import re
import sys

pub = os.path.join(os.path.dirname(os.path.abspath(__file__)), "public")
raw_cmd = re.compile(r"\\(?:frac|partial|nabla|sum|rho|omega|eta|gamma|mathbf|qquad|sqrt|tau|Delta|dot|left|right)\b")
# multi-command fusions only: single words like "nabla" false-positive
# against real prose ("run a blast march"), so require command pairs that
# never occur as english words
fallback_words = (
    "fracpartial",
    "partialrho",
    "nablacdot",
    "mathbf{",
    "frac{partial",
    "frac{d",
    "frac{\\",
)
bad = []
for root, _, files in os.walk(pub):
    for f in files:
        if not f.endswith(".html"):
            continue
        path = os.path.join(root, f)
        src = open(path, encoding="utf-8", errors="ignore").read()
        body = re.sub(r"<script.*?</script>|<style.*?</style>", "", src, flags=re.S)
        text = re.sub(r"<[^>]+>", " ", body)
        # raw signature: match on the tagged text (backslashes are unambiguous)
        raw_hits = sorted(set(raw_cmd.findall(text)))
        # fallback signature: katex flattens words AND strips spaces,
        # so squeeze whitespace and look for command fusions
        squeezed = re.sub(r"\s+", "", text)
        word_hits = [w for w in fallback_words if w in squeezed]
        if raw_hits or word_hits:
            bad.append((os.path.relpath(path, pub), raw_hits[:4], word_hits[:4]))

if bad:
    for page, raws, words in bad:
        print(f"RAW LATEX in {page}: commands={raws} fallback_words={words}")
    sys.exit(1)
print("rendered-text gate: no raw or fallback latex in any built page")
