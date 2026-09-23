<p align="center"><img src="assets/logo.png" width="180" alt="Smiling teal protein ribbon with a coral amino-acid bead"></p>

# Jev Protein Design 🧬

**A toy project that asks Jev to choose a protein sequence, one residue at a time.**

Give it a desired function in ordinary language. At each step, Jev sees the goal and the sequence so far, then chooses among the **20 standard amino acids plus STOP**. The program appends the chosen residue and repeats until Jev selects STOP or the length limit is reached.

> **Toy project, not a biological design tool.** Jev receives no experimental data, structure prediction, folding energy, or binding measurements. Its choices do not establish that a sequence folds, expresses, is safe, or performs the requested function. Treat the output as an illustrative candidate only; any real claim needs independent computational and laboratory validation.

## Quick start

Requires Python 3.10+ and a [TypeSafe AI Jev API key](https://typesafe.ai/). No Python dependencies are needed at runtime.

```bash
git clone https://github.com/JiaweiZhang1997/jev-protein-design.git
cd jev-protein-design
python3 -m pip install -e .
cp .env.example .env
```

Edit **`.env` locally** and replace the placeholder with your own `TYPESAFE_API_KEY`. The `.env` file is ignored by Git. You may instead set `TYPESAFE_API_KEY` in your environment. The key is only sent in the Authorization header to TypeSafe AI; it is never written to FASTA, traces, or the repository.

```bash
jev-protein-design \
  --function "a small, soluble protein with a stable fold" \
  --min-length 8 \
  --max-length 24 \
  --output designs/example.fasta \
  --trace designs/example.json
```

The command prints FASTA-like output. `--count 3` requests three independent runs (duplicates are possible). `--prefix M` starts with an existing sequence. `--help` lists all options. Every added residue requires one API call, so `--max-length` bounds the calls and cost per design.

### Small example cases

These prompts target **simple, countable sequence properties** so you can inspect Jev's decisions without pretending to verify a biological function:

```bash
# Does the chosen sequence tend to include more K/R residues when asked for basic composition?
jev-protein-design --function "a short peptide rich in basic K and R residues" --min-length 6 --max-length 12 --trace designs/basic.json

# Does the choice pattern change when the request asks for glycine-rich composition?
jev-protein-design --function "a short glycine-rich flexible peptide" --min-length 6 --max-length 12 --trace designs/glycine.json
```

Open the JSON traces to compare the goals, prefixes, chosen residues, and simple observed prefix counts. The prompts use property words as instructions to Jev; counting K/R or G in the output is a check of the toy workflow, not evidence of protein function.

Two recorded Jev runs are in [`examples/`](examples/). With a goal of a glycine-rich toy peptide, it produced `GGGGGGGGGGGGGGG`. With a goal of alternating glycine and serine without identical neighbors, it produced `GSGSGSGS`. The second trace shows the prefix changing before each choice. These runs illustrate that Jev can follow simple sequence composition instructions; they also show how easily a naive prompt can collapse into a repetitive sequence. Results may vary between model versions and runs.

### Search for similar annotated proteins

Add `--search` to submit the generated sequence to [NCBI BLASTP](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) against the curated Swiss-Prot database. NCBI asks API users to provide a contact email and limit request frequency, so pass your own address with `--ncbi-email`:

```bash
jev-protein-design \
  --function "a small soluble enzyme-like protein" \
  --min-length 25 --max-length 40 \
  --search --ncbi-email you@example.org \
  --trace designs/with-search.json
```

The search prints up to three similar annotated proteins with a title, accession link, E-value, identity, and query coverage. When available, it also retrieves the known protein's function annotation from [UniProtKB](https://www.uniprot.org/help/api). Check whether that annotation relates to your request, and whether the match spans much of your sequence. A short or low-complexity match can occur by chance; even a strong match does **not** show that the newly generated sequence has the same function. If the sequence is shorter than 15 residues, the program skips the search as too ambiguous. NCBI can take several minutes to return results. The generated sequence and your contact email go to NCBI; hit accessions go to UniProt. The Jev key stays local and is sent only to TypeSafe AI.

## How it works

```text
function request + current sequence
              │
              ▼
    Jev choice: A C D E F G H I K L M N P Q R S T V W Y STOP
              │
        append or finish
              │
              └────────────── repeat
```

Jev is a [typed decision model](https://typesafe.ai/blog/introducing-system-one-models-and-jev), not a sequence generator. This project turns sequence generation into a series of fixed-menu decisions. The state includes the desired function, current prefix, simple counts observed from that prefix, position, and length bounds. These counts are computed by code and change after every residue. The model's choice and the preceding state are recorded at each step; a JSON trace can be saved for inspection. Optionally, a BLAST similarity lookup follows generation; it is evidence for comparison, not a feedback signal used during generation.

STOP is always one of the 21 options. If Jev chooses it before `--min-length`, the program uses the highest-probability amino acid from that same response and continues. At `--max-length`, the program stops regardless of Jev's preference. There is no search over structures, feedback loop from experiments, or guarantee of diversity between runs.

## Install as a Codex Skill

The repository includes [`skills/jev-protein-design/SKILL.md`](skills/jev-protein-design/SKILL.md). After cloning, copy that folder into your Codex skills directory:

```bash
mkdir -p ~/.codex/skills
cp -R skills/jev-protein-design ~/.codex/skills/
```

Then ask Codex to use **`$jev-protein-design`** for an exploratory protein sequence. The Skill guides setup and runs this CLI. Each user supplies their own Jev key locally; the Skill contains no credentials.

## Development

```bash
python3 -m unittest discover -s tests -v
```

The API call uses TypeSafe AI's [`/v1/systemone` choice format](https://api.typesafe.ai/docs). This repository is a small educational experiment, not a validated method for protein engineering.

### Related work

As of September 2026, I did not find a public project using **Jev specifically** to pick protein residues one at a time. The closest idea in Jev's own examples is the [Wikiracing decision loop](https://typesafe.ai/blog/introducing-system-one-models-and-jev): the model repeatedly picks one option from a defined set based on the changing state. For actual protein sequence generation from text, [ProteinDT](https://github.com/chao1224/ProteinDT) and [ProtDAT](https://github.com/GXY0116/ProtDAT) are specialist research projects with different methods and substantially stronger domain grounding. This project is an exploration of Jev's decision interface, not a replacement for them.

## 中文简介

这是一个**玩具项目**：把“生成蛋白质序列”拆成逐位选择，每次由 Jev 在 20 种常规氨基酸和 `STOP` 中选一个。输入的功能描述只是提示，不代表输出序列真的具有该功能。真实用途必须另做结构、表达、功能和安全性验证。API key 只保存在使用者本机的 `.env` 或环境变量中。

## License

MIT
