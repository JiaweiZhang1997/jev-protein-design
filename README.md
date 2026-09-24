<p align="center"><img src="assets/logo.png" width="180" alt="Smiling teal protein ribbon with a coral amino-acid bead"></p>

# Jev Protein Design 🧬

**A toy project that turns a desired protein function into a sequence of Jev decisions.**

Describe a function in ordinary language. Jev chooses one of the **20 standard amino acids plus STOP**, appends it to the sequence, and repeats.

> This is an educational experiment. Generated sequences have no demonstrated biological function; sequence similarity and predicted structures are exploratory checks.

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

The initial goal includes a fixed constraint to avoid long runs of one residue and repeated short motifs. **Only `current_sequence` changes between requests.** The goal, instructions, choices, and length bounds stay fixed. Jev receives no structural or experimental feedback during generation.

Generation ends at STOP or the maximum length. If Jev selects STOP before the minimum length, the program uses its highest-probability amino acid from the same response. A JSON trace records each choice and its preceding sequence.

## Quick start

Requires Python 3.10+ and a [TypeSafe AI Jev API key](https://typesafe.ai/).

```bash
git clone https://github.com/JiaweiZhang1997/jev-protein-design.git
cd jev-protein-design
python3 -m pip install -e .
cp .env.example .env
```

Add your own `TYPESAFE_API_KEY` to **`.env` locally**, or set it in your environment. The `.env` file is ignored by Git, and keys are never included in generated files.

```bash
jev-protein-design \
  --function "a small, soluble protein with a stable fold" \
  --min-length 8 \
  --max-length 24 \
  --output designs/example.fasta \
  --trace designs/example.json
```

Use `--count 3` for three runs, `--prefix M` to extend an existing sequence, or `--help` for all options. Each added residue requires one Jev API call.

## Examples: GFP and RuBisCO

Two full-length candidates were generated from functional goals, then predicted with **Boltz-2.1** and compared with experimental protein structures:

- **Green fluorescent protein, 238 aa:** [GFP sequence P42212](https://www.uniprot.org/uniprotkb/P42212/entry), [PDB 1GFL chain A](https://www.rcsb.org/structure/1GFL).
- **RuBisCO large chain, 475 aa:** [spinach RuBisCO sequence P00875](https://www.uniprot.org/uniprotkb/P00875/entry), [PDB 8QJ0 chain L](https://www.rcsb.org/structure/8QJ0).

![Cartoon views of Jev candidates predicted by Boltz beside experimental PDB protein structures](assets/function-comparison.png)

| Requested function | Length | Sequence identity | TM-score vs PDB |
| --- | ---: | ---: | ---: |
| GFP | 238 aa | 13.4% | 0.299 |
| RuBisCO large chain | 475 aa | 21.7% | 0.259 |

Sequence identity counts identical aligned residues relative to the longer full sequence. TM-scores use the experimental chain length. See [comparison data](examples/function-comparison.json) for alignment details.

The GFP candidate has `EEK` at positions 65–67, where the reference has the chromophore-forming `SYG`. The RuBisCO candidate still contains repeated motifs and an eight-residue identical run despite the fixed anti-repetition goal. **These results do not establish fluorescence or carbon fixation.**

<details>
<summary>GFP: generated and reference sequences</summary>

```fasta
>gfp-jev
MLIVILIFVIVLAVSATITIPISPSIASTISIPSISGPIPAPIVIIIISVAVEAAVAVLVLLLVEEKRNNLALLLPLPDD
SSSNGYYGYYYWSSSYAAGAYYAGSYGASASGALILIIIILLVVVPVPPPVPAVVVAVAPAVVAVPVVAPAAVVVPVAPV
VVAVPVVAARPVAALPAVAVSAISIVVVVPVVAVVVPVSVVAVVVPSPSVSSSLSAAYVSSAAALLAALAALASSPPP

>gfp-reference-P42212
MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQ
HDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNG
IKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
```

</details>

<details>
<summary>RuBisCO: generated and reference sequences</summary>

```fasta
>rubisco-jev
MLILVAGAGLILVAFASVTVPVGGAPGGGLVVIVIIVASGDDREHKNGSAHESAHGSALRHASELHRGGRGAEAHSGHAE
RKAEVHAGGHAREGSHIGDDGGHDSGAFENLSNLSSNLFALHEDGEGHGGGGLEGLELLLLPAEPGAAGAAGGAPVGAGV
GEVIIIIIEGAVAEDDERHKRHHERHDHEHDHSRADEDEHNHHRVDAHHAAEGHGGAGGWGGGGGAAGPGGGGGAGGAGG
GGLAAVEALAALLSLSESSANSAESNLHHSLSLESSSLAHLENLLLSVEEGAGHLPGPLPGPAGPGGAAGGAAAAALVLV
ELVVVLVLVAVLLLVLLPAAGGPGAAAALVVVVPAAGADEGHHHHHAGGGAADEHDHHHDHAHEDDAAEGAAAHDDASHQ
HQSGGGGSAGGAGAGAGGGGAAGAAAAAAAASLNLALLLASALALLSASAENLASALDAHSAGLAASEHSDEAHS

>rubisco-reference-P00875
MSPQTETKASVEFKAGVKDYKLTYYTPEYETLDTDILAAFRVSPQPGVPPEEAGAAVAAESSTGTWTTVWTDGLTNLDRY
KGRCYHIEPVAGEENQYICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPVAYVKTFQGPPHGIQVERD
KLNKYGRPLLGCTIKPKLGLSAKNYGRAVYECLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEALYKAQAETGEIKGHYL
NATAGTCEDMMKRAVFARELGVPIVMHDYLTGGFTANTTLSHYCRDNGLLLHIHRAMHAVIDRQKNHGMHFRVLAKALRL
SGGDHIHSGTVVGKLEGERDITLGFVDLLRDDYTEKDRSRGIYFTQSWVSTPGVLPVASGGIHVWHMPALTEIFGDDSVL
QFGGGTLGHPWGNAPGAVANRVALEACVQARNEGRDLAREGNTIIREATKWSPELAAACEVWKEIKFEFPAMDTV
```

</details>

The [examples folder](examples/) contains FASTA sequences, Jev traces, predicted structures, and experimental reference chains. The [comparison script](scripts/compare_function_examples.py) calculates the metrics, and the [PyMOL script](scripts/render_cartoon_comparison.py) renders the cartoon views.

## Optional: search known proteins

Use local [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) to compare candidates with reviewed UniProtKB/Swiss-Prot sequences. On macOS:

```bash
brew install blast
jev-protein-design-db
jev-protein-design-check --fasta examples/gfp-jev.fasta
```

On other systems, install NCBI BLAST+ before running the database setup. The database stays in the ignored `data/` directory. Add `--search` to a generation command to search automatically after generation. Hits report sequence similarity and reference annotations; they do not confirm function.

## Install as a Codex Skill

After cloning the repository:

```bash
mkdir -p ~/.codex/skills
cp -R skills/jev-protein-design ~/.codex/skills/
```

Ask Codex to use **`$jev-protein-design`**. The [Skill](skills/jev-protein-design/SKILL.md) guides setup, generation, and comparison. Each user supplies their own key locally.

## Development

```bash
python3 -m unittest discover -s tests -v
```

## 中文简介

这是一个用 Jev 逐位选择氨基酸的蛋白设计玩具项目。初始目标包含避免机械重复的约束，之后每次仅更新已生成的序列。示例展示生成序列、Boltz 预测结构与真实蛋白的对照，尚无功能实验验证。

## License

MIT
