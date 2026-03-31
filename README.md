# 🧬 PrimerForge

**Site-Directed Mutagenesis Primer Design Tool**

Designs primers for whole-plasmid inverse PCR mutagenesis with E. coli optimised codons,
nearest-neighbour Tm calculation (SantaLucia 1998 / Owczarzy 2004), and automatic shared
primer deduplication with cost tracking.

---

## 🚀 Deploy for Free in 5 Minutes

### Option A — Vercel (Recommended)

1. Push this folder to a GitHub repository
2. Go to [vercel.com](https://vercel.com) → **Add New Project**
3. Import your GitHub repository
4. Leave all settings as default — Vercel auto-detects Vite
5. Click **Deploy** → your site is live at `https://primerforge.vercel.app` (or your chosen name)

### Option B — Netlify

1. Push this folder to a GitHub repository
2. Go to [netlify.com](https://netlify.com) → **Add new site** → **Import from Git**
3. Set **Build command**: `npm run build`
4. Set **Publish directory**: `dist`
5. Click **Deploy site**

### Option C — Netlify Drag & Drop (no GitHub needed)

1. Run locally first:
   ```bash
   npm install
   npm run build
   ```
2. Go to [netlify.com/drop](https://app.netlify.com/drop)
3. Drag the `dist/` folder onto the page → instantly live

---

## 💻 Run Locally

```bash
npm install
npm run dev
```

Then open http://localhost:5173

---

## Features

- FASTA file upload or paste sequence directly
- Three mutation modes: **both**, **fwd only**, **rev only**
- Automatic shared primer detection across same-site variants (saves ordering cost)
- Warnings for self-dimers, GC% out of range, low Tm
- Cards view + Table view
- Export results as CSV (matches Python script output format)
- SGD cost calculation at SGD 0.18/bp

## Parameters

| Parameter | Default | Description |
|---|---|---|
| Frame shift | 0 | Nucleotide offset if file doesn't start at codon 1 (e.g. -9) |
| Target Tm | 58°C | Tm target for mutation-carrying primer priming region |
| Min anchor | 15 nt | Minimum priming region length |
| Tm priming | Tm − 10°C | Tm target for pure-WT shared primer |
