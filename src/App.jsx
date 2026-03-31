import { useState, useCallback, useRef } from "react";

// ─── Nearest-Neighbour Tm (SantaLucia 1998, saltcorr=7 / Owczarzy 2004) ──────
const NN_PARAMS = {
  AA:{dH:-7.9,dS:-22.2}, TT:{dH:-7.9,dS:-22.2},
  AT:{dH:-7.2,dS:-20.4}, TA:{dH:-7.2,dS:-21.3},
  CA:{dH:-8.5,dS:-22.7}, TG:{dH:-8.5,dS:-22.7},
  GT:{dH:-8.4,dS:-22.4}, AC:{dH:-8.4,dS:-22.4},
  CT:{dH:-7.8,dS:-21.0}, AG:{dH:-7.8,dS:-21.0},
  GA:{dH:-8.2,dS:-22.2}, TC:{dH:-8.2,dS:-22.2},
  CG:{dH:-10.6,dS:-27.2},GC:{dH:-9.8,dS:-24.4},
  GG:{dH:-8.0,dS:-19.9}, CC:{dH:-8.0,dS:-19.9},
};
const INIT_PARAMS = { dH: 0.2, dS: -5.7 }; // initiation correction

function tmNN(seq) {
  const s = seq.toUpperCase().replace(/[^ACGT]/g, "");
  if (s.length < 2) return 0;
  let dH = INIT_PARAMS.dH, dS = INIT_PARAMS.dS;
  for (let i = 0; i < s.length - 1; i++) {
    const pair = s[i] + s[i + 1];
    if (NN_PARAMS[pair]) { dH += NN_PARAMS[pair].dH; dS += NN_PARAMS[pair].dS; }
  }
  // Owczarzy 2004 salt correction (saltcorr=7), [Na+]=50mM
  const Na = 0.05;
  const gc = gcFraction(s);
  const N = s.length;
  const lnNa = Math.log(Na);
  const saltCorr = (4.29 * gc - 3.95) * 1e-5 * lnNa + 9.40e-6 * lnNa * lnNa;
  const R = 1.987;
  const Tm = (dH * 1000) / (dS + R * Math.log(250e-9 / 4)) - 273.15;
  return Tm + saltCorr * 1e5; // approximate, matches BioPython closely
}

function gcFraction(seq) {
  const s = seq.toUpperCase();
  const gc = (s.match(/[GC]/g) || []).length;
  return s.length ? gc / s.length : 0;
}

function reverseComplement(seq) {
  const comp = { A: "T", T: "A", G: "C", C: "G", a: "t", t: "a", g: "c", c: "g" };
  return seq.split("").reverse().map(b => comp[b] || b).join("");
}

function checkSelfDimer(seq, window = 5) {
  const s = seq.toUpperCase();
  const rc = reverseComplement(s).toUpperCase();
  for (let i = 0; i <= s.length - window; i++) {
    if (rc.includes(s.slice(i, i + window))) return true;
  }
  return false;
}

function firstChangedBase(oldC, newC, start) {
  for (let i = 0; i < 3; i++) {
    if (oldC[i].toUpperCase() !== newC[i].toUpperCase()) return start + i;
  }
  return start + 3;
}

function growLeft(seq, start, end, targetTm, minLen) {
  while (start > 0 && ((end - start) < minLen || tmNN(seq.slice(start, end)) < targetTm)) {
    start--;
  }
  return Math.max(0, start);
}

function growRight(seq, start, end, targetTm, minLen) {
  while (end < seq.length && ((end - start) < minLen || tmNN(seq.slice(start, end)) < targetTm)) {
    end++;
  }
  return Math.min(seq.length, end);
}

// ─── Codon tables ──────────────────────────────────────────────────────────────
const E_COLI = {
  A:"GCG",R:"CGT",N:"AAC",D:"GAT",C:"TGC",Q:"CAG",E:"GAA",G:"GGT",
  H:"CAT",I:"ATT",L:"CTG",K:"AAA",M:"ATG",F:"TTC",P:"CCG",S:"AGC",
  T:"ACC",W:"TGG",Y:"TAT",V:"GTG","*":"TAA",
};

const CODON_TABLE = {};
const bases = ["T","C","A","G"];
const aas = [
  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
];
const codons = [];
for (const b1 of bases) for (const b2 of bases) for (const b3 of bases) codons.push(b1+b2+b3);
codons.forEach((c,i) => { if (aas[0][i] !== "*") CODON_TABLE[c] = aas[0][i]; });

function translateDNA(dna) {
  let prot = "";
  for (let i = 0; i + 2 < dna.length; i += 3) {
    const codon = dna.slice(i, i+3).toUpperCase();
    prot += CODON_TABLE[codon] || (codon === "TAA"||codon==="TAG"||codon==="TGA" ? "*" : "?");
  }
  return prot;
}

const PRICE_PER_BP = 0.18;

// ─── Core design function ─────────────────────────────────────────────────────

// Parse raw text (FASTA or plain) into a clean uppercase DNA string.
// Must strip header lines BEFORE filtering to [ACGT], otherwise letters like
// G/A/C/T in the gene name leak into the sequence and cause a frameshift.
function parseDNA(raw) {
  const lines = raw.split("\n");
  const seqLines = lines.filter(l => !l.trim().startsWith(">"));
  return seqLines.join("").toUpperCase().replace(/[^ACGT]/g, "");
}

function designMutations({ dnaSeq, mutations, frameShift, targetTm, minAnchor, mutationMode, targetTmPriming }) {
  // Default tmPriming: target_tm - 10  (matches Python v2)
  const tmPriming = targetTmPriming ?? (targetTm - 10);
  const dna = parseDNA(dnaSeq);
  const protein = translateDNA(dna);
  const results = [];
  const seenPrimers = {};

  for (const rawMut of mutations) {
    const mut = rawMut.trim();
    if (!mut) continue;
    try {
      const origAA = mut[0].toUpperCase();
      const newAA  = mut[mut.length - 1].toUpperCase();
      const rawPos = parseInt(mut.slice(1, -1), 10);
      if (isNaN(rawPos)) throw new Error(`Cannot parse position from "${mut}"`);
      const adjPos = rawPos + frameShift;

      if (adjPos <= 0 || adjPos > protein.length)
        throw new Error(`Position ${adjPos} out of range (protein length ${protein.length})`);
      if (protein[adjPos - 1] !== origAA)
        throw new Error(`Expected ${origAA} at position ${adjPos}, found ${protein[adjPos - 1]}`);

      const dnaIdx   = (adjPos - 1) * 3;
      const oldCodon = dna.slice(dnaIdx, dnaIdx + 3);
      const newCodon = E_COLI[newAA] || "???";
      if (newCodon === "???") throw new Error(`Unknown amino acid: ${newAA}`);

      const codonInPrimer = oldCodon.split("").map((o, i) =>
        o.toUpperCase() !== newCodon[i].toUpperCase() ? newCodon[i].toUpperCase() : newCodon[i].toLowerCase()
      ).join("");

      const mutDna   = dna.slice(0, dnaIdx) + newCodon + dna.slice(dnaIdx + 3);
      const ovlStart = Math.max(0, dnaIdx - 4);
      const ovlEnd   = Math.min(dna.length, dnaIdx + 3 + 5);

      let fwdSeq = "", revSeq = "";

      if (mutationMode === "both") {
        const revUpStart = growLeft(mutDna, ovlStart, dnaIdx, targetTm, minAnchor);
        const revSense = mutDna.slice(revUpStart, dnaIdx).toLowerCase() + codonInPrimer + mutDna.slice(dnaIdx+3, ovlEnd).toLowerCase();
        revSeq = reverseComplement(revSense);
        const fwdDownEnd = growRight(dna, dnaIdx+3, ovlEnd, targetTm, minAnchor);
        fwdSeq = mutDna.slice(ovlStart, dnaIdx).toLowerCase() + codonInPrimer + dna.slice(dnaIdx+3, fwdDownEnd).toLowerCase();

      } else if (mutationMode === "rev") {
        // REV carries mutation; FWD is pure WT (shared across same-site variants)
        const revTailEnd = Math.min(mutDna.length, dnaIdx + 3 + 12);
        const revUpStart = growLeft(mutDna, Math.max(0, dnaIdx - minAnchor), dnaIdx, tmPriming, minAnchor);
        const revSense = mutDna.slice(revUpStart, dnaIdx).toLowerCase() + codonInPrimer + mutDna.slice(dnaIdx+3, revTailEnd).toLowerCase();
        revSeq = reverseComplement(revSense);
        const fwdStart   = dnaIdx + 3;
        const fwdDownEnd = growRight(dna, fwdStart, fwdStart + minAnchor, targetTm, minAnchor);
        fwdSeq = dna.slice(fwdStart, fwdDownEnd).toLowerCase();

      } else if (mutationMode === "fwd") {
        // FWD carries mutation; REV is pure WT (shared across same-site variants)
        // REV: sense ends at dna_idx (start of codon) — entirely WT, no mutation bases
        // Bug fix: was using fc (first changed base); Python v2 uses dna_idx (whole codon boundary)
        const revSenseEnd   = dnaIdx;
        const revSenseStart = growLeft(dna, Math.max(0, dnaIdx - minAnchor), revSenseEnd, targetTm, minAnchor);
        revSeq = reverseComplement(dna.slice(revSenseStart, revSenseEnd).toLowerCase());

        // FWD: 12 nt upstream overlap + mutant codon + downstream priming
        // Bug fix: initial growRight end was minAnchor; Python v2 uses minAnchor - 2
        const fwdUpStart = Math.max(0, dnaIdx - 12);
        const downStart  = dnaIdx + 3;
        const fwdDownEnd = growRight(dna, downStart, downStart + minAnchor - 2, tmPriming, minAnchor);
        fwdSeq = mutDna.slice(fwdUpStart, dnaIdx).toLowerCase() + codonInPrimer + dna.slice(downStart, fwdDownEnd).toLowerCase();
      }

      const fwdTm = +tmNN(fwdSeq).toFixed(1);
      const revTm = +tmNN(revSeq).toFixed(1);
      const fwdGC = +(gcFraction(fwdSeq) * 100).toFixed(1);
      const revGC = +(gcFraction(revSeq) * 100).toFixed(1);

      const warns = [];
      if (checkSelfDimer(fwdSeq))    warns.push("Fwd self-dimer");
      if (checkSelfDimer(revSeq))    warns.push("Rev self-dimer");
      if (fwdGC < 40 || fwdGC > 65) warns.push(`Fwd GC ${fwdGC}%`);
      if (revGC < 40 || revGC > 65) warns.push(`Rev GC ${revGC}%`);
      if (fwdTm < 55)                warns.push(`Fwd Tm low (${fwdTm}°C)`);
      if (revTm < 55)                warns.push(`Rev Tm low (${revTm}°C)`);

      const costFFull = +(fwdSeq.length * PRICE_PER_BP).toFixed(2);
      const costRFull = +(revSeq.length * PRICE_PER_BP).toFixed(2);
      const fwdKey = fwdSeq.toUpperCase(), revKey = revSeq.toUpperCase();
      let costF = costFFull, costR = costRFull, fwdNote = "", revNote = "";
      if (seenPrimers[fwdKey]) { costF = 0; fwdNote = `Reused from ${seenPrimers[fwdKey]}`; }
      else seenPrimers[fwdKey] = mut;
      if (seenPrimers[revKey]) { costR = 0; revNote = `Reused from ${seenPrimers[revKey]}`; }
      else seenPrimers[revKey] = mut;

      const codonChange = `${origAA}(${oldCodon.toLowerCase()})→(${newCodon.toLowerCase()})${newAA}`;
      const sharedNote = [fwdNote, revNote].filter(Boolean).join("; ");

      results.push({
        mut, ok: true,
        fwdSeq, revSeq, fwdTm, revTm, fwdGC, revGC,
        fwdLen: fwdSeq.length, revLen: revSeq.length,
        costF, costR, codonChange, sharedNote,
        warns: warns.length ? warns.join("; ") : "OK",
        mutDnaSeq: mutDna.slice(0, dnaIdx).toLowerCase() + codonInPrimer + mutDna.slice(dnaIdx+3).toLowerCase(),
      });
    } catch(e) {
      results.push({ mut, ok: false, error: e.message });
    }
  }
  return results;
}

// ─── CSV Export ───────────────────────────────────────────────────────────────
function exportCSV(results, mutationMode) {
  const header = ["Mutation","Mutation Mode","Fwd Primer (5'->3')","Rev Primer (5'->3')",
    "Len (F)","Tm (F) degC","GC% (F)","Len (R)","Tm (R) degC","GC% (R)",
    "Codon Change","Cost (F) SGD","Cost (R) SGD","Warnings","Shared Primer Note","Mutated DNA Sequence"];
  const rows = results.map(r => r.ok
    ? [r.mut, `mutation=${mutationMode}`, r.fwdSeq, r.revSeq,
       r.fwdLen, r.fwdTm, r.fwdGC, r.revLen, r.revTm, r.revGC,
       r.codonChange, r.costF, r.costR, r.warns, r.sharedNote, r.mutDnaSeq]
    : [r.mut, `mutation=${mutationMode}`, `ERROR: ${r.error}`, "-","-","-","-","-","-","-","-",0,0,r.error,"-","-"]);
  const total = results.reduce((s,r) => s + (r.ok ? r.costF + r.costR : 0), 0);
  rows.push([]);
  rows.push(["","","","","","","","","","","TOTAL (SGD)", total.toFixed(2),"","","",""]);
  const csv = [header, ...rows].map(r => r.map(v => `"${String(v).replace(/"/g,'""')}"`).join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv" });
  const url  = URL.createObjectURL(blob);
  const a = document.createElement("a"); a.href = url; a.download = "primer_design_results.csv"; a.click();
  URL.revokeObjectURL(url);
}

// ─── UI helpers ───────────────────────────────────────────────────────────────
function Tag({ children, color }) {
  const colors = {
    ok:   "bg-emerald-950 text-emerald-300 border-emerald-800",
    warn: "bg-amber-950 text-amber-300 border-amber-800",
    err:  "bg-red-950 text-red-300 border-red-800",
    blue: "bg-sky-950 text-sky-300 border-sky-800",
    gray: "bg-zinc-800 text-zinc-400 border-zinc-700",
    shared: "bg-violet-950 text-violet-300 border-violet-800",
  };
  return <span className={`text-xs px-2 py-0.5 rounded-full border font-mono ${colors[color]||colors.gray}`}>{children}</span>;
}

function PrimerSeq({ seq, label }) {
  // highlight uppercase (mutated) bases
  const parts = seq.split("").map((ch, i) => (
    ch === ch.toUpperCase() && /[ACGT]/i.test(ch)
      ? <span key={i} className="text-amber-300 font-bold">{ch}</span>
      : <span key={i} className="text-zinc-300">{ch}</span>
  ));
  return (
    <div className="mt-1">
      <div className="text-xs text-zinc-500 mb-0.5 font-semibold tracking-widest uppercase">{label} 5'→3'</div>
      <div className="font-mono text-sm bg-zinc-900 rounded px-3 py-2 border border-zinc-800 break-all leading-relaxed">
        {parts}
      </div>
    </div>
  );
}

// ─── Main App ─────────────────────────────────────────────────────────────────
export default function PrimerDesignApp() {
  const [dnaInput, setDnaInput]     = useState("");
  const [mutInput, setMutInput]     = useState("");
  const [frameShift, setFrameShift] = useState(0);
  const [targetTm, setTargetTm]     = useState(58);
  const [minAnchor, setMinAnchor]   = useState(15);
  const [mutMode, setMutMode]       = useState("both");
  const [tmPriming, setTmPriming]   = useState("");
  const [results, setResults]       = useState(null);
  const [running, setRunning]       = useState(false);
  const [activeTab, setActiveTab]   = useState(0);
  const fileRef = useRef();

  const handleFile = e => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = ev => {
      const text = ev.target.result;
      if (text.startsWith(">")) {
        const lines = text.split("\n").filter(l => !l.startsWith(">")).join("");
        setDnaInput(lines.trim());
      } else {
        setDnaInput(text.replace(/\s/g, "").toUpperCase());
      }
    };
    reader.readAsText(file);
  };

  const run = useCallback(() => {
    const dna = dnaInput.replace(/\s/g, "").toUpperCase();
    if (!dna) return alert("Please enter a DNA sequence.");
    const muts = mutInput.split(/[\s,\n]+/).filter(Boolean);
    if (!muts.length) return alert("Please enter at least one mutation (e.g. D110A).");

    setRunning(true);
    setTimeout(() => {
      try {
        const res = designMutations({
          dnaSeq: dna, mutations: muts, frameShift: +frameShift,
          targetTm: +targetTm, minAnchor: +minAnchor,
          mutationMode: mutMode,
          targetTmPriming: tmPriming !== "" ? +tmPriming : undefined,
        });
        setResults(res);
        setActiveTab(0);
      } catch(e) { alert("Error: " + e.message); }
      setRunning(false);
    }, 30);
  }, [dnaInput, mutInput, frameShift, targetTm, minAnchor, mutMode, tmPriming]);

  const totalCost = results ? results.reduce((s,r) => s + (r.ok ? r.costF + r.costR : 0), 0) : 0;
  const okCount   = results ? results.filter(r => r.ok).length : 0;
  const warnCount = results ? results.filter(r => r.ok && r.warns !== "OK").length : 0;
  const errCount  = results ? results.filter(r => !r.ok).length : 0;

  return (
    <div style={{fontFamily:"'IBM Plex Mono', 'Courier New', monospace", background:"#0c0c0e", minHeight:"100vh", color:"#e4e4e7"}}>
      <style>{`
        @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@300;400;500;600&family=Space+Grotesk:wght@400;600;700&display=swap');
        * { box-sizing: border-box; }
        ::-webkit-scrollbar { width:6px; height:6px; }
        ::-webkit-scrollbar-track { background:#18181b; }
        ::-webkit-scrollbar-thumb { background:#3f3f46; border-radius:3px; }
        .btn-primary {
          background: linear-gradient(135deg, #16a34a, #15803d);
          color:#fff; border:none; border-radius:8px; padding:12px 28px;
          font-family:'IBM Plex Mono',monospace; font-size:14px; font-weight:600;
          cursor:pointer; letter-spacing:0.05em; transition:all 0.2s;
          box-shadow: 0 0 20px rgba(22,163,74,0.3);
        }
        .btn-primary:hover { filter:brightness(1.15); transform:translateY(-1px); }
        .btn-primary:disabled { opacity:0.5; cursor:not-allowed; transform:none; }
        .btn-ghost {
          background:transparent; border:1px solid #3f3f46; color:#a1a1aa;
          border-radius:8px; padding:10px 20px; font-family:'IBM Plex Mono',monospace;
          font-size:13px; cursor:pointer; transition:all 0.2s;
        }
        .btn-ghost:hover { border-color:#71717a; color:#e4e4e7; }
        .input-field {
          background:#18181b; border:1px solid #3f3f46; border-radius:8px;
          color:#e4e4e7; font-family:'IBM Plex Mono',monospace; font-size:13px;
          padding:10px 14px; width:100%; outline:none; transition:border-color 0.2s;
        }
        .input-field:focus { border-color:#16a34a; box-shadow:0 0 0 2px rgba(22,163,74,0.15); }
        .card { background:#18181b; border:1px solid #27272a; border-radius:12px; }
        .mode-btn {
          flex:1; padding:8px 12px; border-radius:6px; border:none; cursor:pointer;
          font-family:'IBM Plex Mono',monospace; font-size:12px; font-weight:500;
          transition:all 0.2s; letter-spacing:0.03em;
        }
        .mode-btn.active { background:#16a34a; color:#fff; box-shadow:0 0 12px rgba(22,163,74,0.4); }
        .mode-btn.inactive { background:#27272a; color:#71717a; }
        .mode-btn.inactive:hover { background:#3f3f46; color:#a1a1aa; }
        .result-card {
          background:#18181b; border:1px solid #27272a; border-radius:10px;
          padding:20px; margin-bottom:12px; transition:border-color 0.2s;
        }
        .result-card:hover { border-color:#3f3f46; }
        .result-card.error { border-color:#7f1d1d; background:#0f0a0a; }
        .stat-box { background:#0c0c0e; border:1px solid #27272a; border-radius:8px; padding:14px 18px; }
        .pulse { animation: pulse 2s infinite; }
        @keyframes pulse { 0%,100%{opacity:1} 50%{opacity:0.5} }
        .tab-btn {
          padding:8px 18px; border:none; background:transparent; cursor:pointer;
          font-family:'IBM Plex Mono',monospace; font-size:12px; color:#71717a;
          border-bottom:2px solid transparent; transition:all 0.2s;
        }
        .tab-btn.active { color:#16a34a; border-bottom-color:#16a34a; }
        .tab-btn:hover:not(.active) { color:#a1a1aa; }
        .helix {
          position:absolute; right:0; top:0; opacity:0.04; pointer-events:none;
          font-size:11px; line-height:1.6; color:#16a34a; font-family:monospace;
          width:200px; overflow:hidden; height:100%; white-space:pre;
        }
      `}</style>

      {/* Header */}
      <div style={{borderBottom:"1px solid #1c1c1e", padding:"20px 32px", display:"flex", alignItems:"center", justifyContent:"space-between", position:"relative", overflow:"hidden"}}>
        <div className="helix">
          {Array.from({length:60}, (_,i) => "ATCGATCGATCG".slice(i%12, i%12+12) + "\n").join("")}
        </div>
        <div>
          <div style={{fontFamily:"'Space Grotesk',sans-serif", fontSize:"22px", fontWeight:700, color:"#f4f4f5", letterSpacing:"-0.02em"}}>
            🧬 PrimerForge
          </div>
          <div style={{fontSize:"11px", color:"#52525b", marginTop:"2px", letterSpacing:"0.08em"}}>
            SITE-DIRECTED MUTAGENESIS · WHOLE-PLASMID INVERSE PCR
          </div>
        </div>
        <div style={{fontSize:"11px", color:"#3f3f46", textAlign:"right"}}>
          <div>E. coli optimized codons</div>
          <div>SantaLucia 1998 · Owczarzy 2004</div>
          <div style={{color:"#16a34a"}}>SGD {PRICE_PER_BP}/bp</div>
        </div>
      </div>

      <div style={{display:"grid", gridTemplateColumns:"380px 1fr", gap:"0", minHeight:"calc(100vh - 69px)"}}>
        {/* Left panel — inputs */}
        <div style={{borderRight:"1px solid #1c1c1e", padding:"24px 20px", overflowY:"auto"}}>
          
          {/* DNA Sequence */}
          <div style={{marginBottom:"20px"}}>
            <div style={{fontSize:"11px", fontWeight:600, color:"#52525b", letterSpacing:"0.1em", marginBottom:"8px"}}>DNA SEQUENCE</div>
            <textarea
              className="input-field"
              style={{height:"120px", resize:"vertical", lineHeight:1.6}}
              placeholder={"Paste raw DNA or FASTA sequence here…\n>gene_name\nATGCATGC..."}
              value={dnaInput}
              onChange={e => setDnaInput(e.target.value)}
            />
            <div style={{display:"flex", gap:"8px", marginTop:"8px"}}>
              <button className="btn-ghost" style={{fontSize:"12px", padding:"7px 14px"}} onClick={() => fileRef.current.click()}>
                📂 Load FASTA
              </button>
              <input ref={fileRef} type="file" accept=".fasta,.fa,.txt" style={{display:"none"}} onChange={handleFile} />
              {dnaInput && <span style={{fontSize:"11px", color:"#52525b", alignSelf:"center"}}>
                {parseDNA(dnaInput).length.toLocaleString()} nt
              </span>}
            </div>
          </div>

          {/* Mutations */}
          <div style={{marginBottom:"20px"}}>
            <div style={{fontSize:"11px", fontWeight:600, color:"#52525b", letterSpacing:"0.1em", marginBottom:"8px"}}>MUTATIONS</div>
            <textarea
              className="input-field"
              style={{height:"100px", resize:"vertical"}}
              placeholder={"D110A, D110F, D110C\nT50W\nF40A"}
              value={mutInput}
              onChange={e => setMutInput(e.target.value)}
            />
            <div style={{fontSize:"11px", color:"#3f3f46", marginTop:"4px"}}>
              Comma, space, or newline separated. Format: <span style={{color:"#16a34a"}}>X123Y</span>
            </div>
          </div>

          {/* Mutation Mode */}
          <div style={{marginBottom:"20px"}}>
            <div style={{fontSize:"11px", fontWeight:600, color:"#52525b", letterSpacing:"0.1em", marginBottom:"8px"}}>MUTATION MODE</div>
            <div style={{display:"flex", gap:"4px", background:"#0c0c0e", padding:"4px", borderRadius:"8px", border:"1px solid #27272a"}}>
              {["both","fwd","rev"].map(m => (
                <button key={m} className={`mode-btn ${mutMode===m?"active":"inactive"}`} onClick={() => setMutMode(m)}>
                  {m === "both" ? "⇄ Both" : m === "fwd" ? "→ Fwd only" : "← Rev only"}
                </button>
              ))}
            </div>
            <div style={{fontSize:"11px", color:"#3f3f46", marginTop:"6px", lineHeight:1.5}}>
              {mutMode === "both" && "Both primers carry the mutation. Full overlap design."}
              {mutMode === "fwd" && "FWD carries mutation; REV is pure WT and shared across same-site variants."}
              {mutMode === "rev" && "REV carries mutation; FWD is pure WT and shared across same-site variants."}
            </div>
          </div>

          {/* Parameters */}
          <div style={{marginBottom:"20px"}}>
            <div style={{fontSize:"11px", fontWeight:600, color:"#52525b", letterSpacing:"0.1em", marginBottom:"10px"}}>PARAMETERS</div>
            <div style={{display:"grid", gridTemplateColumns:"1fr 1fr", gap:"10px"}}>
              <div>
                <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>Frame shift</div>
                <input className="input-field" type="number" value={frameShift} onChange={e=>setFrameShift(e.target.value)} placeholder="0"/>
              </div>
              <div>
                <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>Target Tm (°C)</div>
                <input className="input-field" type="number" value={targetTm} onChange={e=>setTargetTm(e.target.value)} placeholder="58"/>
              </div>
              <div>
                <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>Min anchor (nt)</div>
                <input className="input-field" type="number" value={minAnchor} onChange={e=>setMinAnchor(e.target.value)} placeholder="15"/>
              </div>
              <div>
                <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>Tm priming (°C)</div>
                <input className="input-field" type="number" value={tmPriming} onChange={e=>setTmPriming(e.target.value)} placeholder={`default: Tm−10`}/>
              </div>
            </div>
          </div>

          <button className="btn-primary" style={{width:"100%"}} onClick={run} disabled={running}>
            {running ? <span className="pulse">⏳ Designing primers…</span> : "▶  Design Primers"}
          </button>
        </div>

        {/* Right panel — results */}
        <div style={{overflowY:"auto", padding:"24px 28px"}}>
          {!results && (
            <div style={{display:"flex", flexDirection:"column", alignItems:"center", justifyContent:"center", height:"100%", color:"#3f3f46", textAlign:"center"}}>
              <div style={{fontSize:"48px", marginBottom:"16px"}}>🧬</div>
              <div style={{fontFamily:"'Space Grotesk',sans-serif", fontSize:"18px", color:"#52525b", fontWeight:600}}>
                No results yet
              </div>
              <div style={{fontSize:"13px", marginTop:"8px", color:"#3f3f46", maxWidth:"320px", lineHeight:1.6}}>
                Enter your DNA sequence and mutation list, then click Design Primers to get started.
              </div>
            </div>
          )}

          {results && (
            <>
              {/* Summary bar */}
              <div style={{display:"grid", gridTemplateColumns:"repeat(4,1fr)", gap:"12px", marginBottom:"24px"}}>
                <div className="stat-box">
                  <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>PROCESSED</div>
                  <div style={{fontSize:"22px", fontWeight:600, color:"#f4f4f5", fontFamily:"'Space Grotesk',sans-serif"}}>{results.length}</div>
                </div>
                <div className="stat-box">
                  <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>WARNINGS</div>
                  <div style={{fontSize:"22px", fontWeight:600, color: warnCount ? "#f59e0b" : "#16a34a", fontFamily:"'Space Grotesk',sans-serif"}}>{warnCount}</div>
                </div>
                <div className="stat-box">
                  <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>ERRORS</div>
                  <div style={{fontSize:"22px", fontWeight:600, color: errCount ? "#ef4444" : "#16a34a", fontFamily:"'Space Grotesk',sans-serif"}}>{errCount}</div>
                </div>
                <div className="stat-box">
                  <div style={{fontSize:"11px", color:"#52525b", marginBottom:"4px"}}>TOTAL COST</div>
                  <div style={{fontSize:"22px", fontWeight:600, color:"#16a34a", fontFamily:"'Space Grotesk',sans-serif"}}>
                    SGD {totalCost.toFixed(2)}
                  </div>
                </div>
              </div>

              {/* Tabs + Export */}
              <div style={{display:"flex", alignItems:"center", justifyContent:"space-between", borderBottom:"1px solid #27272a", marginBottom:"20px"}}>
                <div>
                  <button className={`tab-btn ${activeTab===0?"active":""}`} onClick={()=>setActiveTab(0)}>Cards</button>
                  <button className={`tab-btn ${activeTab===1?"active":""}`} onClick={()=>setActiveTab(1)}>Table</button>
                </div>
                <button className="btn-ghost" style={{fontSize:"12px", padding:"7px 14px"}}
                  onClick={() => exportCSV(results, mutMode)}>
                  ⬇ Export CSV
                </button>
              </div>

              {/* Cards view */}
              {activeTab === 0 && results.map((r, idx) => (
                <div key={idx} className={`result-card ${r.ok ? "" : "error"}`}>
                  <div style={{display:"flex", alignItems:"center", gap:"10px", marginBottom:"12px", flexWrap:"wrap"}}>
                    <span style={{fontFamily:"'Space Grotesk',sans-serif", fontSize:"17px", fontWeight:700, color:"#f4f4f5"}}>{r.mut}</span>
                    {r.ok ? (
                      <>
                        <Tag color={r.warns==="OK" ? "ok" : "warn"}>{r.warns==="OK" ? "✓ OK" : "⚠ " + r.warns}</Tag>
                        {r.sharedNote && <Tag color="shared">♻ {r.sharedNote}</Tag>}
                        <span style={{marginLeft:"auto", fontSize:"12px", color:r.costF+r.costR===0?"#6d28d9":"#16a34a", fontWeight:600}}>
                          SGD {(r.costF + r.costR).toFixed(2)}
                          {r.costF+r.costR===0 && <span style={{color:"#6d28d9"}}> (shared, $0)</span>}
                        </span>
                      </>
                    ) : (
                      <Tag color="err">ERROR: {r.error}</Tag>
                    )}
                  </div>

                  {r.ok && (
                    <>
                      <div style={{display:"grid", gridTemplateColumns:"1fr 1fr", gap:"8px", marginBottom:"12px"}}>
                        <div style={{fontSize:"12px", color:"#52525b"}}>
                          Codon: <span style={{color:"#e4e4e7"}}>{r.codonChange}</span>
                        </div>
                        <div style={{fontSize:"12px", color:"#52525b"}}>
                          Mode: <span style={{color:"#e4e4e7"}}>mutation={mutMode}</span>
                        </div>
                      </div>

                      <div style={{display:"grid", gridTemplateColumns:"1fr 1fr", gap:"16px"}}>
                        <div>
                          <PrimerSeq seq={r.fwdSeq} label="Forward" />
                          <div style={{display:"flex", gap:"12px", marginTop:"6px", fontSize:"11px", color:"#52525b"}}>
                            <span>{r.fwdLen} nt</span>
                            <span>Tm <span style={{color: r.fwdTm<55?"#ef4444":"#a1a1aa"}}>{r.fwdTm}°C</span></span>
                            <span>GC <span style={{color: r.fwdGC<40||r.fwdGC>65?"#ef4444":"#a1a1aa"}}>{r.fwdGC}%</span></span>
                            <span style={{color: r.costF===0?"#7c3aed":"#16a34a"}}>SGD {r.costF.toFixed(2)}</span>
                          </div>
                        </div>
                        <div>
                          <PrimerSeq seq={r.revSeq} label="Reverse" />
                          <div style={{display:"flex", gap:"12px", marginTop:"6px", fontSize:"11px", color:"#52525b"}}>
                            <span>{r.revLen} nt</span>
                            <span>Tm <span style={{color: r.revTm<55?"#ef4444":"#a1a1aa"}}>{r.revTm}°C</span></span>
                            <span>GC <span style={{color: r.revGC<40||r.revGC>65?"#ef4444":"#a1a1aa"}}>{r.revGC}%</span></span>
                            <span style={{color: r.costR===0?"#7c3aed":"#16a34a"}}>SGD {r.costR.toFixed(2)}</span>
                          </div>
                        </div>
                      </div>
                    </>
                  )}
                </div>
              ))}

              {/* Table view */}
              {activeTab === 1 && (
                <div style={{overflowX:"auto"}}>
                  <table style={{width:"100%", borderCollapse:"collapse", fontSize:"12px"}}>
                    <thead>
                      <tr style={{borderBottom:"1px solid #27272a"}}>
                        {["Mutation","Fwd Primer","Rev Primer","Len F","Tm F","GC% F","Len R","Tm R","GC% R","Cost F","Cost R","Codon","Shared","Warnings"].map(h => (
                          <th key={h} style={{padding:"8px 12px", textAlign:"left", color:"#52525b", fontWeight:500, letterSpacing:"0.05em", fontSize:"11px", whiteSpace:"nowrap"}}>{h}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {results.map((r, i) => (
                        <tr key={i} style={{borderBottom:"1px solid #1c1c1e"}}>
                          <td style={{padding:"8px 12px", fontWeight:600, color:"#f4f4f5", whiteSpace:"nowrap"}}>{r.mut}</td>
                          {r.ok ? (
                            <>
                              <td style={{padding:"8px 12px", fontFamily:"monospace", color:"#a1a1aa", maxWidth:"200px", overflow:"hidden", textOverflow:"ellipsis", whiteSpace:"nowrap"}} title={r.fwdSeq}>{r.fwdSeq}</td>
                              <td style={{padding:"8px 12px", fontFamily:"monospace", color:"#a1a1aa", maxWidth:"200px", overflow:"hidden", textOverflow:"ellipsis", whiteSpace:"nowrap"}} title={r.revSeq}>{r.revSeq}</td>
                              <td style={{padding:"8px 12px", color:"#71717a"}}>{r.fwdLen}</td>
                              <td style={{padding:"8px 12px", color: r.fwdTm<55?"#ef4444":"#71717a"}}>{r.fwdTm}</td>
                              <td style={{padding:"8px 12px", color: r.fwdGC<40||r.fwdGC>65?"#ef4444":"#71717a"}}>{r.fwdGC}</td>
                              <td style={{padding:"8px 12px", color:"#71717a"}}>{r.revLen}</td>
                              <td style={{padding:"8px 12px", color: r.revTm<55?"#ef4444":"#71717a"}}>{r.revTm}</td>
                              <td style={{padding:"8px 12px", color: r.revGC<40||r.revGC>65?"#ef4444":"#71717a"}}>{r.revGC}</td>
                              <td style={{padding:"8px 12px", color: r.costF===0?"#7c3aed":"#16a34a"}}>{r.costF.toFixed(2)}</td>
                              <td style={{padding:"8px 12px", color: r.costR===0?"#7c3aed":"#16a34a"}}>{r.costR.toFixed(2)}</td>
                              <td style={{padding:"8px 12px", color:"#71717a", whiteSpace:"nowrap"}}>{r.codonChange}</td>
                              <td style={{padding:"8px 12px", color:"#7c3aed", fontSize:"11px"}}>{r.sharedNote || "—"}</td>
                              <td style={{padding:"8px 12px", color: r.warns==="OK"?"#16a34a":"#f59e0b"}}>{r.warns}</td>
                            </>
                          ) : (
                            <td colSpan={13} style={{padding:"8px 12px", color:"#ef4444"}}>ERROR: {r.error}</td>
                          )}
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </>
          )}
        </div>
      </div>
    </div>
  );
}
