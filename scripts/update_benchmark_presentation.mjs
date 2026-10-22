#!/usr/bin/env node
/** Build the complete verified public-acceptance section and update the project deck. */

import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const SKILL = "C:\\Users\\HARDIK\\.codex\\plugins\\cache\\openai-primary-runtime\\presentations\\26.1020.11127\\skills\\presentations";
const PYTHON = "C:\\Users\\HARDIK\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe";
const RUNTIME_NODE_MODULES = "C:\\Users\\HARDIK\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\node\\node_modules";
const DEFAULT_BASE = path.join(ROOT, "deliverables", "LinkerForge_Detailed_Architecture_Data_Provenance_4_Presenters_Final.pptx");
const DEFAULT_FINAL = path.join(ROOT, "deliverables", "LinkerForge_Detailed_Architecture_Data_Provenance_4_Presenters_Linker_Only.pptx");
const BUILD = path.join(ROOT, ".codex", "presentation-build", "multi-protein-benchmark-linker-only");
const CANDIDATE = path.join(BUILD, "candidate-linker-only.pptx");
const CANDIDATE_REPORT = path.join(BUILD, "verified-public-acceptance-run.md");
const PREVIEW_DIR = path.join(BUILD, "previews");
const HIGHEST_CONFIDENCE_IMAGE = path.join(ROOT, ".codex", "presentation-build", "high-confidence-revision", "highest_linker_confidence_3dmol.png");
const ARTIFACT_TOOL = "C:\\Users\\HARDIK\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\node\\node_modules\\@oai\\artifact-tool\\dist\\artifact_tool.mjs";

process.env.RUNTIME_NODE_MODULES ??= RUNTIME_NODE_MODULES;
process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES ??= RUNTIME_NODE_MODULES;

const { FileBlob, PresentationFile } = await import(pathToFileURL(ARTIFACT_TOOL).href);
const utils = await import(pathToFileURL(path.join(SKILL, "container_tools", "artifact_tool_utils.mjs")).href);
const { applyPresentationChartFont, finalizePresentation } = utils;

const C = {
  bg: "091A15", paper: "102720", paperAlt: "143128", ink: "F4F7F5", muted: "AFC0B8", faint: "29433A",
  cyan: "3FC6D6", cyanDark: "167D89", magenta: "E84AAB", green: "2DBD73",
  yellow: "F1C84B", red: "E46760", blue: "4E9BEC", violet: "A88CEB",
};
const FONT = "Arial";
const MONO = "Cascadia Mono";
const W = 1280;

function parseArgs() {
  const parsed = {};
  for (let i = 2; i < process.argv.length; i += 2) parsed[process.argv[i].replace(/^--/, "")] = process.argv[i + 1];
  if (!parsed.summary || !parsed.progress) {
    throw new Error("Required: --summary PATH --progress PATH [--base PATH] [--output PATH] [--report PATH]");
  }
  const summary = path.resolve(parsed.summary);
  return {
    summary,
    progress: path.resolve(parsed.progress),
    base: path.resolve(parsed.base ?? DEFAULT_BASE),
    output: path.resolve(parsed.output ?? DEFAULT_FINAL),
    report: path.resolve(parsed.report ?? path.join(path.dirname(summary), "verified_public_acceptance_run.md")),
  };
}

function shape(slide, geometry, x, y, w, h, fill = "none", lineFill = "none", lineWidth = 0, radius) {
  return slide.shapes.add({
    geometry,
    position: { left: x, top: y, width: w, height: h },
    fill,
    line: { style: "solid", fill: lineFill, width: lineWidth },
    ...(radius ? { borderRadius: radius } : {}),
  });
}

function text(slide, value, x, y, w, h, opts = {}) {
  const box = shape(slide, "textbox", x, y, w, h, opts.fill ?? "none", opts.lineFill ?? "none", opts.lineWidth ?? 0);
  box.text = String(value);
  box.text.style = {
    typeface: opts.typeface ?? FONT,
    fontSize: opts.fontSize ?? 20,
    color: opts.color ?? C.ink,
    bold: opts.bold ?? false,
    italic: opts.italic ?? false,
    alignment: opts.alignment ?? "left",
    verticalAlignment: opts.verticalAlignment ?? "top",
    autoFit: opts.autoFit ?? "shrinkText",
    ...(opts.letterSpacing !== undefined ? { letterSpacing: opts.letterSpacing } : {}),
  };
  if (opts.margin !== undefined) box.text.margin = opts.margin;
  return box;
}

function base(slide, number, section, presenter, source) {
  slide.background.fill = C.bg;
  shape(slide, "rect", 0, 0, W, 8, C.cyan, "none", 0);
  text(slide, "LINKERFORGE", 54, 24, 180, 18, { fontSize: 12, bold: true, color: C.cyan, letterSpacing: 1.5 });
  text(slide, section.toUpperCase(), 234, 24, 650, 18, { fontSize: 12, bold: true, color: C.muted, letterSpacing: 1.2 });
  text(slide, source, 54, 684, 1000, 16, { fontSize: 9, color: C.muted });
  text(slide, `${presenter}  ·  ${String(number).padStart(2, "0")}`, 1080, 681, 150, 18, { fontSize: 10, bold: true, color: C.muted, alignment: "right" });
}

function title(slide, heading, subtitle) {
  text(slide, heading, 54, 60, 1172, 58, { fontSize: 34, bold: true });
  text(slide, subtitle, 56, 120, 1160, 42, { fontSize: 17, color: C.muted });
}

function note(slide, presenter, talk, sources) {
  slide.speakerNotes.text = [
    `Presenter: ${presenter}`, "", talk, "", "Sources and evidence:", ...sources.map((source) => `- ${source}`),
  ].join("\n");
  slide.speakerNotes.setVisible(true);
}

function styleTable(table, rows, columns, fontSize = 11) {
  table.cells.block({ row: 0, column: 0, rowCount: 1, columnCount: columns }).assign({
    fill: C.cyanDark,
    textStyle: { typeface: FONT, fontSize: 12, bold: true, color: "FFFFFF" },
    verticalAlignment: "middle",
  });
  for (let row = 1; row < rows; row += 1) {
    table.cells.block({ row, column: 0, rowCount: 1, columnCount: columns }).assign({
      fill: row % 2 === 0 ? C.paperAlt : C.paper,
      textStyle: { typeface: FONT, fontSize, color: C.ink },
      verticalAlignment: "middle",
    });
  }
  table.borders.assign({ style: "solid", fill: C.faint, width: 1 });
}

function displayPair(pair) {
  return ({ "1UBQ_2GB1": "Ubiquitin–GB1", "1CRN_4PTI": "Crambin–BPTI", "1PGA_1LMB": "Protein G–repressor" })[pair] ?? pair;
}

function value(number, digits = 1) {
  return number === null || number === undefined ? "—" : Number(number).toFixed(digits);
}

function chartValue(number, digits = 2) {
  return Number(Number(number).toFixed(digits));
}

function displaySequence(sequence) {
  return sequence ? sequence : "— direct fusion —";
}

function markdownCell(input) {
  return String(input ?? "—").replaceAll("|", "\\|").replaceAll("\n", " ");
}

async function sha256File(file) {
  return crypto.createHash("sha256").update(await fs.readFile(file)).digest("hex");
}

async function imageBytes(file) {
  const bytes = await fs.readFile(file);
  return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
}

async function pathExists(file) {
  try {
    await fs.access(file);
    return true;
  } catch {
    return false;
  }
}

async function atomicJson(file, object) {
  const temporary = `${file}.tmp`;
  await fs.writeFile(temporary, `${JSON.stringify(object, null, 2)}\n`, "utf8");
  await fs.rename(temporary, file);
}

async function replaceExact(presentation, replacements) {
  const inspected = await presentation.inspect({
    select: "$.slides[*]",
    include: ["id", "index", "elements.id", "elements.text"],
  });
  const found = new Map();
  for (const record of inspected.records) {
    for (const element of record.elements ?? []) {
      if (typeof element.text !== "string") continue;
      for (const [before, after] of Object.entries(replacements)) {
        if (element.text !== before) continue;
        presentation.resolve(element.id).text = after;
        found.set(before, (found.get(before) ?? 0) + 1);
      }
    }
  }
  const missing = Object.keys(replacements).filter((key) => !found.has(key));
  if (missing.length) throw new Error(`Expected source-deck text was not found: ${missing.join(" | ")}`);
}

async function replaceSpeakerSummaryWithLinkerDefinition(presentation, slide) {
  const inspected = await presentation.inspect({
    select: "$.slides[1]",
    include: ["id", "index", "elements.id"],
  });
  const record = inspected.records[0];
  if (!record || record.index !== 1) throw new Error("Could not locate source slide 2");
  presentation.delete((record.elements ?? []).map((element) => element.id));

  base(
    slide,
    2,
    "Biological definition",
    "Presenter 1",
    "Chichili, Kumar & Sivaraman · Protein Science 22 (2013) 153–167 · doi:10.1002/pro.2206",
  );
  slide.background.fill = "000000";
  title(
    slide,
    "Protein linkers connect domains within one polypeptide",
    "Length and amino-acid composition influence domain separation, relative motion, and unwanted interactions",
  );

  text(slide, "DOMAIN A", 86, 200, 230, 24, { fontSize: 13, bold: true, color: C.blue, alignment: "center", letterSpacing: 1.0 });
  shape(slide, "roundRect", 86, 236, 250, 112, "173A4A", C.blue, 2, "rounded-lg");
  text(slide, "Protein domain A", 86, 272, 250, 32, { fontSize: 23, bold: true, alignment: "center", verticalAlignment: "middle" });

  text(slide, "SHORT PEPTIDE SEQUENCE", 428, 200, 424, 24, { fontSize: 13, bold: true, color: C.magenta, alignment: "center", letterSpacing: 0.9 });
  shape(slide, "line", 336, 292, 110, 0, "none", C.muted, 3);
  for (let index = 0; index < 7; index += 1) {
    shape(slide, "ellipse", 454 + index * 49, 266, 38, 38, index % 2 ? "7A245F" : C.magenta, "none", 0);
  }
  text(slide, "example: GGGGS", 430, 319, 420, 24, { fontSize: 15, color: C.muted, alignment: "center" });
  shape(slide, "line", 834, 292, 110, 0, "none", C.muted, 3);

  text(slide, "DOMAIN B", 944, 200, 250, 24, { fontSize: 13, bold: true, color: C.red, alignment: "center", letterSpacing: 1.0 });
  shape(slide, "roundRect", 944, 236, 250, 112, "4A251E", C.red, 2, "rounded-lg");
  text(slide, "Protein domain B", 944, 272, 250, 32, { fontSize: 23, bold: true, alignment: "center", verticalAlignment: "middle" });

  shape(slide, "line", 70, 394, 1140, 0, "none", C.faint, 2);
  text(slide, "SEPARATION", 82, 430, 300, 23, { fontSize: 13, bold: true, color: C.blue, letterSpacing: 1.0 });
  text(slide, "Controls the distance between the connected domains", 82, 466, 300, 64, { fontSize: 20, bold: true });
  shape(slide, "line", 412, 425, 0, 120, "none", C.faint, 2);
  text(slide, "MOBILITY", 458, 430, 300, 23, { fontSize: 13, bold: true, color: C.magenta, letterSpacing: 1.0 });
  text(slide, "Permits or restricts changes in relative domain orientation", 458, 466, 300, 64, { fontSize: 20, bold: true });
  shape(slide, "line", 790, 425, 0, 120, "none", C.faint, 2);
  text(slide, "DOMAIN INTERFERENCE", 836, 430, 340, 23, { fontSize: 13, bold: true, color: C.cyan, letterSpacing: 0.8 });
  text(slide, "Can reduce steric contact and help each domain retain its role", 836, 466, 340, 64, { fontSize: 20, bold: true });

  shape(slide, "line", 70, 570, 1140, 0, "none", C.cyan, 2);
  text(slide, "KEY POINT", 82, 594, 150, 20, { fontSize: 12, bold: true, color: C.cyan, letterSpacing: 1.0 });
  text(slide, "A linker is part of the complete fusion protein. Its useful length and composition therefore depend on the connected domains and the intended biological task.", 230, 585, 950, 54, { fontSize: 18, bold: true });

  note(
    slide,
    "Presenter 1 (major)",
    "Define the biological object before introducing the computational workflow. A protein linker is a peptide segment within one polypeptide chain that connects two domains. The diagram shows the linker as amino-acid residues between Domain A and Domain B. Linkers can change domain separation, relative mobility, and steric interference. The important conclusion is that a linker cannot be interpreted independently of the fusion partners and intended function.",
    [
      "Chichili VPR, Kumar V, Sivaraman J. Linkers in the structural biology of protein–protein interactions. Protein Sci. 2013;22(2):153–167. doi:10.1002/pro.2206.",
      "https://pmc.ncbi.nlm.nih.gov/articles/PMC3588912/",
    ],
  );
}

async function replaceResearchProblemSlide(presentation, slide) {
  const inspected = await presentation.inspect({
    select: "$.slides[2]",
    include: ["id", "index", "elements.id"],
  });
  const record = inspected.records[0];
  if (!record || record.index !== 2) throw new Error("Could not locate source slide 3");
  presentation.delete((record.elements ?? []).map((element) => element.id));

  base(
    slide,
    3,
    "Scientific problem",
    "Presenter 1",
    "Chen, Zaro & Shen · Advanced Drug Delivery Reviews 65 (2013) 1357–1369 · McBride et al. · Proteins 94 (2026) 362–380",
  );
  slide.background.fill = "000000";
  title(
    slide,
    "Why one linker does not fit every fusion",
    "Performance emerges from the partner proteins, attachment sites, linker sequence, and intended function",
  );

  const inputs = [
    ["01", "Partner proteins", "surface shape and terminal geometry", C.blue],
    ["02", "Attachment sites", "junction location and protected residues", C.green],
    ["03", "Linker design", "sequence, length, flexibility, and charge", C.magenta],
    ["04", "Design objective", "separation, orientation, cleavage, or motion", C.yellow],
  ];
  for (const [index, [number, heading, detail, color]] of inputs.entries()) {
    const x = 62 + index * 296;
    text(slide, number, x, 186, 42, 22, { fontSize: 13, bold: true, color });
    text(slide, heading, x + 44, 182, 225, 28, { fontSize: 18, bold: true });
    text(slide, detail, x + 44, 218, 225, 45, { fontSize: 14, color: C.muted });
    if (index < inputs.length - 1) shape(slide, "line", x + 280, 184, 0, 82, "none", C.faint, 1);
  }

  shape(slide, "line", 82, 302, 1116, 0, "none", C.cyan, 2);
  text(slide, "INTERDOMAIN CONFORMATIONAL ENSEMBLE", 82, 330, 1116, 26, { fontSize: 14, bold: true, color: C.cyan, alignment: "center", letterSpacing: 1.0 });
  text(slide, "The complete fusion can sample multiple relative domain positions. A single confident structure may not represent that distribution in solution.", 170, 370, 940, 68, { fontSize: 23, bold: true, alignment: "center" });

  shape(slide, "line", 640, 464, 0, 116, "none", C.faint, 2);
  text(slide, "COMPUTATIONAL EVIDENCE", 82, 472, 500, 22, { fontSize: 13, bold: true, color: C.blue, letterSpacing: 0.9 });
  text(slide, "Structure confidence, PAE, steric clashes, and preservation of each domain", 82, 510, 500, 62, { fontSize: 19, bold: true });
  text(slide, "EXPERIMENTAL EVIDENCE", 698, 472, 500, 22, { fontSize: 13, bold: true, color: C.magenta, letterSpacing: 0.9 });
  text(slide, "Expression, stability, activity, and solution-state measurements such as SAXS or NMR", 698, 510, 500, 62, { fontSize: 19, bold: true });

  shape(slide, "line", 82, 606, 1116, 0, "none", C.magenta, 2);
  text(slide, "LINKERFORGE RESEARCH QUESTION", 82, 624, 270, 20, { fontSize: 12, bold: true, color: C.magenta, letterSpacing: 0.8 });
  text(slide, "Can a pair-conditioned model rank candidates that retain structural plausibility on previously unseen protein families?", 365, 615, 830, 42, { fontSize: 18, bold: true });

  note(
    slide,
    "Presenter 1 (major)",
    "This slide moves from the biological definition to the research problem. Four inputs jointly determine linker behaviour: the two partner proteins, their attachment sites, the candidate linker, and the design objective. The complete fusion may sample many interdomain conformations. Structure prediction can test local plausibility and uncertainty, but experimental measurements remain necessary for expression, stability, activity, and solution ensembles. LinkerForge therefore asks a narrower computational question: can a model rank candidates for a specified pair and generalize to unseen protein families?",
    [
      "Chen X, Zaro JL, Shen WC. Fusion protein linkers: property, design and functionality. Adv Drug Deliv Rev. 2013;65(10):1357–1369. doi:10.1016/j.addr.2012.09.039.",
      "McBride AC et al. Predicting Pose Distribution of Protein Domains Connected by Flexible Linkers Is an Unsolved Problem. Proteins. 2026;94(1):362–380. doi:10.1002/prot.70062.",
    ],
  );
}

async function replacePublishedMethodsSlide(presentation, slide) {
  const inspected = await presentation.inspect({
    select: "$.slides[3]",
    include: ["id", "index", "elements.id"],
  });
  const record = inspected.records[0];
  if (!record || record.index !== 3) throw new Error("Could not locate source slide 4");
  presentation.delete((record.elements ?? []).map((element) => element.id));

  base(
    slide,
    4,
    "Published method landscape",
    "Presenter 1",
    "Primary sources: LINKER (2000), SynLinker (2015), DFLpred (2016), TransDFL (2023), and Li et al. (2026)",
  );
  slide.background.fill = "000000";
  title(
    slide,
    "Published methods address different linker tasks",
    "Generators, residue classifiers, system-specific predictors, and structure models answer different scientific questions",
  );

  const methodRows = [
    ["Method", "Method type", "What it predicts or provides", "Scientific boundary"],
    ["LINKER\n2000", "Library-based generator", "Retrieves PDB loop-derived sequences for a requested linker length", "Does not classify residues or score a specified domain pair"],
    ["SynLinker\n2015", "Database + retrieval/design system", "Searches 2,260 linkers and supports fusion construction and modelling", "Does not learn a general DFL classifier"],
    ["DFLpred\n2016", "Per-residue classifier", "DFL propensity and a DFL versus non-DFL label for each residue", "Detects linker residues; it does not generate or rank fusion linkers"],
    ["TransDFL\n2023", "Transfer-learning per-residue classifier", "DFL versus non-DFL identification across a protein sequence", "Detects natural DFLs; it does not score candidate fusions"],
    ["Li et al.\n2026", "Fusion-specific supervised predictor", "Quantitative activity and high/low activity for one sfGFP–nanobody system", "Not a general DFL detector or universal linker generator"],
  ];
  const methods = slide.tables.add({
    rows: methodRows.length,
    columns: 4,
    left: 54,
    top: 166,
    width: 1172,
    height: 388,
    columnWidths: [150, 224, 380, 418],
    values: methodRows,
  });
  styleTable(methods, methodRows.length, 4, 11.5);

  shape(slide, "roundRect", 54, 574, 1172, 78, C.paperAlt, C.cyanDark, 1, "rounded-lg");
  text(slide, "LINKERFORGE'S DISTINCT TASK", 78, 591, 255, 18, { fontSize: 12, bold: true, color: C.cyan, letterSpacing: 0.8 });
  text(
    slide,
    "Rank candidate linker sequences for a specified pair of fusion partners, then attach matched structural evidence. This is not equivalent to detecting natural linker residues.",
    335,
    585,
    860,
    48,
    { fontSize: 16, bold: true, color: C.ink, verticalAlignment: "middle" },
  );

  note(
    slide,
    "Presenter 1 (major)",
    "Use 'classification' precisely. DFLpred and TransDFL classify each residue as domain-linker or non-linker and provide residue-level propensity. LINKER generates candidate sequences from a PDB-derived loop library. SynLinker is a curated database and retrieval/design resource. Li and colleagues trained a supervised quantitative model on 370 randomized 18-residue linkers in one sfGFP–nanobody fusion and also evaluated high-versus-low activity classification. LinkerForge addresses another task: pair-specific candidate ranking with structural evidence. These methods are therefore complementary rather than interchangeable competitors.",
    [
      "Crasto CJ, Feng JA. LINKER: a program to generate linker sequences for fusion proteins. Protein Eng. 2000;13(5):309–312. doi:10.1093/protein/13.5.309. https://academic.oup.com/peds/article/13/5/309/1564880",
      "Liu C et al. SynLinker: an integrated system for designing linkers and synthetic fusion proteins. Bioinformatics. 2015;31(22):3700–3702. PMID:26227144. https://pubmed.ncbi.nlm.nih.gov/26227144/",
      "Meng F, Kurgan L. DFLpred: High-throughput prediction of disordered flexible linker regions in protein sequences. Bioinformatics. 2016;32(12):i341–i350. https://academic.oup.com/bioinformatics/article/32/12/i341/2289031",
      "Sharma R et al. TransDFL: Identification of disordered flexible linkers in intrinsically disordered proteins using transfer learning. Comput Struct Biotechnol J. 2023;21:496–505. PMID:36272675. https://pubmed.ncbi.nlm.nih.gov/36272675/",
      "Li et al. Multimodal deep learning-assisted rapid linker engineering for synthetic multidomain protein. 2026. PMCID:PMC12995870. https://pmc.ncbi.nlm.nih.gov/articles/PMC12995870/",
    ],
  );
}

async function replaceLiteratureEvidenceSlide(presentation, slide) {
  const inspected = await presentation.inspect({
    select: "$.slides[4]",
    include: ["id", "index", "elements.id"],
  });
  const record = inspected.records[0];
  if (!record || record.index !== 4) throw new Error("Could not locate source slide 5");
  presentation.delete((record.elements ?? []).map((element) => element.id));

  base(
    slide,
    5,
    "Published evidence",
    "Presenter 1",
    "Primary sources: LINKER (2000), SynLinker (2015), Li et al. (2026), and the CASP16 flexible-domain ensemble assessment",
  );
  slide.background.fill = "000000";
  title(
    slide,
    "Evidence from published linker studies",
    "The strength and meaning of each result depend on the task, dataset, and measurement used",
  );

  const evidenceRows = [
    ["Study", "Published finding", "What it establishes"],
    ["2000 · LINKER", "“Linker sequences are generally selected based on intuition.”", "The generator addressed heuristic sequence selection"],
    ["2015 · SynLinker", "“The linker selection is done mostly ad hoc and still based on random test in experiments.”", "The resource organized linker retrieval and fusion construction"],
    ["2026 · Li et al.", "Linker performance was context dependent in one sfGFP–nanobody system", "Measured fluorescence provides strong system-specific experimental evidence"],
    ["2026 · CASP16 ensemble assessment", "“None were close fits to the combined NMR and SAXS data.”", "No submitted ensemble reproduced both distance and orientation evidence. This does not imply that each domain misfolded"],
  ];
  const evidence = slide.tables.add({
    rows: evidenceRows.length,
    columns: 3,
    left: 54,
    top: 166,
    width: 1172,
    height: 400,
    columnWidths: [250, 520, 402],
    values: evidenceRows,
  });
  styleTable(evidence, evidenceRows.length, 3, 12.5);

  shape(slide, "line", 66, 590, 1148, 0, "none", C.cyan, 2);
  text(slide, "RESEARCH GAP", 77, 606, 150, 20, { fontSize: 12, bold: true, color: C.cyan });
  text(
    slide,
    "Choosing a sequence and predicting one static structure cannot establish the solution ensemble or biological outcome. Pair-specific measurements and ensemble-sensitive validation remain necessary.",
    229,
    599,
    972,
    48,
    { fontSize: 15, bold: true, color: C.ink },
  );

  note(
    slide,
    "Presenter 1 (major)",
    "Explain the CASP16 quotation carefully. McBride and colleagues assessed 35 population-weighted conformational ensembles from 25 predictor groups for two Staphylococcal protein A domain-linker-domain constructs. One used the conserved six-residue linker KADNKF and the other used Gly6. The team back-calculated SAXS curves and NMR residual dipolar couplings from each submitted ensemble. SAXS reports global distance and shape information, while NMR RDCs provide residue-specific orientational information between domains. No ensemble matched both types of experimental evidence closely, and none reproduced the experimentally observed difference between the two linker sequences. The result concerns distributions of relative domain poses in solution. It does not show that every individual domain was incorrectly folded.",
    [
      "Crasto CJ, Feng JA. LINKER: a program to generate linker sequences for fusion proteins. Protein Eng. 2000;13(5):309–312. doi:10.1093/protein/13.5.309.",
      "Liu C et al. SynLinker: an integrated system for designing linkers and synthetic fusion proteins. Bioinformatics. 2015;31(22):3700–3702. doi:10.1093/bioinformatics/btv447.",
      "Li et al. Multimodal deep learning-assisted rapid linker engineering for synthetic multidomain protein. Synthetic Biology. 2026. doi:10.1016/j.synbio.2026.02.003.",
      "McBride AC et al. Predicting Pose Distribution of Protein Domains Connected by Flexible Linkers Is an Unsolved Problem. Proteins. 2026;94(1):362–380. doi:10.1002/prot.70062. https://onlinelibrary.wiley.com/doi/10.1002/prot.70062",
    ],
  );
}

async function replaceCASP16ContextSlide(presentation, slide) {
  const inspected = await presentation.inspect({
    select: "$.slides[5]",
    include: ["id", "index", "elements.id"],
  });
  const record = inspected.records[0];
  if (!record || record.index !== 5) throw new Error("Could not locate source slide 6");
  presentation.delete((record.elements ?? []).map((element) => element.id));

  base(
    slide,
    6,
    "CASP16 evidence context",
    "Presenter 1",
    "McBride et al. · Proteins 94 (2026) 362–380 · doi:10.1002/prot.70062",
  );
  slide.background.fill = "000000";
  title(
    slide,
    "CASP16 tested linker-controlled domain motion",
    "A blind community assessment compared predicted pose ensembles with SAXS and NMR measurements",
  );

  const contextRows = [
    ["Assessment component", "What the investigators did", "Why it matters"],
    ["Targets", "Two Staphylococcal protein A domain-linker-domain constructs used KADNKF or GGGGGG", "Changing only the six-residue linker tested whether predictors captured linker-dependent behaviour"],
    ["Predictions", "Twenty-five groups submitted 35 population-weighted ensembles of discrete structures", "The target was a distribution of poses rather than one best PDB model"],
    ["SAXS", "Solution scattering data were measured for both linker variants", "SAXS reports global shape and interdomain distance distributions"],
    ["NMR RDC", "Residual dipolar couplings were measured for the WT KADNKF construct", "RDCs provide residue-specific information about relative domain orientation"],
    ["Evaluation and result", "The team back-calculated SAXS and RDC observables from each ensemble. None fit both closely or reproduced the WT–Gly6 SAXS difference", "A plausible domain fold can still have an incorrect interdomain pose distribution"],
  ];
  const context = slide.tables.add({
    rows: contextRows.length,
    columns: 3,
    left: 54,
    top: 166,
    width: 1172,
    height: 402,
    columnWidths: [220, 525, 427],
    values: contextRows,
  });
  styleTable(context, contextRows.length, 3, 11.7);

  shape(slide, "line", 66, 592, 1148, 0, "none", C.magenta, 2);
  text(slide, "LINKERFORGE IMPLICATION", 77, 609, 205, 20, { fontSize: 12, bold: true, color: C.magenta });
  text(
    slide,
    "Static structure confidence cannot validate flexible interdomain motion. A future dynamics claim requires ensemble sampling followed by experimental restraints such as SAXS, NMR RDC or FRET.",
    285,
    600,
    916,
    48,
    { fontSize: 15, bold: true, color: C.ink },
  );

  note(
    slide,
    "Presenter 1 (major)",
    "CASP16 introduced a blind challenge focused on continuous interdomain pose distributions. The assessors used two versions of a Staphylococcal protein A domain-linker-domain construct: the conserved KADNKF linker and Gly6. Twenty-five groups supplied 35 population-weighted ensembles. SAXS was available for both variants and reports overall shape and distance information. NMR residual dipolar couplings were available for the WT construct and report domain orientation through residue-specific measurements. The assessors back-calculated these observables from every submitted ensemble. No ensemble agreed closely with both evidence types, and no method reproduced the experimentally observed WT-versus-Gly6 SAXS difference. This separates local folding confidence from ensemble-level validation of linker-controlled motion.",
    [
      "McBride AC et al. Predicting Pose Distribution of Protein Domains Connected by Flexible Linkers Is an Unsolved Problem. Proteins. 2026;94(1):362–380. doi:10.1002/prot.70062.",
      "https://onlinelibrary.wiley.com/doi/10.1002/prot.70062",
      "CASP16 prediction and SAXS data: https://pinot3.genomecenter.ucdavis.edu/casp16/UPLOAD_fromPC/UPLOADS_ENSMBL/cleaned/",
      "NMR RDC data: BMRB accession 53002.",
    ],
  );
}

async function replaceComparisonStatusTable(presentation, slide) {
  const inspected = await presentation.inspect({
    select: "$.slides[6]",
    include: ["id", "index", "elements.id", "elements.text"],
  });
  const record = inspected.records[0];
  if (!record || record.index !== 6) throw new Error("Could not locate source slide 7");
  const sourceTable = (record.elements ?? []).find((element) => typeof element.text === "string" && element.text.startsWith("Method\tPrimary output"));
  if (!sourceTable) throw new Error("Could not locate the comparison table on source slide 7");
  presentation.delete((record.elements ?? []).map((element) => element.id));

  base(
    slide,
    7,
    "Model comparison",
    "Presenter 1",
    "Comparison by scientific task and evidence scope; no common superiority benchmark is available",
  );
  slide.background.fill = "000000";
  title(
    slide,
    "Where LinkerForge fits among existing methods",
    "The present distinction concerns workflow and evidence scope; predictive superiority has not been established",
  );

  const comparisonRows = [
    ["Method", "Primary output", "Context represented", "Evidence currently available"],
    ["LINKER / SynLinker", "Generated or retrieved linker candidates", "Linker properties and curated templates", "Generation and retrieval evidence. No matched experimental outcome benchmark"],
    ["DFLpred / APOD / TransDFL", "Per-residue DFL score", "One protein sequence", "Residue-level classification benchmarks for a different scientific task"],
    ["Li et al. (2026)", "Fluorescence prediction and linker customization", "370 linkers in one sfGFP–nanobody system plus AF2 features", "Strong system-specific experimental evidence from measured fluorescence"],
    ["AlphaFold2 / ColabFold", "Structure model and confidence", "Complete supplied sequence", "Computational confidence without functional or ensemble validation"],
    ["LinkerForge target", "Pair-conditioned ranking with a matched evidence bundle", "Both domains, attachment context, candidate, and full fusion", "Three-pair computational acceptance panel. External and experimental evidence remain incomplete"],
  ];
  const comparison = slide.tables.add({
    rows: comparisonRows.length,
    columns: 4,
    left: 54,
    top: 166,
    width: 1172,
    height: 414,
    columnWidths: [220, 250, 320, 382],
    values: comparisonRows,
  });
  styleTable(comparison, comparisonRows.length, 4, 11.5);

  shape(slide, "line", 66, 610, 1148, 0, "none", C.magenta, 2);
  text(slide, "EVIDENCE INTERPRETATION", 77, 626, 190, 20, { fontSize: 12, bold: true, color: C.magenta });
  text(
    slide,
    "Li et al. provides stronger system-specific experimental evidence. CASP16 shows why a confident single structure cannot validate flexible interdomain pose distributions. LinkerForge currently supplies computational evidence only.",
    276,
    619,
    925,
    43,
    { fontSize: 14, bold: true, color: C.ink },
  );

  note(
    slide,
    "Presenter 1 (major)",
    "Frame the comparison around evidence rather than superiority. LINKER and SynLinker provide generation or retrieval evidence. DFLpred and TransDFL provide residue-classification benchmarks for another task. Li et al. reports measured fluorescence for 370 linker variants in one sfGFP–nanobody system, which gives strong system-specific experimental evidence. The CASP16 ensemble assessment shows that computational methods still struggle to reproduce flexible interdomain pose distributions measured jointly by NMR RDC and SAXS. LinkerForge currently provides a three-pair computational acceptance panel and provenance, without equivalent functional, NMR, or SAXS evidence.",
    [
      "Crasto CJ, Feng JA. LINKER: a program to generate linker sequences for fusion proteins. Protein Eng. 2000;13(5):309–312. doi:10.1093/protein/13.5.309.",
      "Liu C et al. SynLinker: an integrated system for designing linkers and synthetic fusion proteins. Bioinformatics. 2015;31(22):3700–3702.",
      "McBride AC et al. Predicting Pose Distribution of Protein Domains Connected by Flexible Linkers Is an Unsolved Problem. Proteins. 2026;94(1):362–380. doi:10.1002/prot.70062.",
      "LinkerForge three-pair public acceptance comparison: the current panel verifies execution and traceability, not superiority.",
    ],
  );
}

function metricCard(slide, x, y, w, label, valueText, explanation, color) {
  shape(slide, "roundRect", x, y, w, 154, C.paper, C.faint, 1, "rounded-lg");
  shape(slide, "rect", x, y, 6, 154, color, "none", 0);
  text(slide, label.toUpperCase(), x + 22, y + 18, w - 42, 20, { fontSize: 12, bold: true, color });
  text(slide, valueText, x + 22, y + 48, w - 42, 35, { fontSize: 25, bold: true });
  text(slide, explanation, x + 22, y + 91, w - 42, 48, { fontSize: 13, color: C.muted });
}

async function rebuildSourceCharts(presentation, originalSlides) {
  const inspected = await presentation.inspect({
    select: "$.slides[*]",
    include: ["index", "elements.id", "elements.kind"],
  });
  const sourceChartIds = inspected.records
    .filter((record) => record.index === 17 || record.index === 18)
    .flatMap((record) => (record.elements ?? []).filter((element) => element.kind === "chart").map((element) => element.id));
  if (sourceChartIds.length !== 2) throw new Error(`Expected two source charts to rebuild; found ${sourceChartIds.length}`);
  presentation.delete(sourceChartIds);

  const regional = originalSlides[17].charts.add("bar", {
    position: { left: 667, top: 172, width: 548, height: 350 },
    title: "Regional mean pLDDT (0–100)",
    categories: ["Domain A", "Whole chain", "Domain B", "Linker"],
    series: [{ name: "Mean pLDDT", values: [91.10, 81.74, 76.90, 37.73], fill: C.cyanDark }],
    hasLegend: false,
    dataLabels: { showValue: true, position: "outEnd", textStyle: { fontSize: 12, fill: C.ink } },
    xAxis: { title: { text: "Region", textStyle: { fontSize: 12, fill: C.muted } }, textStyle: { fontSize: 12, fill: C.ink } },
    yAxis: { min: 0, max: 100, majorUnit: 20, title: { text: "pLDDT", textStyle: { fontSize: 12, fill: C.muted } }, textStyle: { fontSize: 11, fill: C.muted }, majorGridlines: { style: "solid", fill: C.faint, width: 1 } },
    titleTextStyle: { fontSize: 16, bold: true, fill: C.ink },
  });
  applyPresentationChartFont(regional, { fontFamily: FONT });

  const comparison = originalSlides[18].charts.add("bar", {
    position: { left: 57, top: 165, width: 1164, height: 398 },
    title: "Overall and linker-region pLDDT by construct",
    categories: ["GGGGS", "GS10", "(GGGGS)3", "WALDO12", "BIRD18", "EAAAK15"],
    series: [
      { name: "Whole chain", values: [87.03, 81.74, 82.75, 84.09, 80.21, 87.88], fill: C.cyanDark },
      { name: "Linker region", values: [48.57, 37.73, 38.99, 44.42, 44.87, 50.88], fill: C.magenta },
    ],
    hasLegend: true,
    legend: { position: "bottom", textStyle: { fontSize: 12, fill: C.ink } },
    dataLabels: { showValue: true, position: "outEnd", textStyle: { fontSize: 10, fill: C.ink } },
    xAxis: { textStyle: { fontSize: 11, fill: C.ink } },
    yAxis: { min: 0, max: 100, majorUnit: 20, title: { text: "pLDDT", textStyle: { fontSize: 12, fill: C.muted } }, textStyle: { fontSize: 11, fill: C.muted }, majorGridlines: { style: "solid", fill: C.faint, width: 1 } },
    titleTextStyle: { fontSize: 16, bold: true, fill: C.ink },
  });
  applyPresentationChartFont(comparison, { fontFamily: FONT });
}

function resultRows(records) {
  return [["Linker condition", "Sequence", "L", "Fusion\npLDDT", "Linker\npLDDT", "pTM", "Cross-domain\nPAE (Å)", "Displayed\nrank"], ...records.map((row) => [
    row.linker_id,
    displaySequence(row.linker_sequence),
    String(row.linker_length),
    value(row.mean_plddt),
    value(row.linker_mean_plddt),
    value(row.ptm, 2),
    value(row.mean_cross_domain_pae_angstrom),
    String(row.presentation_rank ?? row.descriptive_mean_plddt_rank),
  ])];
}

function addPairSlide(presentation, number, pair, rows, summaryPath) {
  const slide = presentation.slides.add();
  base(slide, number, "Verified public acceptance run", "Presenter 2", `Complete rank-1 results for ${pair}; all six linker-containing conditions retained`);
  title(slide, `${displayPair(pair)}: every linker condition`, `${rows[0].domain_sources.join(" + ")} · five AlphaFold2-ptm models per construct · rank 1 reported`);
  const values = resultRows(rows);
  const table = slide.tables.add({
    rows: values.length,
    columns: 8,
    left: 48,
    top: 172,
    width: 1184,
    height: 392,
    columnWidths: [168, 264, 52, 116, 116, 82, 156, 82],
    values,
  });
  styleTable(table, values.length, 8, 10.5);
  shape(slide, "roundRect", 54, 588, 1172, 66, C.paper, C.faint, 1, "rounded-lg");
  text(slide, "Reading the table", 78, 604, 175, 22, { fontSize: 14, bold: true, color: C.cyan });
  text(slide, "Displayed rank is post-hoc mean pLDDT order across the 18 linker-containing constructs. It is not a linker-success rank; lower cross-domain PAE means more confident relative placement.", 251, 599, 946, 43, { fontSize: 15, bold: true });
  note(slide, "Presenter 2 (major)", `Present all six linker-containing ${displayPair(pair)} conditions. Compare the same metrics across each row, but do not call the top pLDDT row the best biological linker.`, [summaryPath, ...rows[0].domain_sources]);
  return slide;
}

function buildReport(summary, manifest, summaryPath) {
  const sorted = [...summary.records].sort((a, b) => b.mean_plddt - a.mean_plddt);
  const highestOverall = sorted[0];
  const highestWithLinker = sorted.find((row) => row.linker_length > 0);
  const lines = [
    "# Verified public acceptance run",
    "",
    `Generated: ${summary.created_at}`,
    "",
    "## What was run",
    "",
    `- ${summary.domain_pair_count} public protein-domain pairs, ${summary.construct_count} complete fusion constructs, ${summary.models_per_construct} AlphaFold2-ptm models per construct, and ${summary.rankable_model_count} rankable models in total.`,
    `- Predictor: ${manifest.planned_settings.predictor}; MSA mode: ${manifest.planned_settings.msa_mode}; templates: ${manifest.planned_settings.templates}; seed set: ${manifest.planned_settings.seeds.join(", ")}.`,
    "- Every construct is reported. Sorting by mean pLDDT was performed only after prediction and is descriptive.",
    `- Data firewall: public sequences only = ${manifest.data_firewall.public_sequences_only}; protected benchmark labels read = ${manifest.data_firewall.protected_benchmark_labels_read}; used for training = ${manifest.data_firewall.used_for_model_training}.`,
    "- OpenMM, Rosetta, ensemble agreement, reference-domain lDDT/TM-score, automatic site mapping, DSSP, and post-OpenMM change analysis were not run in this acceptance panel.",
    "",
    "## Highest observed confidence examples",
    "",
    `- Highest whole-fusion mean pLDDT: ${displayPair(highestOverall.domain_pair)} with the ${highestOverall.linker_id} control, mean pLDDT ${value(highestOverall.mean_plddt)}, pTM ${value(highestOverall.ptm, 2)}, cross-domain PAE ${value(highestOverall.mean_cross_domain_pae_angstrom)} Å. This direct-fusion control contains no linker and is not a linker recommendation.`,
    `- Highest linker-containing mean pLDDT: ${displayPair(highestWithLinker.domain_pair)} with ${highestWithLinker.linker_id} (${highestWithLinker.linker_sequence}), mean pLDDT ${value(highestWithLinker.mean_plddt)}, linker pLDDT ${value(highestWithLinker.linker_mean_plddt)}, pTM ${value(highestWithLinker.ptm, 2)}, cross-domain PAE ${value(highestWithLinker.mean_cross_domain_pae_angstrom)} Å. The presentation uses this construct as a visible execution example, not as proof of biological superiority.`,
    "",
    "## Complete result matrix",
    "",
    "| Global rank | Protein pair | Construct | Linker | Sequence | Length | Fusion pLDDT | Domain A pLDDT | Linker pLDDT | Domain B pLDDT | pTM | Cross-domain PAE (Å) | Rank-1 model | PDB SHA-256 | Score SHA-256 |",
    "|---:|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|",
  ];
  for (const row of [...summary.records].sort((a, b) => a.descriptive_mean_plddt_rank - b.descriptive_mean_plddt_rank)) {
    lines.push(`| ${row.descriptive_mean_plddt_rank} | ${markdownCell(displayPair(row.domain_pair))} | ${markdownCell(row.construct_id)} | ${markdownCell(row.linker_id)} | ${markdownCell(displaySequence(row.linker_sequence))} | ${row.linker_length} | ${value(row.mean_plddt)} | ${value(row.domain_a_mean_plddt)} | ${value(row.linker_mean_plddt)} | ${value(row.domain_b_mean_plddt)} | ${value(row.ptm, 2)} | ${value(row.mean_cross_domain_pae_angstrom)} | ${markdownCell(row.rank_1_model)} | \`${row.pdb_sha256}\` | \`${row.score_json_sha256}\` |`);
  }
  lines.push(
    "", "## Metric interpretation", "",
    `- **Fusion pLDDT:** ${summary.metric_interpretation.mean_plddt}`,
    `- **Linker pLDDT:** ${summary.metric_interpretation.linker_mean_plddt}`,
    `- **pTM:** ${summary.metric_interpretation.ptm}`,
    `- **Cross-domain PAE:** ${summary.metric_interpretation.mean_cross_domain_pae_angstrom}`,
    "", "## Provenance and reproducibility", "",
    `- Summary: \`${summaryPath}\``,
    `- Panel manifest: \`${summary.provenance.panel_manifest}\` (SHA-256 \`${summary.provenance.panel_manifest_sha256}\`)`,
    `- Existing reference-pair summary: \`${summary.provenance.existing_summary}\` (SHA-256 \`${summary.provenance.existing_summary_sha256}\`)`,
    `- Repository commit recorded by summarizer: \`${summary.provenance.commit_sha}\``,
    `- Summary schema: \`${summary.schema_version}\``,
    "- Each result row retains the exact rank-1 PDB and score-JSON path plus SHA-256 hashes.",
    "", "## Claim boundary", "", summary.claim_boundary,
    "", "## Limitations", "", ...summary.limitations.map((item) => `- ${item}`),
    "", "## Presentation rule", "",
    "The presentation may show the highest observed confidence example only as a transparent implementation example. It must not describe that construct as the best linker, a compatibility probability, a biological winner, or an experimentally validated fusion.",
    "",
  );
  return lines.join("\n");
}

async function main() {
  const input = parseArgs();
  if (await pathExists(input.output)) throw new Error(`Refusing to overwrite completed presentation: ${input.output}`);
  if (await pathExists(input.report)) throw new Error(`Refusing to overwrite completed acceptance report: ${input.report}`);

  const summary = JSON.parse(await fs.readFile(input.summary, "utf8"));
  if (summary.status !== "complete" || summary.construct_count !== 21 || summary.rankable_model_count !== 105) {
    throw new Error("The presentation requires a complete 21-construct, 105-model summary");
  }
  if (summary.records.some((row) => typeof row.linker_sequence !== "string")) {
    throw new Error("Every acceptance record must retain its exact linker sequence");
  }
  const manifest = JSON.parse(await fs.readFile(summary.provenance.panel_manifest, "utf8"));
  if (manifest.record_count !== 14 || manifest.data_firewall.protected_benchmark_labels_read !== false) {
    throw new Error("Unexpected multi-pair manifest or data-firewall state");
  }
  const linkerRecords = [...summary.records]
    .filter((row) => row.linker_length > 0)
    .sort((a, b) => b.mean_plddt - a.mean_plddt)
    .map((row, index) => ({ ...row, presentation_rank: index + 1 }));
  if (linkerRecords.length !== 18) throw new Error(`Expected 18 linker-containing constructs; found ${linkerRecords.length}`);

  await fs.mkdir(BUILD, { recursive: true });
  await fs.mkdir(PREVIEW_DIR, { recursive: true });
  await fs.mkdir(path.dirname(input.output), { recursive: true });
  await fs.mkdir(path.dirname(input.report), { recursive: true });

  const presentation = await PresentationFile.importPptx(await FileBlob.load(input.base));
  if (presentation.slides.items.length !== 25) throw new Error(`Expected 25 source slides; found ${presentation.slides.items.length}`);
  const originalSlides = [...presentation.slides.items];

  await replaceExact(presentation, {
    "Computational prototype · verified results only": "Computational prototype · matched public-sequence evidence",
    "Four-presenter project briefing": "Four-presenter research briefing",
    "Kushagra Srivastrava , Hardik Tyagi\nRam , Kritharith Shankar": "Kushagra Srivastava · Hardik Tyagi\nRam · Kritharith Shankar",
    "PUBLIC ACCEPTANCE CONSTRUCT": "PUBLIC ACCEPTANCE PANEL · 3 PAIRS",
    "A 32-minute briefing with clear handoffs": "A 40-minute briefing with clear handoffs",
    "13 min  ·  Slides 13–20": "19 min  ·  Slides 13–20, 23–30",
    "4 min  ·  Slides 21–25": "6 min  ·  Slides 21–22, 31–34",
    "Verified public acceptance run": "Verified acceptance run: reference pair",
    "Matched structural sensitivity panel": "Reference-pair structural sensitivity panel",
    "Research gap and present scope": "Unresolved research gap",
    "Research question and evidence stages": "Research questions and evidence",
    "LinkerForge workflow and provenance": "Traceable workflow from sequence design to evidence",
    "Proposed multimodal compatibility ranker": "Proposed pair-conditioned compatibility model",
    "How the ranker would learn and make a decision": "Proposed training and evaluation logic",
    "Dataset acquisition and evidence tiers": "Dataset composition and permitted uses",
    "Independent annotation and adjudication": "Annotation protocol and adjudication",
    "Validation layers and the claim firewall": "Validation layers that protect the scientific claim",
    "Candidate generation and ranking status": "Implemented candidate controls and future ranking",
    "Execution evidence by tool": "Evidence contributed by each executed tool",
    "Design and structure review in one workspace": "Research interface for design and structural review",
    "Publication readiness": "Publication readiness and unresolved evidence",
    "Verified scope and next scientific milestone": "Current evidence and next scientific milestone",
    "Direct": "Short GS",
    "no linker": "GGGGS",
    "0 aa control": "5 aa sensitivity",
    "The protected DFL benchmark remains fail-closed. Published or spent labels cannot guide new training, thresholds, filtering, or feature design.": "Training is blocked to prevent leakage. The DP94 screen links 88 of 98 proteins to existing data, and CAID3 results have already been viewed. A new sealed, family-independent test set is required.",
    "142-aa public acceptance + 7-construct panel": "3 public pairs · 18 linker constructs",
    "35 rankable models": "90 rankable models",
    "Presenter 4  ·  23": "Presenter 4  ·  32",
    "Presenter 4  ·  24": "Presenter 4  ·  33",
    "Presenter 4  ·  25": "Presenter 4  ·  34",
    "LinkerForge links exact candidate sequences to matched structural evidence and provenance.": "LinkerForge links exact candidate sequences to matched structural evidence and provenance across three public pairs.",
    "Verified local runs establish execution, not biological success or optimality.": "The 18-construct linker panel verifies execution and traceability, not biological success or optimality.",
  });
  await replaceSpeakerSummaryWithLinkerDefinition(presentation, originalSlides[1]);
  await replaceResearchProblemSlide(presentation, originalSlides[2]);
  await replacePublishedMethodsSlide(presentation, originalSlides[3]);
  await replaceLiteratureEvidenceSlide(presentation, originalSlides[4]);
  await replaceCASP16ContextSlide(presentation, originalSlides[5]);
  await replaceComparisonStatusTable(presentation, originalSlides[6]);
  await rebuildSourceCharts(presentation, originalSlides);
  originalSlides[18].speakerNotes.text = originalSlides[18].speakerNotes.text.replace(
    "The direct fusion has no linker-region value and is plotted as zero only for chart alignment. Do not compare it as if zero were a measured linker pLDDT. The panel is matched for sequence context and prediction settings, but one domain pair cannot establish generalization or superiority.",
    "Only linker-containing conditions are shown. The panel uses matched sequence context and prediction settings, but one domain pair cannot establish generalization or superiority."
  );

  originalSlides[17].speakerNotes.text += `\n\nUpdated acceptance context: this slide is the original Ubiquitin–GB1 reference-pair run. Slides 23–31 report the complete three-pair panel. Source: ${input.summary}`;
  originalSlides[14].speakerNotes.text += "\n\nWhy the protected gate blocks training: this is a data-independence failure, not a software crash. The DP94 candidate screen connected 88 of 98 proteins to existing train, development, or auxiliary proteins. The public CAID3 common set was already viewed in a one-time comparison, so it cannot be reused to choose a new model. The self-contained B1 V2 split passes its internal checks, but it does not create a fresh, untouched external benchmark.";
  originalSlides[19].speakerNotes.text += `\n\nUpdated presentation count: 18 linker-containing constructs and 90 rankable AlphaFold2-ptm models. Source: ${input.summary}`;
  originalSlides[24].speakerNotes.text += `\n\nUpdated conclusion evidence: ${input.summary}`;

  const inserted = [];
  const slide23 = presentation.slides.add();
  inserted.push(slide23);
  base(slide23, 23, "Verified public acceptance run", "Presenter 2", "Public-sequence panel manifest + versioned comparison schema; protected labels were not read");
  title(slide23, "Acceptance run design and verification", "Matched settings across three short public domain pairs; every linker-containing condition is retained");
  metricCard(slide23, 54, 186, 265, "Protein pairs", "3 public pairs", "Ubiquitin–GB1, crambin–BPTI, and protein G–repressor", C.blue);
  metricCard(slide23, 335, 186, 265, "Constructs", "18 linker fusions", "Six matched linker conditions per pair", C.cyan);
  metricCard(slide23, 616, 186, 265, "Structure models", "90 rankable models", "Five AlphaFold2-ptm models per construct; rank 1 summarized", C.magenta);
  metricCard(slide23, 897, 186, 329, "Firewall", "Protected labels: not read", "Public sequences only; no training or threshold selection", C.green);
  const gateRows = [
    ["Gate", "Recorded evidence", "Outcome"],
    ["Matched predictor settings", "AF2-ptm · seed 0 · templates off · same MSA mode", "Pass"],
    ["Completeness", "18/18 linker constructs · 90/90 rankable models", "Pass"],
    ["Sequence identity", "Rank-1 PDB sequence checked against immutable input", "Pass"],
    ["Traceability", "PDB + score JSON paths and SHA-256 per result", "Pass"],
    ["Claim firewall", "No compatibility, function, stability, or optimality label", "Enforced"],
  ];
  const gates = slide23.tables.add({ rows: gateRows.length, columns: 3, left: 54, top: 373, width: 1172, height: 250, columnWidths: [255, 700, 217], values: gateRows });
  styleTable(gates, gateRows.length, 3, 12);
  note(slide23, "Presenter 2 (major)", "Start the expanded acceptance section with what was actually run and what was deliberately excluded. The completion count verifies pipeline execution and traceability only; it does not validate a compatibility ranker.", [input.summary, summary.provenance.panel_manifest]);

  const slide24 = presentation.slides.add();
  inserted.push(slide24);
  const highestWithLinker = linkerRecords[0];
  base(slide24, 24, "Verified public acceptance run", "Presenter 2", "Rank-1 AlphaFold2-ptm model; highest whole-fusion mean pLDDT among the 18 linker-containing constructs");
  title(slide24, "Highest-confidence linker-containing construct", `${displayPair(highestWithLinker.domain_pair)} with ${highestWithLinker.linker_id}: ${highestWithLinker.linker_sequence}`);
  slide24.images.add({
    blob: await imageBytes(HIGHEST_CONFIDENCE_IMAGE),
    contentType: "image/png",
    alt: `3Dmol rendering of ${displayPair(highestWithLinker.domain_pair)} joined by ${highestWithLinker.linker_id}; Domain A blue, linker magenta, Domain B red`,
    fit: "contain",
    position: { left: 48, top: 162, width: 755, height: 432 },
  });
  text(slide24, "DOMAIN A", 104, 579, 110, 20, { fontSize: 12, bold: true, color: C.blue });
  text(slide24, "LINKER", 310, 579, 100, 20, { fontSize: 12, bold: true, color: C.magenta });
  text(slide24, "DOMAIN B", 520, 579, 110, 20, { fontSize: 12, bold: true, color: C.red });
  text(slide24, value(highestWithLinker.mean_plddt), 850, 177, 300, 58, { fontSize: 46, bold: true, color: C.cyan });
  text(slide24, "Fusion mean pLDDT", 850, 231, 310, 24, { fontSize: 16, bold: true, color: C.muted });
  text(slide24, `Linker pLDDT  ${value(highestWithLinker.linker_mean_plddt)}`, 850, 286, 340, 30, { fontSize: 22, bold: true, color: C.magenta });
  text(slide24, `pTM  ${value(highestWithLinker.ptm, 2)}`, 850, 334, 340, 30, { fontSize: 22, bold: true, color: C.blue });
  text(slide24, `Cross-domain PAE  ${value(highestWithLinker.mean_cross_domain_pae_angstrom)} Å`, 850, 382, 360, 30, { fontSize: 22, bold: true, color: C.yellow });
  text(slide24, "Displayed rank 1 of 18 linker-containing constructs", 850, 442, 345, 45, { fontSize: 16, bold: true });
  shape(slide24, "roundRect", 830, 506, 380, 105, C.paper, C.faint, 1, "rounded-lg");
  text(slide24, "Public inputs: PDB 1CRN chain A and PDB 4PTI chain A. Linker: EAAAKEAAAKEAAAK. Complete fusion: 119 residues.", 850, 523, 340, 72, { fontSize: 14, color: C.muted, bold: true });
  text(slide24, "Model confidence supports structural inspection. It does not prove expression, activity, stability, or linker superiority.", 58, 625, 1148, 34, { fontSize: 16, bold: true, color: C.yellow, alignment: "center" });
  note(slide24, "Presenter 2 (major)", "Show this structure as the clearest linker-containing execution example. Explain that this is a descriptive confidence result selected after prediction, not an independently validated winner.", [input.summary, highestWithLinker.pdb, highestWithLinker.score_json, "https://www.nature.com/articles/s41586-021-03819-2"]);

  const grouped = Object.fromEntries([...new Set(linkerRecords.map((row) => row.domain_pair))].sort().map((pair) => [pair, linkerRecords.filter((row) => row.domain_pair === pair)]));
  for (const [offset, pair] of ["1UBQ_2GB1", "1CRN_4PTI", "1PGA_1LMB"].entries()) {
    if (!grouped[pair] || grouped[pair].length !== 6) throw new Error(`Expected six linker-containing records for ${pair}`);
    inserted.push(addPairSlide(presentation, 26 + offset, pair, grouped[pair], input.summary));
  }

  const slide25 = presentation.slides.add();
  inserted.push(slide25);
  base(slide25, 25, "Matched public panel", "Presenter 2", "Rank-1 AlphaFold2-ptm confidence; descriptive top linker-containing row per pair; all 18 linker constructs retained");
  title(slide25, "Structure-confidence comparison across three protein pairs", "The chart uses the highest linker-containing mean pLDDT observed for each pair");
  const perPair = ["1CRN_4PTI", "1PGA_1LMB", "1UBQ_2GB1"].map((pair) => linkerRecords.find((row) => row.domain_pair === pair));
  if (perPair.some((row) => !row)) throw new Error("Missing a linker-containing top row for one or more domain pairs");
  const chart = slide25.charts.add("bar", {
    position: { left: 56, top: 172, width: 1168, height: 365 },
    title: "Top observed linker-containing construct per pair",
    categories: perPair.map((row) => displayPair(row.domain_pair)),
    series: [
      { name: "Whole fusion pLDDT", values: perPair.map((row) => chartValue(row.mean_plddt)), fill: C.cyanDark },
      { name: "Linker pLDDT", values: perPair.map((row) => chartValue(row.linker_mean_plddt)), fill: C.magenta },
    ],
    hasLegend: true,
    legend: { position: "bottom", textStyle: { fontSize: 12, fill: C.ink } },
    dataLabels: { showValue: true, position: "outEnd", textStyle: { fontSize: 11, fill: C.ink } },
    xAxis: { textStyle: { fontSize: 12, fill: C.ink } },
    yAxis: { min: 0, max: 100, majorUnit: 20, title: { text: "pLDDT", textStyle: { fontSize: 12, fill: C.muted } }, textStyle: { fontSize: 11, fill: C.muted }, majorGridlines: { style: "solid", fill: C.faint, width: 1 } },
    titleTextStyle: { fontSize: 16, bold: true, fill: C.ink },
  });
  applyPresentationChartFont(chart, { fontFamily: FONT });
  shape(slide25, "roundRect", 56, 564, 1168, 82, C.paper, C.faint, 1, "rounded-lg");
  text(slide25, "Interpretation", 80, 579, 150, 24, { fontSize: 15, bold: true, color: C.yellow });
  text(slide25, "These values compare model confidence under matched settings. They do not measure compatibility, stability, expression, or activity.", 230, 578, 960, 49, { fontSize: 17, bold: true });
  note(slide25, "Presenter 2 (major)", "Compare the same two confidence summaries across the three public pairs. Each bar is the top linker-containing row observed after prediction. Slides 26 to 28 retain every linker-containing condition.", [input.summary, "https://www.nature.com/articles/s41586-021-03819-2", "https://www.nature.com/articles/s41592-022-01488-1"]);

  const slide29 = presentation.slides.add();
  inserted.push(slide29);
  base(slide29, 29, "Interpretation", "Presenter 2", "Metric definitions are taken from the versioned comparison JSON and AlphaFold/ColabFold output semantics");
  title(slide29, "Interpreting AlphaFold and ColabFold confidence", "All four values come from ColabFold/AlphaFold2-ptm outputs, not from the proposed LinkerForge ranker");
  metricCard(slide29, 54, 184, 560, "Fusion mean pLDDT", "Local-coordinate confidence", "Average over the complete fusion. High values do not establish function, expression, stability, or correct domain orientation.", C.cyan);
  metricCard(slide29, 666, 184, 560, "Linker mean pLDDT", "Inserted-residue confidence", "Low values may be consistent with flexibility, but they do not prove flexibility or incompatibility.", C.magenta);
  metricCard(slide29, 54, 366, 560, "pTM", "Global topology confidence", "A model-level AlphaFold estimate. It is not a linker-specific success probability.", C.blue);
  metricCard(slide29, 666, 366, 560, "Cross-domain PAE", "Relative-placement uncertainty", "Mean bidirectional Domain A↔Domain B PAE. Lower values imply more confident relative placement, not validated activity.", C.yellow);
  shape(slide29, "roundRect", 54, 558, 1172, 86, C.paper, C.red, 1, "rounded-lg");
  text(slide29, "Permitted claim", 78, 576, 160, 22, { fontSize: 15, bold: true, color: C.red });
  text(slide29, summary.claim_boundary, 240, 570, 955, 51, { fontSize: 17, bold: true });
  note(slide29, "Presenter 2 (major)", "This is the slide to use if someone asks which model produced the confidence. Every displayed confidence metric comes from the AlphaFold2-ptm/ColabFold structure prediction. The proposed compatibility ranker remains scientifically untrained.", [input.summary, "https://www.nature.com/articles/s41586-021-03819-2", "https://www.nature.com/articles/s41592-022-01488-1"]);

  const slide30 = presentation.slides.add();
  inserted.push(slide30);
  base(slide30, 30, "Reproducibility", "Presenter 2", "Versioned JSON + CSV + exact PDB/score files + SHA-256 hashes; generated only after full-panel completion");
  title(slide30, "Traceability from input sequence to reported figure", "The run artifacts remain the scientific record; this deck presents a selected view");
  const trail = [
    ["01", "Immutable inputs", "panel.fasta + panel_manifest.json", "Sequence hashes, public sources, matched settings"],
    ["02", "Displayed predictions", "90 score JSON files + 90 PDB models", "Five models per linker construct; rank-1 path retained"],
    ["03", "Versioned summary", "three_pair_comparison.json + .csv", "All 18 displayed linker rows, metrics, limitations"],
    ["04", "Integrity", "Per-row PDB and score SHA-256", "Detects file substitution after reporting"],
    ["05", "Presentation", "Editable charts, tables, and speaker notes", "All numbers derived from the JSON summary"],
  ];
  for (const [i, [number, heading, artifact, meaning]] of trail.entries()) {
    const y = 172 + i * 89;
    shape(slide30, "roundRect", 54, y, 1172, 72, i % 2 ? C.paperAlt : C.paper, C.faint, 1, "rounded-lg");
    shape(slide30, "ellipse", 72, y + 15, 42, 42, i === 4 ? C.magenta : C.cyan, "none", 0);
    text(slide30, number, 72, y + 23, 42, 22, { fontSize: 12, bold: true, color: C.bg, alignment: "center" });
    text(slide30, heading, 136, y + 13, 215, 24, { fontSize: 17, bold: true });
    text(slide30, artifact, 360, y + 12, 430, 26, { fontSize: 15, bold: true, typeface: MONO, color: C.cyan });
    text(slide30, meaning, 807, y + 11, 388, 46, { fontSize: 14, color: C.muted });
  }
  note(slide30, "Presenter 2 (major)", "Use this slide to demonstrate provenance. The presentation shows the linker-containing matrix; the JSON remains the machine-readable authority; the CSV supports analysis; hashes bind each displayed row to its PDB and score JSON.", [input.summary, input.report, summary.provenance.panel_manifest]);

  const slide31 = presentation.slides.add();
  inserted.push(slide31);
  base(slide31, 31, "Future structural validation", "Presenter 4", "Deferred by decision; not executed by this acceptance run; full protocols are documented in the future roadmap");
  title(slide31, "Deferred structural analyses", "These methods add orthogonal evidence later; none can retroactively turn this panel into biological validation");
  const futureRows = [
    ["Deferred method", "Why it is needed", "How it would be run later"],
    ["Ensemble agreement", "Separates one-model confidence from repeatability", "Compare seeds/models and report residue- and domain-level dispersion"],
    ["Reference lDDT / TM-score", "Tests domain preservation beyond RMSD", "Align each predicted domain to the sequence-matched PDB reference"],
    ["Automatic functional-site mapping", "Checks whether protected sites remain spatially accessible", "Map curated residues; calculate exposure, occlusion, and linker proximity"],
    ["DSSP", "Tracks secondary-structure changes near junctions", "Run DSSP on references and predictions; compare junction windows"],
    ["Post-OpenMM change analysis", "Measures drift from the starting prediction", "Compare minimized/replicated trajectories with the same structural metrics"],
    ["Production OpenMM", "Tests longer-timescale mechanical stability", "Pre-register protocol; run replicated GPU trajectories only after resource approval"],
  ];
  const future = slide31.tables.add({ rows: futureRows.length, columns: 3, left: 50, top: 164, width: 1180, height: 454, columnWidths: [260, 420, 500], values: futureRows });
  styleTable(future, futureRows.length, 3, 10.5);
  shape(slide31, "roundRect", 50, 631, 1180, 32, C.paper, C.yellow, 1, "rounded-lg");
  text(slide31, "OpenMM remains on hold. Rosetta and independent predictor families also remain future integrations; no hidden substitute is claimed.", 70, 637, 1140, 20, { fontSize: 13, bold: true, alignment: "center" });
  note(slide31, "Presenter 4 (minor)", "Explain each deferred method as a future evidence layer, not a missing red failure. Production OpenMM is intentionally on hold. Rosetta requires a licensed installation, and other structure-model families require pinned matched environments.", [path.join(ROOT, "docs", "FUTURE_STRUCTURAL_VALIDATION_ROADMAP_2026-10-21.md")]);

  const slideInventory = await presentation.inspect({ select: "$.slides[*]", include: ["id", "index"] });
  const slideIds = slideInventory.records.map((record) => record.id);
  if (slideIds.length !== 34) throw new Error(`Expected 34 slides before reorder; found ${slideIds.length}`);
  const acceptance = slideIds.slice(25);
  presentation.slides.reorder([
    ...slideIds.slice(0, 22),
    acceptance[0],
    acceptance[1],
    acceptance[5],
    acceptance[2],
    acceptance[3],
    acceptance[4],
    acceptance[6],
    acceptance[7],
    acceptance[8],
    ...slideIds.slice(22, 25),
  ]);

  for (const slide of presentation.slides.items) {
    slide.background.fill = "000000";
  }

  await fs.writeFile(CANDIDATE_REPORT, buildReport(summary, manifest, input.summary), "utf8");
  await (await PresentationFile.exportPptx(presentation)).save(CANDIDATE);
  const outputStem = path.basename(input.output, path.extname(input.output)).replace(/[^A-Za-z0-9._-]+/g, "-");
  const result = await finalizePresentation({
    explicitTotalSlideCount: 34,
    workspaceDir: ROOT,
    candidatePath: CANDIDATE,
    finalPath: input.output,
    pythonExecutable: PYTHON,
    integrityValidatorPath: path.join(SKILL, "container_tools", "inspect_presentation_package_integrity.py"),
    layoutValidatorPath: path.join(SKILL, "container_tools", "inspect_presentation_layout_geometry.py"),
    layoutArgs: [
      "--expected-slide-size-emu", "12192000,6858000", "--validate-bullet-geometry", "--validate-heading-fit",
      ...[4, 5, 6, 7, 13, 22, 23, 26, 27, 28, 31].flatMap((number) => ["--require-native-table-slide", String(number)]),
    ],
    requiredNativeTableOwnerSlides: [4, 5, 6, 7, 13, 22, 23, 26, 27, 28, 31],
    requiredNativeChartOwnerSlides: [18, 19, 25],
    materializeLiteralChartWorkbooks: true,
    fontPolicy: { basis: "reference", families: [FONT, MONO], referencePath: input.base, referenceSha256: await sha256File(input.base) },
    verifyArtifactToolImport: true,
    receiptPath: path.join(BUILD, `${outputStem}.final.validation.json`),
  });

  await fs.copyFile(CANDIDATE_REPORT, input.report);
  for (const [index, slide] of presentation.slides.items.entries()) {
    const png = await presentation.export({ slide, format: "png", scale: 1 });
    await fs.writeFile(path.join(PREVIEW_DIR, `slide-${String(index + 1).padStart(2, "0")}.png`), new Uint8Array(await png.arrayBuffer()));
  }
  const montage = await presentation.export({ format: "webp", montage: true });
  const montagePath = path.join(BUILD, "linker-only-montage.webp");
  await fs.writeFile(montagePath, new Uint8Array(await montage.arrayBuffer()));

  const progress = JSON.parse(await fs.readFile(input.progress, "utf8"));
  Object.assign(progress, {
    state: "complete",
    overall_percent: 100.0,
    acceptance_report: input.report,
    presentation: input.output,
    presentation_sha256: await sha256File(input.output),
    presentation_validation: result,
    presentation_montage: montagePath,
    updated_utc: new Date().toISOString(),
  });
  await atomicJson(input.progress, progress);
  process.stdout.write(`${JSON.stringify({ output: input.output, report: input.report, slides: 34, result })}\n`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
