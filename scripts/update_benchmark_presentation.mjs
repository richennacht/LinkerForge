#!/usr/bin/env node
/** Append verified multi-pair benchmark and evidence-scope slides to the project deck. */

import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

const ROOT = "M:\\ChatGPT\\LinkerForge 2";
const SKILL = "C:\\Users\\HARDIK\\.codex\\plugins\\cache\\openai-primary-runtime\\presentations\\26.1020.11127\\skills\\presentations";
const PYTHON = "C:\\Users\\HARDIK\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe";
const DEFAULT_BASE = path.join(ROOT, "deliverables", "LinkerForge_Detailed_Architecture_Data_Provenance_4_Presenters_Final.pptx");
const DEFAULT_FINAL = path.join(ROOT, "deliverables", "LinkerForge_Detailed_Architecture_Data_Provenance_4_Presenters_Multi_Protein_Benchmark.pptx");
const BUILD = path.join(ROOT, ".codex", "presentation-build", "multi-protein-benchmark");
const CANDIDATE = path.join(BUILD, "candidate-multi-protein-benchmark.pptx");
const PREVIEW_DIR = path.join(BUILD, "previews");
const ARTIFACT_TOOL = "C:\\Users\\HARDIK\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\node\\node_modules\\@oai\\artifact-tool\\dist\\artifact_tool.mjs";

const { FileBlob, PresentationFile } = await import(pathToFileURL(ARTIFACT_TOOL).href);

const utils = await import(pathToFileURL(path.join(SKILL, "container_tools", "artifact_tool_utils.mjs")).href);
const { applyPresentationChartFont, finalizePresentation } = utils;

const C = {
  bg: "091A15", paper: "102720", ink: "F4F7F5", muted: "AFC0B8", faint: "29433A",
  cyan: "3FC6D6", cyanDark: "167D89", magenta: "E84AAB", green: "2DBD73",
  yellow: "F1C84B", red: "E46760", blue: "4E9BEC", violet: "A88CEB",
};
const FONT = "Arial";
const MONO = "Cascadia Mono";
const W = 1280;
const H = 720;

function args() {
  const parsed = {};
  for (let i = 2; i < process.argv.length; i += 2) parsed[process.argv[i].replace(/^--/, "")] = process.argv[i + 1];
  if (!parsed.summary || !parsed.progress) throw new Error("Required: --summary PATH --progress PATH [--base PATH] [--output PATH]");
  return {
    summary: path.resolve(parsed.summary),
    progress: path.resolve(parsed.progress),
    base: path.resolve(parsed.base ?? DEFAULT_BASE),
    output: path.resolve(parsed.output ?? DEFAULT_FINAL),
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

function styleTable(table, rows, columns, widths) {
  table.cells.block({ row: 0, column: 0, rowCount: 1, columnCount: columns }).assign({
    fill: C.cyanDark,
    textStyle: { typeface: FONT, fontSize: 13, bold: true, color: "FFFFFF" },
    verticalAlignment: "middle",
  });
  if (rows > 1) table.cells.block({ row: 1, column: 0, rowCount: rows - 1, columnCount: columns }).assign({
    fill: C.paper,
    textStyle: { typeface: FONT, fontSize: widths ? 11 : 12, color: C.ink },
    verticalAlignment: "middle",
  });
  table.borders.assign({ style: "solid", fill: C.faint, width: 1 });
}

function displayPair(pair) {
  return ({ "1UBQ_2GB1": "Ubiquitin–GB1", "1CRN_4PTI": "Crambin–BPTI", "1PGA_1LMB": "Protein G–repressor" })[pair] ?? pair;
}

function value(number, digits = 1) {
  return number === null || number === undefined ? "—" : Number(number).toFixed(digits);
}

async function atomicJson(file, object) {
  const temporary = `${file}.tmp`;
  await fs.writeFile(temporary, `${JSON.stringify(object, null, 2)}\n`, "utf8");
  await fs.rename(temporary, file);
}

async function main() {
  const input = args();
  const summary = JSON.parse(await fs.readFile(input.summary, "utf8"));
  if (summary.status !== "complete" || summary.construct_count !== 21 || summary.rankable_model_count !== 105) {
    throw new Error("The presentation requires a complete 21-construct, 105-model summary");
  }
  await fs.mkdir(BUILD, { recursive: true });
  await fs.mkdir(PREVIEW_DIR, { recursive: true });
  await fs.mkdir(path.dirname(input.output), { recursive: true });

  const presentation = await PresentationFile.importPptx(await FileBlob.load(input.base));
  if (presentation.slides.items.length !== 25) throw new Error(`Expected 25 source slides; found ${presentation.slides.items.length}`);

  const perPair = summary.descriptive_top_per_domain_pair;
  const slide26 = presentation.slides.add();
  base(slide26, 26, "Multi-protein benchmark", "Presenter 2", "Durable public panel summary; rank-1 AlphaFold2-ptm confidence; all 21 constructs retained");
  title(slide26, "Matched public panel across three protein pairs", "21 complete fusions · seven linker conditions per pair · 105 rankable AlphaFold2-ptm models");
  const chart = slide26.charts.add("bar", {
    position: { left: 56, top: 172, width: 1168, height: 365 },
    title: "Highest-confidence linker-containing construct per pair (descriptive)",
    categories: perPair.map((row) => displayPair(row.domain_pair)),
    series: [
      { name: "Whole fusion pLDDT", values: perPair.map((row) => row.mean_plddt), fill: C.cyanDark },
      { name: "Linker pLDDT", values: perPair.map((row) => row.linker_mean_plddt), fill: C.magenta },
    ],
    hasLegend: true,
    legend: { position: "bottom", textStyle: { typeface: FONT, fontSize: 12, color: C.ink } },
    dataLabels: { showValue: true, position: "outEnd", textStyle: { typeface: FONT, fontSize: 11, color: C.ink } },
    xAxis: { textStyle: { typeface: FONT, fontSize: 12, color: C.ink } },
    yAxis: { minimumScale: 0, maximumScale: 100, majorUnit: 20, title: { text: "pLDDT", textStyle: { typeface: FONT, fontSize: 12, color: C.muted } }, textStyle: { typeface: FONT, fontSize: 11, color: C.muted }, majorGridlines: { style: "solid", fill: C.faint, width: 1 } },
    titleTextStyle: { typeface: FONT, fontSize: 16, bold: true, color: C.ink },
  });
  applyPresentationChartFont(chart, { fontFamily: FONT });
  shape(slide26, "roundRect", 56, 564, 1168, 82, C.paper, C.faint, 1, "rounded-lg");
  text(slide26, "Interpretation boundary", 80, 579, 212, 24, { fontSize: 15, bold: true, color: C.yellow });
  text(slide26, "These are post-hoc AlphaFold confidence summaries, not compatibility probabilities, biological winners, or evidence that one linker is generally superior.", 289, 578, 910, 49, { fontSize: 17, bold: true });
  note(slide26, "Presenter 2 (major)", "This extends the prior single-pair panel to three deliberately short public pairs while keeping settings matched. Explain that pLDDT measures predicted local coordinate confidence. We retain every construct and display a descriptive top linker-containing row per pair only for implementation demonstration.", [input.summary, "https://www.nature.com/articles/s41586-021-03819-2", "https://www.nature.com/articles/s41592-022-01488-1"]);

  const top = summary.descriptive_top_by_mean_plddt;
  const slide27 = presentation.slides.add();
  base(slide27, 27, "Descriptive ranking", "Presenter 2", "Sorted after prediction by rank-1 whole-fusion pLDDT; not an independent test or biological selection rule");
  title(slide27, "Highest structure-confidence examples", "The complete CSV and JSON retain every construct, metric, file hash, setting, and limitation");
  const rows = [["Rank", "Protein pair", "Linker", "Whole\npLDDT", "Linker\npLDDT", "pTM", "Cross-domain\nPAE (Å)"], ...top.map((row) => [
    String(row.descriptive_mean_plddt_rank), displayPair(row.domain_pair), row.linker_id,
    value(row.mean_plddt), value(row.linker_mean_plddt), value(row.ptm, 2), value(row.mean_cross_domain_pae_angstrom),
  ])];
  const table = slide27.tables.add({ rows: rows.length, columns: 7, left: 54, top: 169, width: 1172, height: 386, columnWidths: [70, 206, 180, 134, 134, 104, 190], values: rows });
  styleTable(table, rows.length, 7, true);
  shape(slide27, "roundRect", 54, 581, 1172, 68, C.paper, C.faint, 1, "rounded-lg");
  text(slide27, "Use in the demo", 78, 597, 170, 24, { fontSize: 15, bold: true, color: C.cyan });
  text(slide27, "Show the highest-confidence public example as an execution result, then show its linker pLDDT, pTM, and cross-domain PAE so uncertainty remains visible.", 248, 594, 946, 43, { fontSize: 16, bold: true });
  note(slide27, "Presenter 2 (major)", "Read the first row as the highest mean pLDDT observed in this small public panel, not as the best linker. Point out the associated linker confidence, global pTM, and cross-domain PAE. This is transparent example selection for a presentation, not prospective validation.", [input.summary]);

  const slide28 = presentation.slides.add();
  base(slide28, 28, "Evidence alternatives", "Presenter 4", "Alternatives narrow specific evidence gaps; they do not convert missing validation into a positive claim");
  title(slide28, "Implemented alternatives and future method scope", "Red is reserved for an actual failed check; unavailable methods remain future integrations");
  const altRows = [
    ["Evidence need", "Available project-level evidence", "Status / claim boundary"],
    ["Rosetta loop closure / clash", "Cα break-overlap screen + OpenMM all-atom minimization", "Implemented alternative; not Rosetta KIC sampling"],
    ["Long-timescale stability", "3 × 50 ps CUDA OpenMM replicas", "Execution smoke test; production MD is future scope"],
    ["External linker detectors", "DFLpred public sensitivity panel", "Implemented comparator; not compatibility validation"],
    ["Independent predictor families", "ESMFold / OpenFold / Boltz integration plan", "Future pinned, matched sensitivity analysis"],
    ["Linker-specific external tools", "LINKER-Pred / TransDFL / DisoFLAG protocol", "Future environment + firewall-safe target set"],
    ["Biological success", "No computational substitute", "Outside current scope; no success claim"],
  ];
  const alt = slide28.tables.add({ rows: altRows.length, columns: 3, left: 54, top: 165, width: 1172, height: 424, columnWidths: [250, 456, 466], values: altRows });
  styleTable(alt, altRows.length, 3);
  shape(slide28, "roundRect", 54, 610, 1172, 49, C.paper, C.violet, 1, "rounded-lg");
  text(slide28, "Future scope is an honest status, not a failure. Each later method must use the same public targets, scoring rules, and protected-data firewall.", 77, 622, 1124, 28, { fontSize: 16, bold: true, alignment: "center" });
  note(slide28, "Presenter 4 (minor)", "Explain that OpenMM minimization addresses atomic strain but does not reproduce Rosetta loop-closure sampling. Other predictor families and linker detectors are useful sensitivity analyses, but none independently validates fusion compatibility. Experimental success has no computational replacement and is outside the requested scope.", ["https://docs.openmm.org/latest/userguide/application/02_running_sims.html", "https://github.com/facebookresearch/esm", "https://github.com/aqlaboratory/openfold", "https://pcrgwd.ucd.ie/linker_pred/download/", "https://github.com/YihePang/DisoFLAG"]);

  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(CANDIDATE);
  const result = await finalizePresentation({
    explicitTotalSlideCount: 28,
    workspaceDir: ROOT,
    candidatePath: CANDIDATE,
    finalPath: input.output,
    pythonExecutable: PYTHON,
    integrityValidatorPath: path.join(SKILL, "container_tools", "inspect_presentation_package_integrity.py"),
    layoutValidatorPath: path.join(SKILL, "container_tools", "inspect_presentation_layout_geometry.py"),
    layoutArgs: [
      "--expected-slide-size-emu", "12192000,6858000", "--validate-bullet-geometry", "--validate-heading-fit",
      "--require-native-table-slide", "6", "--require-native-table-slide", "7", "--require-native-table-slide", "13",
      "--require-native-table-slide", "22", "--require-native-table-slide", "27", "--require-native-table-slide", "28",
    ],
    requiredNativeTableOwnerSlides: [6, 7, 13, 22, 27, 28],
    requiredNativeChartOwnerSlides: [18, 19, 26],
    materializeLiteralChartWorkbooks: true,
    fontPolicy: { basis: "design", families: [FONT, MONO] },
    verifyArtifactToolImport: true,
    receiptPath: path.join(BUILD, "multi-protein-benchmark.final.validation.json"),
  });

  for (const [index, slide] of presentation.slides.items.entries()) {
    const png = await presentation.export({ slide, format: "png", scale: 1 });
    await fs.writeFile(path.join(PREVIEW_DIR, `slide-${String(index + 1).padStart(2, "0")}.png`), new Uint8Array(await png.arrayBuffer()));
  }
  const montage = await presentation.export({ format: "webp", montage: true });
  await fs.writeFile(path.join(BUILD, "multi-protein-benchmark-montage.webp"), new Uint8Array(await montage.arrayBuffer()));
  const progress = JSON.parse(await fs.readFile(input.progress, "utf8"));
  Object.assign(progress, { state: "complete", overall_percent: 100.0, presentation: input.output, presentation_validation: result, updated_utc: new Date().toISOString() });
  await atomicJson(input.progress, progress);
  process.stdout.write(`${JSON.stringify({ output: input.output, slides: 28, result })}\n`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
