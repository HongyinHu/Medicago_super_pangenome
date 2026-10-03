import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const [inputTsv, outputXlsx, previewPng] = process.argv.slice(2);
if (!inputTsv || !outputXlsx || !previewPng) {
  throw new Error(
    "Usage: node build_supplementary_workbook.mjs <gene_category.tsv> <output.xlsx> <preview.png>",
  );
}

function parseTsv(text) {
  const lines = text.replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/);
  const header = lines[0].split("\t");
  return lines.slice(1).map((line) => {
    const fields = line.split("\t");
    return Object.fromEntries(header.map((name, i) => [name, fields[i] ?? ""]));
  });
}

const rows = parseTsv(await fs.readFile(inputTsv, "utf8"));
if (rows.length !== 18) throw new Error(`Expected 18 genomes, found ${rows.length}`);
for (const r of rows) {
  const fourCategorySum = [
    r.Core_genes_n,
    r.Softcore_genes_n,
    r.Dispensable_genes_n,
    r.Private_genes_n,
  ].reduce((sum, value) => sum + Number(value), 0);
  if (fourCategorySum !== Number(r.Genes_in_orthogroups_n)) {
    throw new Error(`Four-category sum mismatch for ${r.sample}`);
  }
  if (
    Number(r.Genes_in_orthogroups_n) + Number(r.Unassigned_genes_n) !==
    Number(r.Total_genes_n)
  ) {
    throw new Error(`Assigned + unassigned mismatch for ${r.sample}`);
  }
}
console.log("builder: input parsed");

const workbook = Workbook.create();
const main = workbook.worksheets.add("Supplementary Table 11");
const qc = workbook.worksheets.add("QC");
const definitions = workbook.worksheets.add("Definitions");
console.log("builder: workbook created");

main.showGridLines = false;
qc.showGridLines = false;
definitions.showGridLines = false;

main.mergeCells("A1:G1");
main.getRange("A1:G1").values = [[
  "Supplementary Table 11 | Numbers of genes assigned to pangenome prevalence categories in each Medicago genome.",
]];
main.getRange("A3:G3").values = [[
  "Entry",
  "Sample ID",
  "Species",
  "Core genes (n)",
  "Softcore genes (n)",
  "Dispensable genes (n)",
  "Private genes (n)",
]];
main.getRange("A4:G21").values = rows.map((r) => [
  Number(r.Entry),
  r.Sample_ID,
  r.Species,
  Number(r.Core_genes_n),
  Number(r.Softcore_genes_n),
  Number(r.Dispensable_genes_n),
  Number(r.Private_genes_n),
]);

const titleFormat = {
  font: { name: "Times New Roman", size: 12, bold: true, color: "#111111" },
  verticalAlignment: "center",
  horizontalAlignment: "left",
};
main.getRange("A1:G1").format = titleFormat;
main.getRange("A1:G1").format.rowHeight = 28;
main.getRange("A3:G3").format = {
  font: { name: "Times New Roman", size: 10, bold: true, color: "#111111" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
  wrapText: true,
  borders: {
    top: { style: "medium", color: "#111111" },
    bottom: { style: "thin", color: "#111111" },
  },
};
main.getRange("A3:G3").format.rowHeight = 34;
main.getRange("A4:G21").format = {
  font: { name: "Times New Roman", size: 10, color: "#111111" },
  verticalAlignment: "center",
};
main.getRange("A4:B21").format.horizontalAlignment = "center";
main.getRange("D4:G21").format.horizontalAlignment = "right";
main.getRange("D4:G21").format.numberFormat = "#,##0";
main.getRange("C4:C21").format.font = {
  name: "Times New Roman",
  size: 10,
  italic: true,
  color: "#111111",
};
main.getRange("A21:G21").format.borders = {
  bottom: { style: "medium", color: "#111111" },
};
main.getRange("A4:G21").format.rowHeight = 20;
main.getRange("A1:A21").format.columnWidth = 8;
main.getRange("B1:B21").format.columnWidth = 13;
main.getRange("C1:C21").format.columnWidth = 35;
main.getRange("D1:G21").format.columnWidth = 19;
main.freezePanes.freezeRows(3);
console.log("builder: main sheet formatted");

qc.getRange("A1:I1").values = [[
  "Entry",
  "Sample ID",
  "OrthoFinder ID",
  "Genes in orthogroups (n)",
  "Unassigned genes (n)",
  "Total genes (n)",
  "Four-category sum (n)",
  "Assigned difference",
  "Total difference",
]];
qc.getRange("A2:F19").values = rows.map((r) => [
  Number(r.Entry),
  r.Sample_ID,
  r.sample,
  Number(r.Genes_in_orthogroups_n),
  Number(r.Unassigned_genes_n),
  Number(r.Total_genes_n),
]);
console.log("builder: QC values written");
for (let excelRow = 2; excelRow <= 19; excelRow += 1) {
  qc.getRange(`G${excelRow}:I${excelRow}`).formulas = [[
    `=SUM('Supplementary Table 11'!D${excelRow + 2}:G${excelRow + 2})`,
    `=G${excelRow}-D${excelRow}`,
    `=D${excelRow}+E${excelRow}-F${excelRow}`,
  ]];
}
console.log("builder: QC formulas written");
qc.getRange("A1:I1").format = {
  fill: "#365F91",
  font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
  wrapText: true,
};
console.log("builder: QC header formatted");
qc.getRange("A2:I19").format.font = { name: "Arial", size: 9, color: "#111111" };
qc.getRange("A2:C19").format.horizontalAlignment = "left";
qc.getRange("D2:I19").format.horizontalAlignment = "right";
qc.getRange("D2:I19").format.numberFormat = "#,##0";
console.log("builder: QC body formatted");
qc.getRange("A1:A19").format.columnWidth = 8;
qc.getRange("B1:C19").format.columnWidth = 16;
qc.getRange("D1:G19").format.columnWidth = 24;
qc.getRange("H1:I19").format.columnWidth = 18;
qc.getRange("A1:I1").format.rowHeight = 34;
qc.freezePanes.freezeRows(1);
console.log("builder: QC sheet formatted");

definitions.getRange("A1:B8").values = [
  ["Item", "Definition"],
  ["Core", "Orthogroups present in all 18 genomes (frequency = 18)."],
  ["Softcore", "Orthogroups present in 15-17 genomes."],
  ["Dispensable", "Orthogroups present in 2-14 genomes."],
  ["Private", "Orthogroups present in exactly one genome (frequency = 1)."],
  ["Counting unit", "Number of genes contributed by each genome to orthogroups in the indicated family category."],
  ["Unassigned genes", "Reported only on the QC sheet and excluded from Core, Softcore, Dispensable, and Private."],
  ["Source", "OrthoFinder 2.5.4 Orthogroups.GeneCount.tsv, Orthogroups_UnassignedGenes.tsv, and Statistics_PerSpecies.tsv."],
];
definitions.getRange("A1:B1").format = {
  fill: "#365F91",
  font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" },
};
definitions.getRange("A2:B8").format = {
  font: { name: "Arial", size: 9, color: "#111111" },
  verticalAlignment: "top",
  wrapText: true,
};
definitions.getRange("A1:A8").format.columnWidth = 22;
definitions.getRange("B1:B8").format.columnWidth = 78;
definitions.getRange("A2:B8").format.rowHeight = 32;
definitions.freezePanes.freezeRows(1);
console.log("builder: definitions sheet formatted");

const keyInspection = await workbook.inspect({
  kind: "table",
  range: "'Supplementary Table 11'!A1:G21",
  include: "values,formulas",
  tableMaxRows: 25,
  tableMaxCols: 10,
  maxChars: 12000,
});
console.log(keyInspection.ndjson);

const qcInspection = await workbook.inspect({
  kind: "table",
  range: "QC!A1:I19",
  include: "values,formulas",
  tableMaxRows: 22,
  tableMaxCols: 12,
  maxChars: 12000,
});
console.log(qcInspection.ndjson);

const errorScan = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errorScan.ndjson);

await fs.mkdir(path.dirname(outputXlsx), { recursive: true });
console.log("builder: output directory ready");

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
console.log("builder: workbook exported");
await xlsx.save(outputXlsx);
console.log("builder: workbook saved");

if (process.env.SKIP_ARTIFACT_RENDER !== "1") {
  const preview = await workbook.render({
    sheetName: "Supplementary Table 11",
    range: "A1:G21",
    scale: 1.5,
    format: "png",
  });
  console.log("builder: preview rendered");
  await fs.writeFile(previewPng, new Uint8Array(await preview.arrayBuffer()));
  console.log("builder: preview saved");
} else {
  console.log("builder: artifact preview skipped; use Excel-native preview QA");
}

console.log(JSON.stringify({ rows: rows.length, outputXlsx, previewPng }));
