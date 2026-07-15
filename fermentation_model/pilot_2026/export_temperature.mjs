import fs from "node:fs";
import path from "node:path";
import MDBReader from "mdb-reader";

function csvCell(value) {
  if (value === null || value === undefined) return "";
  let text;
  if (value instanceof Date) text = value.toISOString();
  else if (typeof value === "object") text = JSON.stringify(value);
  else text = String(value);
  return /[",\n\r]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

const [inputDir, outputCsv] = process.argv.slice(2);
if (!inputDir || !outputCsv) {
  throw new Error("Usage: node export_temperature.mjs INPUT_DIR OUTPUT_CSV");
}

const files = fs.readdirSync(inputDir)
  .filter((name) => name.toLowerCase().endsWith(".mdb"))
  .sort();

const records = [];
const allColumns = new Set(["source_file", "source_table", "source_row"]);
for (const file of files) {
  const reader = new MDBReader(fs.readFileSync(path.join(inputDir, file)));
  for (const tableName of reader.getTableNames().sort()) {
    const table = reader.getTable(tableName);
    const rows = table.getData({ rowOffset: 0, rowLimit: Number.MAX_SAFE_INTEGER });
    rows.forEach((row, index) => {
      const record = { source_file: file, source_table: tableName, source_row: index + 1, ...row };
      Object.keys(record).forEach((key) => allColumns.add(key));
      records.push(record);
    });
  }
}

const columns = [...allColumns];
const lines = [columns.map(csvCell).join(",")];
for (const record of records) {
  lines.push(columns.map((column) => csvCell(record[column])).join(","));
}
fs.writeFileSync(outputCsv, `${lines.join("\n")}\n`, "utf8");
process.stdout.write(JSON.stringify({ files: files.length, rows: records.length, columns }, null, 2));
