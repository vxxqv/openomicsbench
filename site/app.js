const catalog = document.querySelector("#catalog");
const status = document.querySelector("#catalog-status");
const search = document.querySelector("#search");
const filters = [...document.querySelectorAll(".filter")];
const labels = {
  rnaseq: "Bulk RNA-seq",
  sequences: "Sequence",
  singlecell: "Single-cell",
  intervals: "Genomic peaks"
};
let records = [];
let family = "all";

function showRecords() {
  const query = search.value.trim().toLocaleLowerCase();
  const visible = records.filter(item => (family === "all" || item.family === family) && item.haystack.includes(query));
  status.textContent = `${visible.length} of ${records.length} benchmarks`;
  catalog.replaceChildren();
  for (const item of visible) {
    const card = document.createElement("a");
    card.className = "catalog-card";
    card.href = `https://github.com/vxxqv/openomicsbench/tree/main/datasets/${item.family}/${encodeURIComponent(item.id)}`;
    const top = document.createElement("div");
    top.className = "card-top";
    const category = document.createElement("span");
    category.textContent = labels[item.family];
    const id = document.createElement("span");
    id.className = "card-id";
    id.textContent = item.id;
    top.append(category, id);
    const title = document.createElement("h3");
    title.textContent = item.title;
    const detail = document.createElement("p");
    detail.className = "card-detail";
    detail.textContent = item.organism || item.assay.replaceAll("_", " ");
    const bottom = document.createElement("div");
    bottom.className = "card-bottom";
    const rights = document.createElement("span");
    rights.textContent = `Rights ${item.rights}`;
    const link = document.createElement("span");
    link.textContent = "View files ↗";
    bottom.append(rights, link);
    card.append(top, title, detail, bottom);
    catalog.append(card);
  }
}

search.addEventListener("input", showRecords);
for (const button of filters) {
  button.addEventListener("click", () => {
    family = button.dataset.family;
    for (const option of filters) {
      const selected = option === button;
      option.classList.toggle("active", selected);
      option.setAttribute("aria-pressed", String(selected));
    }
    showRecords();
  });
}

fetch("collection.json")
  .then(response => {
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    return response.json();
  })
  .then(data => {
    const groups = [
      ["rnaseq", data.objects],
      ["sequences", data.sequences],
      ["singlecell", data.single_cell],
      ["intervals", data.intervals]
    ];
    records = groups.flatMap(([name, items]) => items.map(item => ({
      ...item,
      family: name,
      haystack: [item.id, item.title, item.organism, item.assay, item.archetype, item.source_accession]
        .filter(Boolean).join(" ").toLocaleLowerCase()
    })));
    showRecords();
  })
  .catch(() => {
    status.textContent = "The catalog is unavailable right now. Browse all datasets on GitHub.";
    const link = document.createElement("a");
    link.href = "https://github.com/vxxqv/openomicsbench/tree/main/datasets";
    link.textContent = "Browse datasets ↗";
    catalog.replaceChildren(link);
  });
