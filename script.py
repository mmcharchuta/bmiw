from Bio import SeqIO

def load_genbank_file(file_path):
    try:
        with open(file_path, 'r') as handle:
            return list(SeqIO.parse(handle, "genbank"))
    except Exception as e:
        print(f"Error loading file {file_path}: {e}")
        return []

def calculate_gene_statistics(records):
    gene_lengths = []
    for record in records:
        for feature in record.features:
            if feature.type == "gene":
                start = int(feature.location.start)
                end = int(feature.location.end)
                gene_lengths.append(abs(end - start))
    gene_count = len(gene_lengths)
    avg_length = sum(gene_lengths) / gene_count if gene_count > 0 else 0
    return gene_count, avg_length

def extract_polymerase_genes(records):
    polymerase_genes = []
    keywords = ["RNA polymerase", "DNA polymerase"]
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and "product" in feature.qualifiers:
                product = feature.qualifiers["product"][0].lower()
                if any(keyword.lower() in product for keyword in keywords):
                    polymerase_genes.append({
                        "location": feature.location,
                        "product": product,
                        "sequence": feature.extract(record.seq)
                    })
    return polymerase_genes

def save_polymerase_genes(polymerase_genes, output_file):
    with open(output_file, 'w') as f:
        for gene in polymerase_genes:
            f.write(f"> {gene['product']} | {gene['location']}\n")
            f.write(f"{gene['sequence']}\n\n")

def main():
    gbff_file = "assembly.gbff"
    gbk_file = "assembly.gbk"

    records_bakta = load_genbank_file(gbff_file)
    records_prokka = load_genbank_file(gbk_file)

    bakta_stats = calculate_gene_statistics(records_bakta)
    prokka_stats = calculate_gene_statistics(records_prokka)

    print("Gene statistics:")
    print(f"Bakta - Number of genes: {bakta_stats[0]}, Average length: {bakta_stats[1]:.2f} bp")
    print(f"Prokka - Number of genes: {prokka_stats[0]}, Average length: {prokka_stats[1]:.2f} bp")

    bakta_polymerases = extract_polymerase_genes(records_bakta)
    prokka_polymerases = extract_polymerase_genes(records_prokka)

    print(f"Number of polymerase genes in Bakta: {len(bakta_polymerases)}")
    print(f"Number of polymerase genes in Prokka: {len(prokka_polymerases)}")

    save_polymerase_genes(bakta_polymerases, "bakta_polymerases.txt")
    save_polymerase_genes(prokka_polymerases, "prokka_polymerases.txt")

    print("Polymerase genes saved to bakta_polymerases.txt and prokka_polymerases.txt")

if __name__ == "__main__":
    main()
